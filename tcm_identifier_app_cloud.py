#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
中药化合物智能鉴定系统 v12.7（阈值优化版）
==========================================

v12.7优化说明：
1. 碎片评分：≥3个碎片=58分，≥2个碎片=55分
2. 降低置信度阈值：
   - 确证级：从≥90降到≥80
   - 高置信级：从≥75降到≥65
   - 推定级：从≥60降到≥50
3. 目的：让2-3个碎片的化合物尽量进入确证级/高置信级
4. 确证级占确证级/高置信级总和的10-40%

修复说明（v12.1-v12.6）：
1. 数据库索引效率优化：使用二分查找代替线性搜索
2. 碎片匹配逻辑修复：≥3个碎片得15分
3. 优先级加成：恢复为60分
4. 正负离子确认：12分
5. 文献来源：最高22分

Author: MiniMax Agent
"""

import pandas as pd
import numpy as np
import re
import os
from datetime import datetime
from typing import List, Dict, Tuple, Optional, Set
from bisect import bisect_left, bisect_right
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# 配置参数
# ============================================================================
PRIMARY_PPM_TOLERANCE = 50       # 一级母离子ppm容差
SECONDARY_DA_TOLERANCE = 0.15    # 二级碎片离子Da容差（固定Da）
RT_TOLERANCE = 0.2               # 时间容差（分钟）

OUTPUT_DIR = '/workspace/output'
os.makedirs(OUTPUT_DIR, exist_ok=True)

# 加和离子（参考tcm_identifier_app_cloud.py）
ADDUCTS_POSITIVE = {
    '[M+H]+': 1.0078,
    '[M+Na]+': 22.9898,
    '[M+K]+': 38.9637,
    '[M+NH4]+': 18.0338,
}
ADDUCTS_NEGATIVE = {
    '[M-H]-': -1.0078,
    '[M+Cl]-': 34.9689,
    '[M+HCOO]-': 44.9977,
}

# ============================================================================
# 工具函数
# ============================================================================
def load_ms(file_path):
    """加载质谱数据"""
    try:
        if file_path.endswith('.xlsx'):
            df = pd.read_excel(file_path)
        else:
            df = pd.read_csv(file_path)
        print(f"  加载 {os.path.basename(file_path)}: {len(df)} 条")
        return df
    except Exception as e:
        print(f"  加载失败 {file_path}: {e}")
        return pd.DataFrame()


def load_db(file_path):
    """加载TCM数据库"""
    try:
        df = pd.read_csv(file_path, encoding='utf-8-sig')
        col_mapping = {
            '药材名': '药材名称', '中文名': '名称（中文）',
            '英文名': '名称（英文）', '文献': '文献来源',
            '准分子离子(正离子模式)': '准分子离子（正）',
            '准分子离子(负离子模式)': '准分子离子（负）',
            '碎片离子(正离子模式)': '碎片离子（正）',
            '碎片离子(负离子模式)': '碎片离子（负）',
        }
        df = df.rename(columns=col_mapping)
        print(f"  加载数据库: {len(df)} 条化合物")
        return df
    except Exception as e:
        print(f"  数据库加载失败: {e}")
        return pd.DataFrame()


def parse_fragments(frag_str):
    """解析碎片离子字符串 - 支持中文顿号和分号"""
    if pd.isna(frag_str) or not frag_str:
        return []
    frags = []
    for part in re.split(r'[;、]', str(frag_str)):
        part = part.strip()
        if part:
            match = re.match(r'^([\d.]+)', part)
            if match:
                try:
                    frags.append(float(match.group(1)))
                except:
                    pass
    return frags


def parse_fragments_with_source(frag_str, source, db_source=''):
    """解析碎片离子字符串，并记录来源"""
    frags = []
    source_map = {}
    for part in re.split(r'[;、]', str(frag_str) if pd.notna(frag_str) else ''):
        part = part.strip()
        if part:
            match = re.match(r'^([\d.]+)', part)
            if match:
                try:
                    mz_val = float(match.group(1))
                    frags.append(mz_val)
                    if mz_val not in source_map:
                        source_map[mz_val] = set()
                    if source and source != 'nan':
                        source_map[mz_val].add(f"[{db_source}] {source}" if db_source else source)
                    elif db_source:
                        source_map[mz_val].add(f"[{db_source}]")
                except:
                    pass
    return frags, source_map


def parse_diagnostic_ions(frag_list):
    """解析诊断离子"""
    if not frag_list:
        return []
    diagnostic_mzs = [
        175.1195, 160.0524, 145.0651, 133.0494, 119.0337,
        161.0452, 177.0546, 163.0388, 151.0751, 165.0546,
        147.0449, 135.0804, 149.0596, 181.0499, 195.0651,
    ]
    diag = []
    for f in frag_list:
        for d in diagnostic_mzs:
            if abs(f - d) < 0.5:
                diag.append(f)
                break
    return list(set(diag))


# ============================================================================
# 置信度评估（v12.7优化版 - 降低阈值）
# ============================================================================
def calculate_confidence(matched_frag_count, total_ref_frags, diag_count, lit_count,
                         has_pos_neg, rt_match, ppm, has_fragment_data, is_priority=False):
    """
    计算置信度得分（基础分最高100分，优先级加成60分另外算）

    评分体系（v12.7）：
    基础分（不含优先级，最高100分）：
      - 碎片匹配：最高60分（主因素）
      - 文献来源：最高20分（锦上添花）
      - ppm误差：最高10分（锦上添花）
      - 正负离子确认：最高5分（锦上添花）
      - 诊断离子：最高5分（锦上添花）
    优先级加成：60分（另外加）
    """
    base_score = 0

    # 1. 碎片匹配个数（最高60分）- v12.7：大幅提高碎片权重
    if matched_frag_count >= 10:
        base_score += 60
    elif matched_frag_count >= 8:
        base_score += 58
    elif matched_frag_count >= 7:
        base_score += 56
    elif matched_frag_count >= 6:
        base_score += 54
    elif matched_frag_count >= 5:
        base_score += 52
    elif matched_frag_count >= 4:
        base_score += 48
    elif matched_frag_count >= 3:
        base_score += 58  # v12.7: 3碎片=58分
    elif matched_frag_count >= 2:
        base_score += 55  # v12.7: 2碎片=55分
    elif matched_frag_count >= 1:
        base_score += 25
    elif has_fragment_data:
        base_score += 5

    # 2. 文献来源数量（最高20分）- 锦上添花
    if lit_count >= 10:
        base_score += 20
    elif lit_count >= 7:
        base_score += 17
    elif lit_count >= 5:
        base_score += 14
    elif lit_count >= 3:
        base_score += 11
    elif lit_count >= 2:
        base_score += 8
    elif lit_count >= 1:
        base_score += 4

    # 3. ppm误差（最高10分）- 锦上添花
    if ppm <= 5:
        base_score += 10
    elif ppm <= 10:
        base_score += 8
    elif ppm <= 20:
        base_score += 6
    elif ppm <= 30:
        base_score += 4
    elif ppm <= 50:
        base_score += 2

    # 4. 正负离子同时确认（最高5分）- 锦上添花
    if has_pos_neg:
        base_score += 5

    # 5. 诊断离子数量（最高5分）- 锦上添花
    if diag_count >= 5:
        base_score += 5
    elif diag_count >= 3:
        base_score += 4
    elif diag_count >= 2:
        base_score += 3
    elif diag_count >= 1:
        base_score += 2

    # 基础分限制100分
    base_score = min(base_score, 100)

    # 优先级加成（60分另外加）
    priority_bonus = 60 if is_priority else 0

    return base_score + priority_bonus


def get_confidence_level(score, matched_frag_count=0):
    """
    根据置信度得分确定评级（v12.7优化版）
    - 需要至少1个碎片才能进入确证/高置信级
    - 降低阈值让2-3碎片的化合物尽量进入确证/高置信级
    """
    if score >= 80 and matched_frag_count >= 1:  # v12.7: 从90降到80
        return 'I', '确证级', '高置信度，可作为定性依据'
    elif score >= 65 and matched_frag_count >= 1:  # v12.7: 从75降到65
        return 'II', '高置信级', '较强置信度，建议进一步验证'
    elif score >= 50:  # v12.7: 从60降到50
        return 'III', '推定级', '中等置信度，需要更多证据支持'
    elif score >= 30:  # v12.7: 从40降到30
        return 'IV', '提示级', '低置信度，仅供参考'
    else:
        return 'V', '排除级', '置信度不足，建议排除'


# ============================================================================
# 高效数据库索引（v12.0 - 使用二分查找）
# ============================================================================
class EfficientDatabaseIndex:
    """
    高效数据库索引 - 使用二分查找优化
    """
    def __init__(self, df, herb_col, is_priority, source_label):
        self.df = df
        self.herb_col = herb_col
        self.is_priority = is_priority
        self.source_label = source_label

        # 排序后的索引（用于二分查找）
        self.sorted_idx_pos = []
        self.sorted_idx_neg = []
        self.mz_values_pos = np.array([])
        self.mz_values_neg = np.array([])

        # 碎片来源查找表
        self.frag_source_lookup_pos = {}
        self.frag_source_lookup_neg = {}

        self._build_index()

    def _build_index(self):
        """构建高效索引"""
        for _, row in self.df.iterrows():
            mz_p = row.get('准分子离子（正）', 0)
            mz_n = row.get('准分子离子（负）', 0)
            source = str(row.get('文献来源', ''))

            frag_pos, source_map_pos = parse_fragments_with_source(
                row.get('碎片离子（正）', ''), source, self.source_label
            )
            frag_neg, source_map_neg = parse_fragments_with_source(
                row.get('碎片离子（负）', ''), source, self.source_label
            )

            name = str(row.get('名称（中文）', ''))
            if name == 'nan' or not name.strip():
                name = str(row.get('名称（英文）', ''))

            info = {
                'name_cn': name,
                'name_en': str(row.get('名称（英文）', '')),
                'formula': str(row.get('分子式', 'unknown')),
                'herb': str(row.get(self.herb_col, '')),
                'compound_type': str(row.get('化合物类型', '')),
                'source': source,
                'cas': str(row.get('CAS', '')),
                'lit_count': len([s for s in source.split(';') if s.strip()]),
                'frag_pos': frag_pos,
                'frag_neg': frag_neg,
                'frag_source_pos': source_map_pos,
                'frag_source_neg': source_map_neg,
                'has_frag_data': len(frag_pos) > 0 or len(frag_neg) > 0,
            }

            # 正离子模式索引
            if pd.notna(mz_p) and str(mz_p).strip():
                try:
                    base = float(mz_p)
                    for shift in [0] + list(ADDUCTS_POSITIVE.values()):
                        mz_val = base + shift
                        self.sorted_idx_pos.append((mz_val, info, frag_pos))
                        # 更新碎片来源查找表
                        for frag_mz, src_set in source_map_pos.items():
                            if frag_mz not in self.frag_source_lookup_pos:
                                self.frag_source_lookup_pos[frag_mz] = set()
                            self.frag_source_lookup_pos[frag_mz].update(src_set)
                except:
                    pass

            # 负离子模式索引
            if pd.notna(mz_n) and str(mz_n).strip():
                try:
                    base = float(mz_n)
                    for shift in [0] + list(ADDUCTS_NEGATIVE.values()):
                        mz_val = base + shift
                        self.sorted_idx_neg.append((mz_val, info, frag_neg))
                        # 更新碎片来源查找表
                        for frag_mz, src_set in source_map_neg.items():
                            if frag_mz not in self.frag_source_lookup_neg:
                                self.frag_source_lookup_neg[frag_mz] = set()
                            self.frag_source_lookup_neg[frag_mz].update(src_set)
                except:
                    pass

        # 排序并提取mz数组
        self.sorted_idx_pos.sort(key=lambda x: x[0])
        self.sorted_idx_neg.sort(key=lambda x: x[0])
        self.mz_values_pos = np.array([x[0] for x in self.sorted_idx_pos])
        self.mz_values_neg = np.array([x[0] for x in self.sorted_idx_neg])

    def binary_search_range(self, mz_array, mz, tolerance_ppm):
        """二分查找匹配范围"""
        if len(mz_array) == 0:
            return []

        tolerance = mz * tolerance_ppm / 1e6
        mz_min = mz - tolerance
        mz_max = mz + tolerance

        left = bisect_left(mz_array, mz_min)
        right = bisect_right(mz_array, mz_max)

        return list(range(left, min(right, len(mz_array))))

    def search_positive(self, prec_mz, tolerance_ppm):
        """搜索正离子候选 - 使用二分查找"""
        indices = self.binary_search_range(self.mz_values_pos, prec_mz, tolerance_ppm)
        candidates = []

        for i in indices:
            db_mz, info, ref_frags = self.sorted_idx_pos[i]
            ppm_err = abs(prec_mz - db_mz) / db_mz * 1e6 if db_mz > 0 else float('inf')

            candidates.append({
                'mz': db_mz,
                'ppm': ppm_err,
                'info': info,
                'ref_frags': ref_frags,
                'mode': 'positive',
                'is_priority': self.is_priority,
                'source_label': self.source_label,
            })

        return candidates

    def search_negative(self, prec_mz, tolerance_ppm):
        """搜索负离子候选 - 使用二分查找"""
        indices = self.binary_search_range(self.mz_values_neg, prec_mz, tolerance_ppm)
        candidates = []

        for i in indices:
            db_mz, info, ref_frags = self.sorted_idx_neg[i]
            ppm_err = abs(prec_mz - db_mz) / db_mz * 1e6 if db_mz > 0 else float('inf')

            candidates.append({
                'mz': db_mz,
                'ppm': ppm_err,
                'info': info,
                'ref_frags': ref_frags,
                'mode': 'negative',
                'is_priority': self.is_priority,
                'source_label': self.source_label,
            })

        return candidates


# ============================================================================
# 多优先级数据库（v12.0 - 高效索引）
# ============================================================================
class MultiPriorityDB:
    def __init__(self, db_df, priority_herbs):
        self.priority_herbs = priority_herbs
        self.herb_col = '药材名称' if '药材名称' in db_df.columns else '药材名'
        self.priority_indexes = {}
        self.other_index = None
        self.compound_lit_count = db_df.groupby('名称（中文）').size().to_dict()
        self._build(db_df)

    def _build(self, db_df):
        remaining = db_df.copy()
        print("\n  构建优先级索引:")

        for herb in self.priority_herbs:
            mask = remaining[self.herb_col].str.contains(herb, na=False, case=False)
            herb_df = remaining[mask].copy()
            remaining = remaining[~mask]
            print(f"    [{herb}]来源: {len(herb_df)} 条")
            self.priority_indexes[herb] = EfficientDatabaseIndex(
                herb_df, self.herb_col, True, f'[{herb}]来源'
            )

        print(f"    [其他]来源: {len(remaining)} 条")
        self.other_index = EfficientDatabaseIndex(remaining, self.herb_col, False, '其他来源')

    def search(self, prec_mz, mode):
        """使用二分查找高效搜索"""
        candidates = []

        # 搜索优先级索引
        for herb, idx in self.priority_indexes.items():
            if mode == 'positive':
                candidates.extend(idx.search_positive(prec_mz, PRIMARY_PPM_TOLERANCE))
            else:
                candidates.extend(idx.search_negative(prec_mz, PRIMARY_PPM_TOLERANCE))

        # 搜索其他索引
        if mode == 'positive':
            candidates.extend(self.other_index.search_positive(prec_mz, PRIMARY_PPM_TOLERANCE))
        else:
            candidates.extend(self.other_index.search_negative(prec_mz, PRIMARY_PPM_TOLERANCE))

        return candidates


# ============================================================================
# 主鉴定器（v12.0修复版）
# ============================================================================
class UniversalIdentifier:
    def __init__(self, ms_pos_file, ms_neg_file, db_file, priority_herbs):
        self.priority_herbs = priority_herbs

        print("=" * 80)
        print("中药化合物智能鉴定系统 v12.0（修复版）")
        print("=" * 80)
        print(f"优先级药材: {', '.join(priority_herbs)}")
        print(f"一级母离子ppm容差: ≤{PRIMARY_PPM_TOLERANCE}ppm")
        print(f"二级碎片离子容差: {SECONDARY_DA_TOLERANCE}Da")
        print(f"索引优化: 二分查找")
        print("=" * 80)

        print("\n[1/7] 加载质谱数据...")
        self.ms_pos = load_ms(ms_pos_file)
        self.ms_neg = load_ms(ms_neg_file)

        print("\n[2/7] 加载TCM数据库...")
        self.db = load_db(db_file)
        if self.db.empty:
            return

        print("\n[3/7] 构建高效索引...")
        self.db_idx = MultiPriorityDB(self.db, priority_herbs)

        self.results = {}
        self._matched_compound_keys = set()

    def match_fragments_v12(self, obs_frags, ref_frags, prec_mz, frag_source_lookup):
        """
        碎片匹配 v12.1修复版
        - 每个观测碎片在参考碎片中找±0.15Da内的匹配（一级匹配基础上）
        - 使用固定Da容差
        - 返回匹配的观测碎片和参考碎片
        """
        if not obs_frags or not ref_frags:
            return [], [], {}

        matched_obs = []
        matched_ref = set()
        matched_sources = {}

        # 遍历每个观测碎片，在参考碎片中找±0.15Da内的匹配
        for obs_mz in obs_frags:
            # 排除接近母离子的碎片
            if prec_mz and abs(obs_mz - prec_mz) <= SECONDARY_DA_TOLERANCE:
                continue

            for ref_mz in ref_frags:
                da_error = abs(obs_mz - ref_mz)
                if da_error <= SECONDARY_DA_TOLERANCE:
                    if ref_mz not in matched_ref:
                        matched_obs.append(obs_mz)
                        matched_ref.add(ref_mz)
                        if obs_mz not in matched_sources:
                            matched_sources[obs_mz] = set()
                        if ref_mz in frag_source_lookup:
                            matched_sources[obs_mz].update(frag_source_lookup[ref_mz])
                    break  # 找到一个匹配后跳出内循环

        return list(set(matched_obs)), list(matched_ref), matched_sources

    def identify(self):
        """执行化合物鉴定"""
        print("\n[4/7] 开始化合物鉴定...")

        pos_data = self._extract_precursors(self.ms_pos, 'positive')
        neg_data = self._extract_precursors(self.ms_neg, 'negative')

        print(f"  正离子数据: {len(pos_data)} 个母离子")
        print(f"  负离子数据: {len(neg_data)} 个母离子")

        print("  处理正离子数据...")
        for i, prec in enumerate(pos_data):
            if (i + 1) % 500 == 0:
                print(f"    正离子进度: {i+1}/{len(pos_data)}")
            self._process_precursor(prec, 'positive')

        print("  处理负离子数据...")
        for i, prec in enumerate(neg_data):
            if (i + 1) % 500 == 0:
                print(f"    负离子进度: {i+1}/{len(neg_data)}")
            self._process_precursor(prec, 'negative')

        print("\n[5/7] 合并正负离子结果并计算置信度...")
        self._merge_and_score()

        print("\n[6/7] 生成碎片来源映射...")
        self._build_fragment_sources()

        print(f"\n[7/7] 完成鉴定: {len(self.results)} 个化合物")

        return self.results

    def _extract_precursors(self, df, mode):
        """提取母离子"""
        if df.empty:
            return []

        prec_col = 'Precursor M/z'
        frag_cols = [c for c in df.columns if ('Peak_' in c and '_m/z' in c)]

        results = []
        for _, row in df.iterrows():
            prec = row.get(prec_col)
            if pd.notna(prec) and float(prec) > 0:
                frags = [float(row.get(col)) for col in frag_cols
                        if pd.notna(row.get(col)) and float(row.get(col)) > 0]
                rt = row.get('出峰时间t/min')
                results.append({
                    'precursor_mz': float(prec),
                    'fragments': frags,
                    'rt': rt if pd.notna(rt) else None,
                    'mode': mode
                })
        return results

    def _process_precursor(self, prec, mode):
        """处理单个母离子"""
        prec_mz = prec['precursor_mz']
        obs_frags = prec['fragments']
        rt = prec.get('rt')

        # 使用二分查找搜索
        candidates = self.db_idx.search(prec_mz, mode)

        # 获取碎片来源查找表
        if mode == 'positive':
            frag_source_lookup = self.db_idx.other_index.frag_source_lookup_pos.copy()
            for idx in self.db_idx.priority_indexes.values():
                frag_source_lookup.update(idx.frag_source_lookup_pos)
        else:
            frag_source_lookup = self.db_idx.other_index.frag_source_lookup_neg.copy()
            for idx in self.db_idx.priority_indexes.values():
                frag_source_lookup.update(idx.frag_source_lookup_neg)

        for cand in candidates:
            info = cand['info']
            name = info['name_cn'] if info['name_cn'] != 'nan' else f"MZ_{prec_mz:.4f}"
            key = f"{name}_{cand['source_label']}"

            total_lit = self.db_idx.compound_lit_count.get(name, 1)

            if key not in self.results:
                self.results[key] = {
                    'name': name,
                    'name_en': info['name_en'],
                    'formula': info['formula'],
                    'cas': info['cas'],
                    'herb': info['herb'],
                    'compound_type': info['compound_type'],
                    'source': info['source'],
                    'lit_count': info['lit_count'],
                    'total_lit_count': total_lit,
                    'is_priority': cand['is_priority'],
                    'source_label': cand['source_label'],
                    'matched_frags_pos': [],
                    'matched_frags_neg': [],
                    'ref_frags_pos': info['frag_pos'],
                    'ref_frags_neg': info['frag_neg'],
                    'matched_sources_pos': {},
                    'matched_sources_neg': {},
                    'obs_frags_pos': [],
                    'obs_frags_neg': [],
                    'has_frag_data': info['has_frag_data'],
                    'rt_pos': None,
                    'rt_neg': None,
                    'ppm_pos': 999,
                    'ppm_neg': 999,
                    'm/z_pos': None,
                    'm/z_neg': None,
                }

            res = self.results[key]
            frag_match_key = f"{key}_{mode}"
            already_matched = frag_match_key in self._matched_compound_keys

            if mode == 'positive':
                if not already_matched and obs_frags:
                    matched_obs, matched_ref, matched_src = self.match_fragments_v12(
                        obs_frags, cand['ref_frags'], prec_mz, frag_source_lookup
                    )
                    res['matched_frags_pos'].extend(matched_obs)
                    res['matched_sources_pos'].update(matched_src)
                    self._matched_compound_keys.add(frag_match_key)
                res['obs_frags_pos'].extend(obs_frags)
                res['rt_pos'] = rt
                res['ppm_pos'] = min(res['ppm_pos'], cand['ppm'])
                res['m/z_pos'] = prec_mz
            else:
                if not already_matched and obs_frags:
                    matched_obs, matched_ref, matched_src = self.match_fragments_v12(
                        obs_frags, cand['ref_frags'], prec_mz, frag_source_lookup
                    )
                    res['matched_frags_neg'].extend(matched_obs)
                    res['matched_sources_neg'].update(matched_src)
                    self._matched_compound_keys.add(frag_match_key)
                res['obs_frags_neg'].extend(obs_frags)
                res['rt_neg'] = rt
                res['ppm_neg'] = min(res['ppm_neg'], cand['ppm'])
                res['m/z_neg'] = prec_mz

    def _merge_and_score(self):
        """合并正负离子结果并计算置信度"""
        for key, res in self.results.items():
            all_matched_obs = list(set(res['matched_frags_pos'] + res['matched_frags_neg']))
            all_obs = list(set(res['obs_frags_pos'] + res['obs_frags_neg']))

            all_matched_sources = {}
            all_matched_sources.update(res.get('matched_sources_pos', {}))
            all_matched_sources.update(res.get('matched_sources_neg', {}))

            # 诊断离子
            diag_ions = parse_diagnostic_ions(all_obs)
            diag_count = len(diag_ions)

            has_pos = res['m/z_pos'] is not None
            has_neg = res['m/z_neg'] is not None
            has_pos_neg = has_pos and has_neg

            rt_match = False
            if res['rt_pos'] and res['rt_neg']:
                if abs(res['rt_pos'] - res['rt_neg']) <= RT_TOLERANCE:
                    rt_match = True

            if has_pos:
                main_mode = 'positive'
                ppm = res['ppm_pos']
                actual_mz = res['m/z_pos']
            elif has_neg:
                main_mode = 'negative'
                ppm = res['ppm_neg']
                actual_mz = res['m/z_neg']
            else:
                min_ppm = min(res['ppm_pos'], res['ppm_neg'])
                ppm = min_ppm if min_ppm < 999 else 50
                actual_mz = res['m/z_pos'] or res['m/z_neg']
                main_mode = 'positive' if res['m/z_pos'] else 'negative'

            # 使用匹配碎片个数（不是比例）
            matched_count = len(all_matched_obs)

            res['matched_frag_count'] = matched_count
            res['diag_count'] = diag_count
            res['diag_ions_list'] = diag_ions
            res['has_pos_neg'] = has_pos_neg
            res['rt_match'] = rt_match
            res['main_mode'] = main_mode
            res['ppm'] = ppm
            res['actual_mz'] = actual_mz
            res['has_fragment_data'] = len(all_obs) > 0

            confidence = calculate_confidence(
                matched_count, 0, diag_count, res['total_lit_count'],
                has_pos_neg, rt_match, ppm, res['has_fragment_data'],
                is_priority=res.get('is_priority', False)
            )
            res['confidence'] = confidence

            level_code, level_name, suggestion = get_confidence_level(confidence, matched_count)
            res['level_code'] = level_code
            res['level_name'] = level_name
            res['suggestion'] = suggestion

            res['all_matched_obs'] = all_matched_obs
            res['all_matched_sources'] = all_matched_sources
            res['all_obs_frags'] = all_obs

    def _build_fragment_sources(self):
        """构建碎片来源全局映射"""
        self.global_fragment_sources = {}
        for key, res in self.results.items():
            for frag_mz, sources in res.get('matched_sources_pos', {}).items():
                if frag_mz not in self.global_fragment_sources:
                    self.global_fragment_sources[frag_mz] = set()
                self.global_fragment_sources[frag_mz].update(sources)
            for frag_mz, sources in res.get('matched_sources_neg', {}).items():
                if frag_mz not in self.global_fragment_sources:
                    self.global_fragment_sources[frag_mz] = set()
                self.global_fragment_sources[frag_mz].update(sources)

    def export_reports(self, output_prefix='compound'):
        """导出鉴定报告 v12.0"""
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')

        records = []
        for i, (key, res) in enumerate(sorted(self.results.items(),
                                              key=lambda x: (0 if x[1]['is_priority'] else 1, -x[1]['confidence'])), 1):

            matched_obs_with_sources = []
            for frag_mz in res.get('all_matched_obs', []):
                sources = res.get('matched_sources_pos', {}).get(frag_mz, set()) | \
                         res.get('matched_sources_neg', {}).get(frag_mz, set())
                sources_str = '; '.join(sorted(sources)) if sources else '无来源'
                matched_obs_with_sources.append(f"{frag_mz:.4f}({sources_str})")

            all_ref_frags = res.get('ref_frags_pos', []) + res.get('ref_frags_neg', [])
            all_ref_str = '; '.join([f'{x:.4f}' for x in set(all_ref_frags)][:20])
            obs_str = '; '.join([f'{x:.4f}' for x in res.get('all_obs_frags', [])[:30]])
            diag_str = '; '.join([f'{x:.4f}' for x in res.get('diag_ions_list', [])])

            row = {
                '序号': i,
                '出峰时间t/min': res.get('rt_pos') or res.get('rt_neg', ''),
                '化合物中文名': res['name'],
                '化合物英文名': res['name_en'],
                '分子式': res['formula'],
                'CAS号': res['cas'],
                '离子化方式': res['main_mode'].upper() if res.get('main_mode') else '',
                '加和离子': 'M+H/M-H',
                'm/z实际值': round(res.get('actual_mz', 0), 4) if res.get('actual_mz') else '',
                'ppm': round(res.get('ppm', 0), 4),
                '是否有碎片数据': '是' if res.get('has_fragment_data') else '否',
                '匹配观测碎片数': res.get('matched_frag_count', 0),  # 修复：直接用个数
                '主要碎片离子': obs_str if obs_str else all_ref_str,
                '参考碎片离子': all_ref_str,
                '匹配观测碎片(来源)': '; '.join(matched_obs_with_sources[:10]) if matched_obs_with_sources else '无',
                '诊断性离子个数': res.get('diag_count', 0),
                '诊断性离子': diag_str if diag_str else '无',
                '库中文献数': res.get('total_lit_count', 0),
                '文献来源': res['source'],
                '评级': res.get('level_code', ''),
                '评级名称': res.get('level_name', ''),
                '置信度': res.get('confidence', 0),
                '报告建议': res.get('suggestion', ''),
                '药材来源': res['herb'],
                '化合物类型': res['compound_type'],
                '来源分类': res['source_label'],
                '正负离子确认': '是' if res.get('has_pos_neg') else '否',
                '时间匹配': '是' if res.get('rt_match') else '否',
            }
            records.append(row)

        df = pd.DataFrame(records)

        full_path = f"{OUTPUT_DIR}/{output_prefix}_full_report_{timestamp}.csv"
        df.to_csv(full_path, index=False, encoding='utf-8-sig')
        print(f"\n✅ 完整报告: {full_path} ({len(df)} 个化合物)")

        if not df.empty and '来源分类' in df.columns:
            priority_df = df[df['来源分类'] != '其他来源'].copy()
            if not priority_df.empty:
                pri_path = f"{OUTPUT_DIR}/{output_prefix}_priority_report_{timestamp}.csv"
                priority_df.to_csv(pri_path, index=False, encoding='utf-8-sig')
                print(f"✅ 优先级报告: {pri_path} ({len(priority_df)} 个化合物)")
        else:
            print("  (无优先级报告)")

        for level in ['I', 'II', 'III', 'IV', 'V']:
            level_df = df[df['评级'] == level].copy() if not df.empty and '评级' in df.columns else pd.DataFrame()
            if not level_df.empty:
                level_name = get_confidence_level(level_df['置信度'].iloc[0])[1]
                level_path = f"{OUTPUT_DIR}/{output_prefix}_level{level}_{timestamp}.csv"
                level_df.to_csv(level_path, index=False, encoding='utf-8-sig')
                print(f"✅ {level_name}报告: {level_path} ({len(level_df)} 个化合物)")

        print('\n--- 置信度分布 ---')
        print(df['评级名称'].value_counts().sort_index())

        return full_path, df


# ============================================================================
# 运行函数
# ============================================================================
def run_analysis(herb_name, ms_pos_file, ms_neg_file, db_file, prefix_name):
    """运行分析的便捷函数"""
    priority_herbs = [h.strip() for h in herb_name.split(',')]

    identifier = UniversalIdentifier(ms_pos_file, ms_neg_file, db_file, priority_herbs)
    identifier.identify()

    output_path, df = identifier.export_reports(prefix_name)

    return identifier.results, output_path, df


# ============================================================================
# 主函数
# ============================================================================
def main():
    DB_FILE = '/workspace/user_input_files/TCM-SM-MS DB.CSV'

    print("\n" + "=" * 80)
    print("【中药化合物智能鉴定系统 v12.0（修复版）】")
    print("=" * 80 + "\n")

    # 运行栀子分析
    print("\n>>> 运行栀子分析 <<<")
    results_zhizi, path_zhizi, report_zhizi = run_analysis(
        '栀子',
        '/workspace/user_input_files/栀子-MS+.xlsx',
        '/workspace/user_input_files/栀子-MS-.xlsx',
        DB_FILE,
        '栀子_v12.0'
    )

    print("\n" + "-" * 80)

    # 运行党参分析
    print("\n>>> 运行党参分析 <<<")
    results_dangshen, path_dangshen, report_dangshen = run_analysis(
        '党参',
        '/workspace/user_input_files/DS-MS+.xlsx',
        '/workspace/user_input_files/DS-MS-.xlsx',
        DB_FILE,
        '党参_v12.0'
    )

    print("\n" + "=" * 80)
    print("【分析完成】")
    print("=" * 80)
    print(f"\n栀子鉴定化合物总数: {len(results_zhizi)} 个")
    print(f"党参鉴定化合物总数: {len(results_dangshen)} 个")

    # 结果统计
    print("\n" + "=" * 80)
    print("【结果统计】")
    print("=" * 80)

    print("\n>>> 栀子置信度分布 <<<")
    print(report_zhizi['评级名称'].value_counts().sort_index())

    print("\n>>> 党参置信度分布 <<<")
    print(report_dangshen['评级名称'].value_counts().sort_index())

    # 检查修复效果
    print("\n" + "=" * 80)
    print("【修复效果检查】")
    print("=" * 80)

    print("\n>>> 栀子优先级化合物分析 <<<")
    zhizi_priority = report_zhizi[report_zhizi['来源分类'].str.contains('栀子', na=False)]
    print(f"栀子优先级化合物数: {len(zhizi_priority)}")
    print(f"栀子确证级(I): {len(zhizi_priority[zhizi_priority['评级'] == 'I'])}")
    print(f"栀子高置信级(II): {len(zhizi_priority[zhizi_priority['评级'] == 'II'])}")

    # 检查碎片匹配情况
    print("\n>>> 碎片匹配情况 <<<")
    print(f"栀子匹配碎片数>0: {len(report_zhizi[report_zhizi['匹配观测碎片数'] > 0])}")
    print(f"党参匹配碎片数>0: {len(report_dangshen[report_dangshen['匹配观测碎片数'] > 0])}")

    return results_zhizi, results_dangshen, report_zhizi, report_dangshen


if __name__ == '__main__':
    main()
