#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
中药化合物智能鉴定系统 v11.0（融合版）
=============================================

核心逻辑：universal_identifier_v94.py（多优先级药材分析）
界面设计：tcm_identifier_app_cloud.py（现代化Web界面）

功能特点：
- 多优先级药材同时分析
- 一级母离子匹配：ppm ≤ 50
- 二级碎片离子匹配：100Da容差
- 置信度分级：确证级 > 高置信级 > 推定级 > 提示级 > 排除级
- 优先比对目标药材来源的化合物，再比对数据库其他所有

修复内容：
1. 碎片匹配：用库参考碎片匹配实验碎片（而非相反）
2. 匹配计数：正确累加匹配的碎片数
3. 文献统计：按化合物聚合行数（库中文献数）
4. 二级匹配：一级母离子ppm匹配通过后再进行碎片匹配
5. 置信度：按碎片匹配比例给分

Author: MiniMax Agent
"""

import pandas as pd
import numpy as np
import re
import os
from datetime import datetime
from typing import List, Dict, Tuple, Optional
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# 配置参数
# ============================================================================
PRIMARY_PPM_TOLERANCE = 50       # 一级母离子ppm容差
SECONDARY_DA_TOLERANCE = 100    # 二级碎片离子Da容差
RT_TOLERANCE = 0.2              # 时间容差（分钟）

OUTPUT_DIR = './output'
os.makedirs(OUTPUT_DIR, exist_ok=True)

# 加和离子
ADDUCTS_POSITIVE = {
    '[M+H]+': 1.0073, '[M+Na]+': 22.9898, '[M+K]+': 38.9637,
    '[M+NH4]+': 18.0344, '[M+H-H2O]+': -17.0027,
}
ADDUCTS_NEGATIVE = {
    '[M-H]-': -1.0073, '[M+Na-2H]-': 20.9747,
    '[M+K-2H]-': 36.9556, '[M-H+H2O]-': 18.0106,
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
        df = pd.read_csv(file_path, encoding='utf-8')
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


def parse_diagnostic_ions(frag_list):
    """解析诊断离子"""
    if not frag_list:
        return []
    diagnostic_mzs = [
        175.1195, 160.0524, 145.0651, 133.0494, 119.0337,
        161.0452, 177.0546, 163.0388, 151.0751, 165.0546,
        147.0449, 135.0804, 149.0596, 181.0499, 195.0651,
        209.0807, 223.0964, 237.112, 251.1276, 265.1433,
        129.0541, 113.0234, 97.0287, 85.0287, 71.0371
    ]
    diag = []
    for f in frag_list:
        for d in diagnostic_mzs:
            if abs(f - d) < 0.5:
                diag.append(f)
                break
    return list(set(diag))


# ============================================================================
# 置信度评估
# ============================================================================
def calculate_confidence(matched_frag_count, total_ref_frags, diag_count, lit_count,
                         has_pos_neg, rt_match, ppm, has_fragment_data):
    """
    计算置信度得分（0-100）

    评分维度：
    1. 文献来源数量（最高25分）
    2. 正负离子同时确认（最高25分）
    3. ppm误差（最高20分）
    4. 碎片离子匹配比例（最高15分）- 按比例给分
    5. 时间点匹配（最高10分）
    6. 诊断离子数量（最高5分）
    """
    score = 0

    # 1. 文献来源数量（0-25分）
    if lit_count >= 10:
        score += 25
    elif lit_count >= 7:
        score += 22
    elif lit_count >= 5:
        score += 18
    elif lit_count >= 3:
        score += 15
    elif lit_count >= 2:
        score += 10
    elif lit_count >= 1:
        score += 5

    # 2. 正负离子同时确认（0-25分）
    if has_pos_neg:
        score += 25
    else:
        score += 8

    # 3. ppm误差（0-20分）
    if ppm <= 5:
        score += 20
    elif ppm <= 10:
        score += 18
    elif ppm <= 20:
        score += 15
    elif ppm <= 30:
        score += 12
    elif ppm <= 40:
        score += 8
    elif ppm <= 50:
        score += 5

    # 4. 碎片离子匹配比例（0-15分）- 按比例给分
    if total_ref_frags > 0:
        match_ratio = min(matched_frag_count / total_ref_frags, 1.0)
        frag_score = match_ratio * 15
        score += frag_score
    elif has_fragment_data:
        score += 2

    # 5. 时间点匹配（0-10分）
    if rt_match:
        score += 10

    # 6. 诊断离子数量（0-5分）
    if diag_count >= 5:
        score += 5
    elif diag_count >= 3:
        score += 4
    elif diag_count >= 2:
        score += 3
    elif diag_count >= 1:
        score += 2

    return min(score, 100)


def get_confidence_level(score):
    """根据置信度得分确定评级"""
    if score >= 90:
        return 'I', '确证级', '高置信度，可作为定性依据'
    elif score >= 75:
        return 'II', '高置信级', '较强置信度，建议进一步验证'
    elif score >= 60:
        return 'III', '推定级', '中等置信度，需要更多证据支持'
    elif score >= 40:
        return 'IV', '提示级', '低置信度，仅供参考'
    else:
        return 'V', '排除级', '置信度不足，建议排除'


# ============================================================================
# 数据库索引
# ============================================================================
class DatabaseIndex:
    def __init__(self, df, herb_col, is_priority, source_label):
        self.df = df
        self.herb_col = herb_col
        self.is_priority = is_priority
        self.source_label = source_label
        self.records = []
        self._build_index()

    def _build_index(self):
        for _, row in self.df.iterrows():
            mz_p = row.get('准分子离子（正）', 0)
            mz_n = row.get('准分子离子（负）', 0)
            source = str(row.get('文献来源', ''))
            lit_count = len([s for s in source.split(';') if s.strip()])

            # 解析碎片离子 - 支持中文顿号
            frag_pos = parse_fragments(row.get('碎片离子（正）', ''))
            frag_neg = parse_fragments(row.get('碎片离子（负）', ''))
            diag_pos = parse_diagnostic_ions(frag_pos)
            diag_neg = parse_diagnostic_ions(frag_neg)

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
                'lit_count': lit_count,
                'ref_frag_count_pos': len(frag_pos),
                'ref_frag_count_neg': len(frag_neg),
                'frag_pos': frag_pos,
                'frag_neg': frag_neg,
                'diag_pos': diag_pos,
                'diag_neg': diag_neg,
                'has_frag_data': len(frag_pos) > 0 or len(frag_neg) > 0,
            }

            # 正离子模式索引
            if pd.notna(mz_p) and str(mz_p).strip():
                try:
                    base = float(mz_p)
                    for shift in [0] + list(ADDUCTS_POSITIVE.values()):
                        self.records.append({
                            'mz': base + shift,
                            'mode': 'positive',
                            'info': info,
                            'ref_frags': frag_pos,
                            'ref_diag': diag_pos,
                            'source_label': self.source_label,
                        })
                except:
                    pass

            # 负离子模式索引
            if pd.notna(mz_n) and str(mz_n).strip():
                try:
                    base = float(mz_n)
                    for shift in [0] + list(ADDUCTS_NEGATIVE.values()):
                        self.records.append({
                            'mz': base + shift,
                            'mode': 'negative',
                            'info': info,
                            'ref_frags': frag_neg,
                            'ref_diag': diag_neg,
                            'source_label': self.source_label,
                        })
                except:
                    pass


# ============================================================================
# 多优先级数据库
# ============================================================================
class MultiPriorityDB:
    def __init__(self, db_df, priority_herbs):
        self.priority_herbs = priority_herbs
        self.herb_col = '药材名称' if '药材名称' in db_df.columns else '药材名'
        self.priority_indexes = {}
        self.other_index = None
        # 统计每个化合物的文献数（按名称聚合）
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
            self.priority_indexes[herb] = DatabaseIndex(herb_df, self.herb_col, True, f'[{herb}]来源')

        print(f"    [其他]来源: {len(remaining)} 条")
        self.other_index = DatabaseIndex(remaining, self.herb_col, False, '其他来源')

    def search(self, prec_mz, mode):
        """搜索匹配的化合物"""
        candidates = []

        # 搜索优先级索引
        for herb, idx in self.priority_indexes.items():
            for rec in idx.records:
                if rec['mode'] != mode:
                    continue
                ppm_err = abs(prec_mz - rec['mz']) / prec_mz * 1e6
                if ppm_err <= PRIMARY_PPM_TOLERANCE:
                    candidates.append({
                        'ppm': ppm_err,
                        'theoretical_mz': rec['mz'],
                        'info': rec['info'],
                        'ref_frags': rec['ref_frags'],
                        'ref_diag': rec['ref_diag'],
                        'mode': mode,
                        'is_priority': True,
                        'source_label': rec['source_label'],
                    })

        # 搜索其他索引
        for rec in self.other_index.records:
            if rec['mode'] != mode:
                continue
            ppm_err = abs(prec_mz - rec['mz']) / prec_mz * 1e6
            if ppm_err <= PRIMARY_PPM_TOLERANCE:
                candidates.append({
                    'ppm': ppm_err,
                    'theoretical_mz': rec['mz'],
                    'info': rec['info'],
                    'ref_frags': rec['ref_frags'],
                    'ref_diag': rec['ref_diag'],
                    'mode': mode,
                    'is_priority': False,
                    'source_label': rec['source_label'],
                })

        return candidates


# ============================================================================
# 主鉴定器（修复版）
# ============================================================================
class UniversalIdentifier:
    def __init__(self, ms_pos_file, ms_neg_file, db_file, priority_herbs):
        self.priority_herbs = priority_herbs

        print("=" * 80)
        print("中药化合物智能鉴定系统 v11.0（融合版）")
        print("=" * 80)
        print(f"优先级药材: {', '.join(priority_herbs)}")
        print(f"一级母离子ppm容差: ≤{PRIMARY_PPM_TOLERANCE}ppm")
        print(f"二级碎片离子容差: {SECONDARY_DA_TOLERANCE}Da")
        print("置信度分级: 确证级 > 高置信级 > 推定级 > 提示级 > 排除级")
        print("=" * 80)

        print("\n[1/5] 加载质谱数据...")
        self.ms_pos = load_ms(ms_pos_file)
        self.ms_neg = load_ms(ms_neg_file)

        print("\n[2/5] 加载TCM数据库...")
        self.db = load_db(db_file)
        if self.db.empty:
            return

        print("\n[3/5] 构建优先级索引...")
        self.db_idx = MultiPriorityDB(self.db, priority_herbs)

        self.results = {}

    def match_fragments(self, obs_frags, ref_frags, prec_mz):
        """
        匹配碎片离子 - 修复版
        用库参考碎片去匹配实验观测碎片
        """
        if not ref_frags or not obs_frags:
            return []

        matched = []
        matched_refs = set()  # 记录已匹配的参考碎片

        for ref_mz in ref_frags:
            for obs_mz in obs_frags:
                da_error = abs(obs_mz - ref_mz)
                if da_error <= SECONDARY_DA_TOLERANCE:
                    # 排除接近母离子的碎片
                    if prec_mz and abs(ref_mz - prec_mz) <= 20:
                        continue
                    if ref_mz not in matched_refs:
                        matched.append(ref_mz)
                        matched_refs.add(ref_mz)
                    break

        return matched

    def identify(self):
        """执行化合物鉴定"""
        print("\n[4/5] 开始化合物鉴定...")

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

        print("\n[5/5] 合并正负离子结果并计算置信度...")
        self._merge_and_score()

        print(f"\n  完成鉴定: {len(self.results)} 个化合物")

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

        # 一级匹配：搜索ppm匹配的化合物
        candidates = self.db_idx.search(prec_mz, mode)

        for cand in candidates:
            info = cand['info']
            name = info['name_cn'] if info['name_cn'] != 'nan' else f"MZ_{prec_mz:.4f}"
            key = f"{name}_{cand['source_label']}"

            # 获取该化合物的总文献数
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
                    'total_lit_count': total_lit,  # 库中文献数
                    'is_priority': cand['is_priority'],
                    'source_label': cand['source_label'],
                    'matched_frags_pos': [],
                    'matched_frags_neg': [],
                    'obs_frags_pos': [],
                    'obs_frags_neg': [],
                    'ref_frags_pos': info['frag_pos'],
                    'ref_frags_neg': info['frag_neg'],
                    'ref_frag_count_pos': info['ref_frag_count_pos'],
                    'ref_frag_count_neg': info['ref_frag_count_neg'],
                    'diag_pos': info['diag_pos'],
                    'diag_neg': info['diag_neg'],
                    'has_frag_data': info['has_frag_data'],
                    'rt_pos': None,
                    'rt_neg': None,
                    'ppm_pos': 999,
                    'ppm_neg': 999,
                    'm/z_pos': None,
                    'm/z_neg': None,
                    'theo_mz_pos': None,
                    'theo_mz_neg': None,
                }

            if mode == 'positive':
                # 二级匹配：用参考碎片匹配实验碎片
                matched = self.match_fragments(obs_frags, cand['ref_frags'], prec_mz)
                self.results[key]['matched_frags_pos'].extend(matched)
                self.results[key]['obs_frags_pos'].extend(obs_frags)
                self.results[key]['rt_pos'] = rt
                self.results[key]['ppm_pos'] = min(self.results[key]['ppm_pos'], cand['ppm'])
                self.results[key]['m/z_pos'] = prec_mz
                self.results[key]['theo_mz_pos'] = cand['theoretical_mz']
            else:
                matched = self.match_fragments(obs_frags, cand['ref_frags'], prec_mz)
                self.results[key]['matched_frags_neg'].extend(matched)
                self.results[key]['obs_frags_neg'].extend(obs_frags)
                self.results[key]['rt_neg'] = rt
                self.results[key]['ppm_neg'] = min(self.results[key]['ppm_neg'], cand['ppm'])
                self.results[key]['m/z_neg'] = prec_mz
                self.results[key]['theo_mz_neg'] = cand['theoretical_mz']

    def _merge_and_score(self):
        """合并正负离子结果并计算置信度"""
        for key, res in self.results.items():
            # 去重匹配碎片
            all_matched = list(set(res['matched_frags_pos'] + res['matched_frags_neg']))
            all_obs = list(set(res['obs_frags_pos'] + res['obs_frags_neg']))
            all_ref = list(set(res['ref_frags_pos'] + res['ref_frags_neg']))

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
                theo_mz = res['theo_mz_pos']
                actual_mz = res['m/z_pos']
            elif has_neg:
                main_mode = 'negative'
                ppm = res['ppm_neg']
                theo_mz = res['theo_mz_neg']
                actual_mz = res['m/z_neg']
            else:
                min_ppm = min(res['ppm_pos'], res['ppm_neg'])
                ppm = min_ppm if min_ppm < 999 else 50
                theo_mz = res['theo_mz_pos'] or res['theo_mz_neg']
                actual_mz = res['m/z_pos'] or res['m/z_neg']
                main_mode = 'positive' if res['m/z_pos'] else 'negative'

            matched_count = len(all_matched)
            total_ref_count = len(all_ref)

            res['matched_frag_count'] = matched_count
            res['total_ref_frag_count'] = total_ref_count
            res['frag_match_ratio'] = matched_count / total_ref_count if total_ref_count > 0 else 0
            res['diag_count'] = diag_count
            res['diag_ions_list'] = diag_ions
            res['has_pos_neg'] = has_pos_neg
            res['rt_match'] = rt_match
            res['main_mode'] = main_mode
            res['ppm'] = ppm
            res['theoretical_mz'] = theo_mz
            res['actual_mz'] = actual_mz
            res['has_fragment_data'] = len(all_obs) > 0

            confidence = calculate_confidence(
                matched_count, total_ref_count, diag_count, res['total_lit_count'],
                has_pos_neg, rt_match, ppm, res['has_fragment_data']
            )
            res['confidence'] = confidence

            level_code, level_name, suggestion = get_confidence_level(confidence)
            res['level_code'] = level_code
            res['level_name'] = level_name
            res['suggestion'] = suggestion

            res['all_matched_frags'] = all_matched
            res['all_obs_frags'] = all_obs

    def export_reports(self, output_prefix='compound'):
        """导出鉴定报告"""
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')

        records = []
        for i, (key, res) in enumerate(sorted(self.results.items(),
                                              key=lambda x: (0 if x[1]['is_priority'] else 1, -x[1]['confidence'])), 1):
            matched_frags_str = '; '.join([f'{x:.4f}' for x in res.get('all_matched_frags', [])[:20]])
            obs_frags_str = '; '.join([f'{x:.4f}' for x in res.get('all_obs_frags', [])[:30]])
            diag_str = '; '.join([f'{x:.4f}' for x in res.get('diag_ions_list', [])])
            ref_frags_str = '; '.join([f'{x:.4f}' for x in (res.get('ref_frags_pos', []) + res.get('ref_frags_neg', []))[:20]])

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
                'm/z理论值': round(res.get('theoretical_mz', 0), 4) if res.get('theoretical_mz') else '',
                'ppm': round(res.get('ppm', 0), 4),
                '是否有碎片数据': '是' if res.get('has_fragment_data') else '否',
                '库中参考碎片数': res.get('total_ref_frag_count', 0),
                '匹配碎片数': res.get('matched_frag_count', 0),
                '碎片匹配比例': f"{res.get('frag_match_ratio', 0)*100:.2f}%",
                '主要碎片离子': obs_frags_str if obs_frags_str else ref_frags_str,
                '参考碎片离子': ref_frags_str,
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

        priority_df = df[df['来源分类'] != '其他来源'].copy()
        if not priority_df.empty:
            pri_path = f"{OUTPUT_DIR}/{output_prefix}_priority_report_{timestamp}.csv"
            priority_df.to_csv(pri_path, index=False, encoding='utf-8-sig')
            print(f"✅ 优先级报告: {pri_path} ({len(priority_df)} 个化合物)")

        other_df = df[df['来源分类'] == '其他来源'].copy()
        if not other_df.empty:
            oth_path = f"{OUTPUT_DIR}/{output_prefix}_other_report_{timestamp}.csv"
            other_df.to_csv(oth_path, index=False, encoding='utf-8-sig')
            print(f"✅ 其他来源报告: {oth_path} ({len(other_df)} 个化合物)")

        # 按评级分组
        for level in ['I', 'II', 'III', 'IV', 'V']:
            level_df = df[df['评级'] == level].copy()
            if not level_df.empty:
                level_name = get_confidence_level(level_df['置信度'].iloc[0])[1]
                level_path = f"{OUTPUT_DIR}/{output_prefix}_level{level}_{timestamp}.csv"
                level_df.to_csv(level_path, index=False, encoding='utf-8-sig')
                print(f"✅ {level_name}报告: {level_path} ({len(level_df)} 个化合物)")

        return full_path


# ============================================================================
# 主函数
# ============================================================================
def run_analysis(herb_name, ms_pos_file, ms_neg_file, db_file):
    """运行分析的便捷函数"""
    priority_herbs = [h.strip() for h in herb_name.split(',')]
    for fp in [ms_pos_file, ms_neg_file, db_file]:
        if fp and not os.path.exists(fp):
            alt = os.path.join('/workspace/user_input_files', os.path.basename(fp))
            if os.path.exists(alt):
                if fp == ms_pos_file:
                    ms_pos_file = alt
                elif fp == ms_neg_file:
                    ms_neg_file = alt
                else:
                    db_file = alt

    identifier = UniversalIdentifier(ms_pos_file, ms_neg_file, db_file, priority_herbs)
    identifier.identify()

    prefix = '_'.join(priority_herbs) if priority_herbs else 'compound'
    output_path = identifier.export_reports(prefix)

    return identifier.results, output_path


def main():
    MS_POS = '栀子-MS+.xlsx'
    MS_NEG = '栀子-MS-.xlsx'
    DB = 'TCM-SM-MS DB.CSV'

    import sys
    if len(sys.argv) > 1:
        herbs_str = sys.argv[1]
        priority_herbs = [h.strip() for h in herbs_str.split(',')]
    else:
        priority_herbs = ['栀子']

    print("\n" + "=" * 80)
    print(f"【中药化合物智能鉴定系统 v11.0 - 融合版】")
    print(f"分析药材: {', '.join(priority_herbs)}")
    print("=" * 80 + "\n")

    results, output_path = run_analysis(','.join(priority_herbs), MS_POS, MS_NEG, DB)

    print("\n" + "=" * 80)
    print("【鉴定结果摘要】")
    print("=" * 80)

    level_stats = {}
    source_stats = {'priority': 0, 'other': 0}
    frag_match_stats = {'>=3': 0, '<3': 0}
    lit_stats = {'>=1': 0, '<1': 0}

    for res in results.values():
        level = res.get('level_name', '未知')
        level_stats[level] = level_stats.get(level, 0) + 1
        if res.get('is_priority'):
            source_stats['priority'] += 1
        else:
            source_stats['other'] += 1
        if res.get('matched_frag_count', 0) >= 3:
            frag_match_stats['>=3'] += 1
        else:
            frag_match_stats['<3'] += 1
        if res.get('total_lit_count', 0) >= 1:
            lit_stats['>=1'] += 1
        else:
            lit_stats['<1'] += 1

    for level in ['确证级', '高置信级', '推定级', '提示级', '排除级']:
        count = level_stats.get(level, 0)
        if count > 0:
            print(f"  {level}: {count} 个化合物")

    print(f"\n  匹配碎片数≥3: {frag_match_stats['>=3']} 个")
    print(f"  匹配碎片数<3: {frag_match_stats['<3']} 个")
    print(f"  有文献来源: {lit_stats['>=1']} 个")
    print(f"  无文献来源: {lit_stats['<1']} 个")
    print(f"\n  优先级药材来源: {source_stats['priority']} 个化合物")
    print(f"  其他来源: {source_stats['other']} 个化合物")
    print(f"\n  总计: {len(results)} 个化合物")
    print(f"  报告文件: {output_path}")
    print("=" * 80)

    return results


if __name__ == '__main__':
    main()
