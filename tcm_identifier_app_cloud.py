# -*- coding: utf-8 -*-
"""
中药化合物智能鉴定平台 - GitHub主程序 v1.0
============================================================
功能说明：
- 整合中药化合物鉴定核心算法
- 详细碎片比对报告：列出所有碎片离子（检测碎片+库中匹配碎片）
- 碎片映射到库中位置（药材来源、文献来源）
- 支持多种加和离子扩展匹配

核心逻辑：
1. 一级匹配：前体离子ppm容差匹配
2. 二级匹配：碎片离子100Da容差匹配（不使用ppm）

报告输出：
- 检测碎片离子：实验数据中的碎片离子
- 库中匹配碎片离子：数据库中与检测碎片匹配上的参考碎片
- 碎片来源映射：每个碎片对应的药材来源和文献来源
"""

import pandas as pd
import numpy as np
import os
import re
from datetime import datetime
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')


# ============================================================================
# 配置参数
# ============================================================================

# 一级匹配容差（ppm）
PRIMARY_PPM_TOLERANCE = 50

# 二级碎片匹配容差（Da）- 用户核心要求
SECONDARY_DA_TOLERANCE = 100

# 保留时间容差（min）
RT_TOLERANCE = 0.2

# 输出目录
OUTPUT_DIR = '/workspace/output'
os.makedirs(OUTPUT_DIR, exist_ok=True)


# ============================================================================
# 加和离子扩展常量
# ============================================================================

ADDUCTS_POSITIVE = {
    '[M+H]+': 1.0073,
    '[M+Na]+': 22.9898,
    '[M+K]+': 38.9637,
    '[M+NH4]+': 18.0344,
    '[M+H-H2O]+': -17.0027,
}

ADDUCTS_NEGATIVE = {
    '[M-H]-': -1.0073,
    '[M+Na-2H]-': 20.9747,
    '[M+K-2H]-': 36.9556,
    '[M-H+H2O]-': 18.0106,
}


# ============================================================================
# 辅助函数
# ============================================================================

def normalize_formula(formula):
    """标准化分子式"""
    if pd.isna(formula) or not formula:
        return 'unknown'
    formula = str(formula)
    subscripts = {'₀': '0', '₁': '1', '₂': '2', '₃': '3', '₄': '4',
                  '₅': '5', '₆': '6', '₇': '7', '₈': '8', '₉': '9',
                  'ₐ': 'a', 'ₑ': 'e', 'ₕ': 'h', 'ₖ': 'k', 'ₙ': 'n',
                  'ₒ': 'o', 'ₚ': 'p', 'ₛ': 's', 'ₜ': 't', 'ᵢ': 'i', 'ᵣ': 'r'}
    for sub, norm in subscripts.items():
        formula = formula.replace(sub, norm)
    return formula.strip()


def parse_fragments(frag_str):
    """解析碎片离子字符串"""
    if pd.isna(frag_str) or not frag_str:
        return []
    frags = []
    for part in re.split(r'[;,、；\s\n\r]+', str(frag_str)):
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


def load_database(db_path):
    """加载主数据库"""
    if not os.path.exists(db_path):
        print(f"错误: 数据库文件不存在 - {db_path}")
        return pd.DataFrame()

    try:
        if db_path.lower().endswith('.csv'):
            df = pd.read_csv(db_path, encoding='utf-8-sig')
        else:
            df = pd.read_excel(db_path)

        # 列名标准化
        col_mapping = {
            '药材名': '药材名称',
            '文献': '文献来源',
            '中文名': '名称（中文）',
            '英文名': '名称（英文）',
            '保留时间(min)': '保留时间',
            '中丢失': '中性丢失',
            '准分子离子(正离子模式)': '准分子离子（正）',
            '碎片离子(正离子模式)': '碎片离子（正）',
            '准分子离子(负离子模式)': '准分子离子（负）',
            '碎片离子(负离子模式)': '碎片离子（负）'
        }
        df = df.rename(columns=col_mapping)

        print(f"  主数据库加载成功: {len(df)} 条记录")
        return df
    except Exception as e:
        print(f"  数据库加载失败: {e}")
        return pd.DataFrame()


def load_ms_data(file_path):
    """加载质谱数据"""
    if not file_path or not os.path.exists(file_path):
        print(f"  质谱文件不存在或未提供: {file_path}")
        return pd.DataFrame()

    try:
        if file_path.endswith('.xlsx'):
            df = pd.read_excel(file_path)
        else:
            df = pd.read_csv(file_path)
        print(f"  质谱数据加载成功: {len(df)} 条记录")
        return df
    except Exception as e:
        print(f"  质谱数据加载失败: {e}")
        return pd.DataFrame()


# ============================================================================
# 数据库索引构建（支持加和离子扩展）
# ============================================================================

class DatabaseIndex:
    """
    数据库索引类 - 融合版
    支持加和离子扩展，同时记录碎片来源映射
    """
    def __init__(self, df, herb_col, is_priority, source_label):
        self.df = df
        self.herb_col = herb_col
        self.is_priority = is_priority
        self.source_label = source_label
        self.records = []
        self._build_index()

    def _build_index(self):
        """构建索引 - 同时添加原始m/z和扩展加和离子"""
        for _, row in self.df.iterrows():
            mz_p = row.get('准分子离子（正）', 0)
            mz_n = row.get('准分子离子（负）', 0)
            source = str(row.get('文献来源', ''))
            herb = str(row.get(self.herb_col, ''))
            name = str(row.get('名称（中文）', ''))
            if name == 'nan' or not name.strip():
                name = str(row.get('名称（英文）', ''))

            lit_count = len([s for s in source.split(';') if s.strip()])

            frag_pos = parse_fragments(row.get('碎片离子（正）', ''))
            frag_neg = parse_fragments(row.get('碎片离子（负）', ''))
            diag_pos = parse_diagnostic_ions(frag_pos)
            diag_neg = parse_diagnostic_ions(frag_neg)

            info = {
                'name_cn': name,
                'name_en': str(row.get('名称（英文）', '')),
                'formula': normalize_formula(str(row.get('分子式', 'unknown'))),
                'herb': herb,
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

            # 正离子索引 - 原始值 + 加和离子扩展
            if pd.notna(mz_p) and str(mz_p).strip():
                try:
                    base = float(mz_p)
                    # 添加原始值
                    self._add_record(base, 'positive', info, frag_pos, diag_pos)
                    # 添加加和离子扩展
                    for shift in ADDUCTS_POSITIVE.values():
                        self._add_record(base + shift, 'positive', info, frag_pos, diag_pos)
                except:
                    pass

            # 负离子索引 - 原始值 + 加和离子扩展
            if pd.notna(mz_n) and str(mz_n).strip():
                try:
                    base = float(mz_n)
                    self._add_record(base, 'negative', info, frag_neg, diag_neg)
                    for shift in ADDUCTS_NEGATIVE.values():
                        self._add_record(base + shift, 'negative', info, frag_neg, diag_neg)
                except:
                    pass

    def _add_record(self, mz, mode, info, ref_frags, ref_diag):
        """添加索引记录"""
        self.records.append({
            'mz': mz,
            'mode': mode,
            'info': info,
            'ref_frags': ref_frags,
            'ref_diag': ref_diag,
            'source_label': self.source_label,
        })


class MultiPriorityDB:
    """多优先级数据库管理"""
    def __init__(self, db_df, priority_herbs):
        self.priority_herbs = priority_herbs
        self.herb_col = '药材名称' if '药材名称' in db_df.columns else '药材名'
        self.priority_indexes = {}
        self.other_index = None
        self.compound_lit_count = db_df.groupby('名称（中文）').size().to_dict()
        self._build(db_df)

    def _build(self, db_df):
        """构建优先级索引"""
        remaining = db_df.copy()
        print('\n  构建优先级索引（加和离子扩展模式）:')

        for herb in self.priority_herbs:
            mask = remaining[self.herb_col].str.contains(herb, na=False, case=False)
            herb_df = remaining[mask].copy()
            remaining = remaining[~mask]
            print(f'    [{herb}]来源: {len(herb_df)} 条')
            self.priority_indexes[herb] = DatabaseIndex(herb_df, self.herb_col, True, f'[{herb}]来源')

        print(f'    [其他]来源: {len(remaining)} 条')
        self.other_index = DatabaseIndex(remaining, self.herb_col, False, '其他来源')

    def search(self, prec_mz, mode):
        """搜索匹配的候选化合物"""
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
# 碎片匹配核心逻辑（来自v11.0原始逻辑）
# ============================================================================

class TCMCompoundIdentifier:
    """
    中药化合物智能鉴定系统 - GitHub主程序 v1.0

    核心设计：
    1. 一级匹配：前体离子ppm容差匹配
    2. 二级匹配：碎片离子Da容差匹配（不使用ppm）
    3. 详细碎片报告：列出所有检测碎片和库中匹配碎片，映射到库中位置
    """

    def __init__(self, db_path, ms_pos_path, ms_neg_path, priority_herbs):
        """初始化鉴定系统"""
        self.priority_herbs = priority_herbs
        print('=' * 80)
        print('中药化合物智能鉴定系统 - GitHub主程序 v1.0')
        print('=' * 80)
        print(f'优先级药材: {", ".join(priority_herbs)}')
        print(f'一级母离子ppm容差: <={PRIMARY_PPM_TOLERANCE}ppm')
        print(f'二级碎片离子容差: {SECONDARY_DA_TOLERANCE}Da（用户核心要求）')
        print('加和离子扩展: [M+H]+, [M+Na]+, [M+K]+, [M+NH4]+, [M-H]-等')
        print('=' * 80)

        print('\n[1/5] 加载数据库...')
        self.db = load_database(db_path)
        if self.db.empty:
            print("错误: 数据库为空，程序退出")
            # 初始化空数据以避免属性错误
            self.ms_pos = pd.DataFrame()
            self.ms_neg = pd.DataFrame()
            self.db_idx = None
            self.results = {}
            self._frag_match_cache = {}
            return

        print('\n[2/5] 加载质谱数据...')
        self.ms_pos = load_ms_data(ms_pos_path)
        self.ms_neg = load_ms_data(ms_neg_path)

        print('\n[3/5] 构建优先级数据库索引（加和离子扩展）...')
        self.db_idx = MultiPriorityDB(self.db, priority_herbs)

        self.results = {}
        self._frag_match_cache = {}  # 碎片匹配缓存

    def match_fragments(self, obs_frags, ref_frags, prec_mz):
        """
        碎片匹配核心逻辑（来自v11.0原始逻辑）

        用户核心要求：使用Da容差，不使用ppm

        输入：
        - obs_frags: 观测到的碎片离子列表
        - ref_frags: 数据库中的参考碎片离子列表
        - prec_mz: 前体离子m/z

        输出：
        - matched_obs: 匹配上的观测碎片
        - matched_ref: 匹配上的参考碎片
        - match_details: 匹配详情（包含来源映射）
        """
        matched_obs = []
        matched_ref = []
        match_details = []

        if not ref_frags or not obs_frags:
            return matched_obs, matched_ref, match_details

        matched_ref_set = set()
        matched_obs_set = set()

        for ref_mz in ref_frags:
            for obs_mz in obs_frags:
                # 使用Da容差进行碎片匹配（用户核心要求）
                da_error = abs(obs_mz - ref_mz)
                if da_error <= SECONDARY_DA_TOLERANCE:
                    # 排除母离子本身
                    if prec_mz and abs(ref_mz - prec_mz) <= 20:
                        continue

                    # 避免重复匹配
                    if ref_mz not in matched_ref_set:
                        matched_ref.append(ref_mz)
                        matched_ref_set.add(ref_mz)
                    if obs_mz not in matched_obs_set:
                        matched_obs.append(obs_mz)
                        matched_obs_set.add(obs_mz)

                    # 记录匹配详情
                    match_details.append({
                        'obs_mz': obs_mz,
                        'ref_mz': ref_mz,
                        'da_error': da_error,
                        'match_type': 'exact' if da_error == 0 else 'tolerance'
                    })
                    break

        return matched_obs, matched_ref, match_details

    def identify(self):
        """执行化合物鉴定"""
        print('\n[4/5] 开始化合物鉴定...')
        pos_data = self._extract_precursors(self.ms_pos, 'positive')
        neg_data = self._extract_precursors(self.ms_neg, 'negative')
        print(f'  正离子数据: {len(pos_data)} 个母离子')
        print(f'  负离子数据: {len(neg_data)} 个母离子')

        print('  处理正离子数据...')
        for i, prec in enumerate(pos_data):
            if (i + 1) % 500 == 0:
                print(f'    正离子进度: {i+1}/{len(pos_data)}')
            self._process_precursor(prec, 'positive')

        print('  处理负离子数据...')
        for i, prec in enumerate(neg_data):
            if (i + 1) % 500 == 0:
                print(f'    负离子进度: {i+1}/{len(neg_data)}')
            self._process_precursor(prec, 'negative')

        print('\n[5/5] 合并正负离子结果并计算置信度...')
        self._merge_and_score()
        print(f'\n  完成鉴定: {len(self.results)} 个化合物')
        return self.results

    def _extract_precursors(self, df, mode):
        """提取母离子和碎片离子"""
        if df.empty:
            return []

        prec_col = None
        for col in df.columns:
            if 'Precursor' in col and 'M/z' in col:
                prec_col = col
                break
        if prec_col is None:
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
        candidates = self.db_idx.search(prec_mz, mode)

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
                    # 碎片数据
                    'matched_frags_pos': [],
                    'matched_frags_neg': [],
                    'obs_frags_pos': [],
                    'obs_frags_neg': [],
                    'ref_frags_pos': info['frag_pos'],
                    'ref_frags_neg': info['frag_neg'],
                    'ref_frag_count_pos': info['ref_frag_count_pos'],
                    'ref_frag_count_neg': info['ref_frag_count_neg'],
                    # 碎片匹配详情（包含来源映射）
                    'frag_match_details': [],
                    # 诊断离子
                    'diag_pos': info['diag_pos'],
                    'diag_neg': info['diag_neg'],
                    'has_frag_data': info['has_frag_data'],
                    # RT和ppm
                    'rt_pos': None,
                    'rt_neg': None,
                    'ppm_pos': 999,
                    'ppm_neg': 999,
                    'm/z_pos': None,
                    'm/z_neg': None,
                    'theo_mz_pos': None,
                    'theo_mz_neg': None,
                }

            res = self.results[key]
            frag_match_key = f"{key}_{mode}"
            already_matched = frag_match_key in self._frag_match_cache

            if mode == 'positive':
                if not already_matched:
                    matched_obs, matched_ref, match_details = self.match_fragments(
                        obs_frags, cand['ref_frags'], prec_mz)
                    existing = set(res['matched_frags_pos'])
                    new_frags = [f for f in matched_obs if f not in existing]
                    res['matched_frags_pos'].extend(new_frags)

                    # 保存碎片匹配详情
                    for detail in match_details:
                        detail['mode'] = 'positive'
                        detail['herb'] = info['herb']
                        detail['source'] = info['source']
                    res['frag_match_details'].extend(match_details)

                    self._frag_match_cache[frag_match_key] = True
                res['obs_frags_pos'].extend(obs_frags)
                res['rt_pos'] = rt
                res['ppm_pos'] = min(res['ppm_pos'], cand['ppm'])
                res['m/z_pos'] = prec_mz
                res['theo_mz_pos'] = cand['theoretical_mz']
            else:
                if not already_matched:
                    matched_obs, matched_ref, match_details = self.match_fragments(
                        obs_frags, cand['ref_frags'], prec_mz)
                    existing = set(res['matched_frags_neg'])
                    new_frags = [f for f in matched_obs if f not in existing]
                    res['matched_frags_neg'].extend(new_frags)

                    for detail in match_details:
                        detail['mode'] = 'negative'
                        detail['herb'] = info['herb']
                        detail['source'] = info['source']
                    res['frag_match_details'].extend(match_details)

                    self._frag_match_cache[frag_match_key] = True
                res['obs_frags_neg'].extend(obs_frags)
                res['rt_neg'] = rt
                res['ppm_neg'] = min(res['ppm_neg'], cand['ppm'])
                res['m/z_neg'] = prec_mz
                res['theo_mz_neg'] = cand['theoretical_mz']

    def _merge_and_score(self):
        """合并结果并计算置信度"""
        for key, res in self.results.items():
            # 合并正负离子碎片
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

            # 选择主模式
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
            has_fragment_data = len(all_obs) > 0

            # 计算置信度
            confidence = self._calculate_confidence(
                matched_count, total_ref_count, diag_count, res['total_lit_count'],
                has_pos_neg, rt_match, ppm, has_fragment_data,
                res.get('is_priority', False)
            )

            level_code, level_name, suggestion = self._get_confidence_level(confidence)

            # 更新结果
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
            res['confidence'] = confidence
            res['level_code'] = level_code
            res['level_name'] = level_name
            res['suggestion'] = suggestion
            res['all_matched_frags'] = all_matched
            res['all_obs_frags'] = all_obs
            res['all_ref_frags'] = all_ref

    def _calculate_confidence(self, matched_frag_count, total_ref_frags, diag_count, lit_count,
                              has_pos_neg, rt_match, ppm, has_fragment_data, is_priority=False):
        """计算置信度得分"""
        score = 0

        # 文献数量得分
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

        # 正负离子确认
        if has_pos_neg:
            score += 25
        else:
            score += 8

        # ppm精确度
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

        # 碎片匹配
        if total_ref_frags > 0:
            match_ratio = min(matched_frag_count / total_ref_frags, 1.0)
            score += match_ratio * 15
        elif has_fragment_data:
            score += 2

        # RT匹配
        if rt_match:
            score += 10

        # 诊断离子
        if diag_count >= 5:
            score += 5
        elif diag_count >= 3:
            score += 4
        elif diag_count >= 2:
            score += 3
        elif diag_count >= 1:
            score += 2

        # 优先级加成
        if is_priority:
            score += 60

        return min(score, 100)

    def _get_confidence_level(self, score):
        """获取置信度等级"""
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

    def _format_fragment_mapping(self, res):
        """
        格式化碎片映射报告（核心功能）

        输出两行详细信息：
        1. 检测碎片离子：实验数据中检测到的碎片
        2. 库中匹配碎片离子：与检测碎片匹配上的数据库参考碎片

        每个碎片映射到库中位置（药材来源）
        """
        lines = []

        # 第一行：检测碎片离子
        obs_frags = sorted(res.get('all_obs_frags', []))
        if obs_frags:
            obs_str = '; '.join([f'{x:.4f}' for x in obs_frags])
            lines.append(f"【检测碎片】{obs_str}")

        # 第二行：库中匹配碎片离子（与检测碎片匹配上的参考碎片）
        matched_frags = sorted(res.get('all_matched_frags', []))
        if matched_frags:
            # 构建碎片到来源的映射
            frag_sources = defaultdict(set)
            for detail in res.get('frag_match_details', []):
                ref_mz = detail['ref_mz']
                herb = detail.get('herb', '')
                if herb and herb != 'nan':
                    frag_sources[ref_mz].add(herb)

            # 格式化匹配碎片及来源
            frag_detail_list = []
            for frag in matched_frags:
                sources = frag_sources.get(frag, set())
                if sources:
                    src_str = '; '.join(list(sources)[:3])  # 最多显示3个来源
                    frag_detail_list.append(f"{frag:.4f}({src_str})")
                else:
                    frag_detail_list.append(f"{frag:.4f}")

            matched_str = '; '.join(frag_detail_list)
            lines.append(f"【库中匹配碎片】{matched_str}")

        return '\n'.join(lines)

    def export_reports(self, output_prefix='compound'):
        """
        导出鉴定报告

        报告包含：
        1. 基础信息
        2. 碎片比对详情（含来源映射）
        3. 检测碎片和库中匹配碎片两行显示
        """
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')

        records = []
        for i, (key, res) in enumerate(sorted(self.results.items(),
                                              key=lambda x: (0 if x[1]['is_priority'] else 1, -x[1]['confidence'])), 1):
            # 格式化碎片报告
            frag_report = self._format_fragment_mapping(res)

            # 获取检测碎片列表
            obs_frags = sorted(res.get('all_obs_frags', []))
            obs_frags_str = '; '.join([f'{x:.4f}' for x in obs_frags[:30]]) if obs_frags else '无'

            # 获取库中匹配碎片列表
            matched_frags = sorted(res.get('all_matched_frags', []))
            matched_frags_str = '; '.join([f'{x:.4f}' for x in matched_frags[:20]]) if matched_frags else '无'

            # 获取所有参考碎片
            ref_frags = sorted(res.get('all_ref_frags', []))
            ref_frags_str = '; '.join([f'{x:.4f}' for x in ref_frags[:20]]) if ref_frags else '无'

            # 诊断离子
            diag_str = '; '.join([f'{x:.4f}' for x in res.get('diag_ions_list', [])])

            row = {
                '序号': i,
                '出峰时间t/min': res.get('rt_pos') or res.get('rt_neg', ''),
                '化合物中文名': res['name'],
                '化合物英文名': res['name_en'],
                '分子式': res['formula'],
                'CAS号': res['cas'],
                '离子化方式': res['main_mode'].upper() if res.get('main_mode') else '',
                'm/z实际值': round(res.get('actual_mz', 0), 4) if res.get('actual_mz') else '',
                'm/z理论值': round(res.get('theoretical_mz', 0), 4) if res.get('theoretical_mz') else '',
                'ppm': round(res.get('ppm', 0), 4),
                # 碎片信息
                '检测碎片离子(m/z)': obs_frags_str,
                '库中匹配碎片离子(m/z)': matched_frags_str,
                '库中参考碎片离子(m/z)': ref_frags_str,
                '匹配碎片数': res.get('matched_frag_count', 0),
                '碎片匹配比例': f"{res.get('frag_match_ratio', 0)*100:.2f}%",
                # 详细碎片报告
                '碎片比对详情': frag_report,
                # 其他信息
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

        # 完整报告
        full_path = f"{OUTPUT_DIR}/{output_prefix}_full_report_{timestamp}.csv"
        df.to_csv(full_path, index=False, encoding='utf-8-sig')
        print(f'\n完整报告: {full_path} ({len(df)} 个化合物)')

        # 优先级报告
        priority_df = df[df['来源分类'] != '其他来源'].copy()
        if not priority_df.empty:
            pri_path = f"{OUTPUT_DIR}/{output_prefix}_priority_report_{timestamp}.csv"
            priority_df.to_csv(pri_path, index=False, encoding='utf-8-sig')
            print(f'优先级报告: {pri_path} ({len(priority_df)} 个化合物)')

        # 置信度分级报告
        for level in ['I', 'II', 'III', 'IV', 'V']:
            level_df = df[df['评级'] == level].copy()
            if not level_df.empty:
                level_name = self._get_confidence_level(level_df['置信度'].iloc[0])[1]
                level_path = f"{OUTPUT_DIR}/{output_prefix}_level{level}_{timestamp}.csv"
                level_df.to_csv(level_path, index=False, encoding='utf-8-sig')
                print(f'{level_name}报告: {level_path} ({len(level_df)} 个化合物)')

        # 打印置信度分布
        print('\n--- 置信度分布 ---')
        print(df['评级名称'].value_counts().sort_index())

        return full_path, df


def run_gardenia_analysis():
    """运行栀子分析"""
    DB_FILE = '/workspace/user_input_files/TCM-SM-MS DB.CSV'
    MS_POS_FILE = '/workspace/user_input_files/栀子-MS+.xlsx'
    MS_NEG_FILE = '/workspace/user_input_files/栀子-MS-.xlsx'

    print('\n\n' + '=' * 80)
    print('【栀子化合物鉴定分析 - GitHub主程序 v1.0】')
    print('=' * 80)

    identifier = TCMCompoundIdentifier(
        db_path=DB_FILE,
        ms_pos_path=MS_POS_FILE,
        ms_neg_path=MS_NEG_FILE,
        priority_herbs=['栀子']
    )

    identifier.identify()
    output_path, report_df = identifier.export_reports('栀子_github_v1')

    print('\n\n' + '=' * 80)
    print('【分析完成】')
    print('=' * 80)
    print(f'报告路径: {output_path}')
    print(f'鉴定化合物总数: {len(report_df)}')

    return identifier.results, output_path, report_df


# ============================================================================
# 主程序入口
# ============================================================================

if __name__ == '__main__':
    results, path, report = run_gardenia_analysis()
