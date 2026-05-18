#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
中药化合物智能鉴定平台 v2.0（新评分体系版）
==========================================

评分体系：
- 基础分100分：碎片匹配60分 + 文献来源30分 + ppm误差5分 + 正负离子5分
- 优先级加成：50分
- 诊断离子：不列入评分

评级标准：
- 来源药材：ppm≤50 AND 碎片数≥2 → 确证级
- 非来源药材：碎片≥20 AND 文献数≥5 AND ppm≤15 → 确证级

Author: MiniMax Agent
"""

import streamlit as st
import pandas as pd
import numpy as np
import re
import os
import time
import tempfile
from datetime import datetime
from io import BytesIO
from typing import List, Dict, Tuple, Optional, Set
from bisect import bisect_left, bisect_right
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# 登录验证常量
# ============================================================================
VALID_USERNAME = "ZY"
VALID_PASSWORD = "513513"

# ============================================================================
# 登录页面函数
# ============================================================================
def login_page():
    """显示登录页面"""
    st.markdown("""
    <style>
        .login-container {
            max-width: 450px;
            margin: 100px auto;
            padding: 2.5rem;
            background: white;
            border-radius: 24px;
            box-shadow: 0 20px 60px rgba(0, 0, 0, 0.15);
        }
        .login-title {
            text-align: center;
            font-size: 1.8rem;
            font-weight: 700;
            margin-bottom: 2rem;
            background: linear-gradient(135deg, #059669, #0891b2);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
        }
        .login-subtitle {
            text-align: center;
            color: #666;
            margin-bottom: 2rem;
            font-size: 0.95rem;
        }
    </style>
    <div class="login-container">
        <div class="login-title">🌿 中药化合物智能鉴定平台</div>
        <div class="login-subtitle">v2.0 | 新评分体系版</div>
    </div>
    """, unsafe_allow_html=True)

    with st.form("login_form", clear_on_submit=True):
        st.markdown("---")
        username = st.text_input("👤 用户名", placeholder="请输入用户名", label_visibility="collapsed")
        password = st.text_input("🔐 密码", type="password", placeholder="请输入密码", label_visibility="collapsed")
        st.markdown("---")
        submitted = st.form_submit_button("🚀 登录", use_container_width=True)

        if submitted:
            if username == VALID_USERNAME and password == VALID_PASSWORD:
                st.session_state.logged_in = True
                st.session_state.username = username
                st.success("✅ 登录成功！正在跳转...")
                time.sleep(0.5)
                st.rerun()
            else:
                st.error("❌ 用户名或密码错误，请重试")

def logout_button():
    """显示登出按钮"""
    if st.sidebar.button("🚪 登出", use_container_width=True, type="primary"):
        st.session_state.logged_in = False
        st.session_state.pop('username', None)
        st.rerun()

def check_login():
    """检查登录状态"""
    if 'logged_in' not in st.session_state:
        st.session_state.logged_in = False
    if not st.session_state.logged_in:
        login_page()
        return False
    return True

# ============================================================================
# 内置数据库配置
# ============================================================================
_CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
BUILTIN_DB_PATH = os.path.join(_CURRENT_DIR, 'TCM-SM-MS DB.CSV')

# ============================================================================
# 配置参数
# ============================================================================
PRIMARY_PPM_TOLERANCE = 50
SECONDARY_DA_TOLERANCE = 0.15
RT_TOLERANCE = 0.2

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
        return df
    except Exception as e:
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
        return df
    except Exception as e:
        return pd.DataFrame()

def parse_fragments_with_source(frag_str, source, db_source='', min_mz=None):
    """解析碎片离子字符串

    Args:
        frag_str: 碎片离子字符串
        source: 文献来源
        db_source: 数据库来源标签
        min_mz: 最小m/z阈值，如果为None则不限制
    """
    frags = []
    source_map = {}
    for part in re.split(r'[;、]', str(frag_str) if pd.notna(frag_str) else ''):
        part = part.strip()
        if part:
            match = re.match(r'^([\d.]+)', part)
            if match:
                try:
                    mz_val = float(match.group(1))
                    # 【修改】如果设置了min_mz，则过滤掉m/z < min_mz的碎片
                    if min_mz is not None and mz_val < min_mz:
                        continue
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

# ============================================================================
# 新评分体系（v2.0）
# ============================================================================
def calculate_confidence_v2(matched_frag_count, total_ref_frags, lit_count,
                           has_pos_neg, ppm, is_priority=False):
    """
    计算置信度得分（新评分体系v2.0）

    基础分100分：
    - 碎片匹配：最高60分
    - 文献来源：最高30分
    - ppm误差：最高5分
    - 正负离子确认：最高5分

    优先级加成：50分
    """
    base_score = 0

    if matched_frag_count == 0:
        return base_score

    # 1. 碎片匹配（最高60分）
    if matched_frag_count >= 10:
        base_score += 60
    elif matched_frag_count >= 8:
        base_score += 54
    elif matched_frag_count >= 6:
        base_score += 48
    elif matched_frag_count >= 5:
        base_score += 42
    elif matched_frag_count >= 4:
        base_score += 36
    elif matched_frag_count >= 3:
        base_score += 30
    elif matched_frag_count >= 2:
        base_score += 24
    elif matched_frag_count >= 1:
        base_score += 18

    # 2. 文献来源（最高30分）
    if lit_count >= 15:
        base_score += 30
    elif lit_count >= 12:
        base_score += 27
    elif lit_count >= 10:
        base_score += 24
    elif lit_count >= 7:
        base_score += 21
    elif lit_count >= 5:
        base_score += 18
    elif lit_count >= 3:
        base_score += 15
    elif lit_count >= 2:
        base_score += 12
    elif lit_count >= 1:
        base_score += 9

    # 3. ppm误差（最高5分）
    if ppm <= 10:
        base_score += 5
    elif ppm <= 30:
        base_score += 3
    elif ppm <= 50:
        base_score += 1

    # 4. 正负离子同时确认（最高5分）
    if has_pos_neg:
        base_score += 5

    base_score = min(base_score, 100)
    priority_bonus = 50 if is_priority else 0
    return base_score + priority_bonus

def get_confidence_level_v2(matched_frag_count=0, ppm=999, lit_count=0, is_source_herb_match=False):
    """
    根据条件确定评级（新评分体系v2.0）

    - 来源药材：ppm≤50 AND 碎片数≥2 → 确证级
    - 非来源药材：碎片≥20 AND 文献数≥5 AND ppm≤15 → 确证级
    """
    source_match = is_source_herb_match

    # 来源药材标准（宽松）
    if source_match:
        if matched_frag_count >= 2 and ppm <= 50:
            return 'I', '确证级', '高置信度，来源药材匹配'
        if matched_frag_count >= 1 and ppm <= 50:
            return 'II', '高置信级', '来源药材匹配'

    # 非来源药材标准（严格）
    if not source_match:
        if matched_frag_count >= 20 and lit_count >= 5 and ppm <= 15:
            return 'I', '确证级', '高置信度，可作为定性依据'
        if matched_frag_count >= 10 and lit_count >= 2 and ppm <= 50:
            return 'II', '高置信级', '较强置信度，建议进一步验证'

    # 通用标准
    if matched_frag_count >= 2:
        return 'III', '推定级', '中等置信度，需要更多证据支持'
    if matched_frag_count >= 1:
        return 'IV', '提示级', '低置信度，仅供参考'
    return 'V', '排除级', '置信度不足，建议排除'

# ============================================================================
# 高效数据库索引
# ============================================================================
class EfficientDatabaseIndex:
    def __init__(self, df, herb_col, is_priority, source_label):
        self.df = df
        self.herb_col = herb_col
        self.is_priority = is_priority
        self.source_label = source_label

        self.sorted_idx_pos = []
        self.sorted_idx_neg = []
        self.mz_values_pos = np.array([])
        self.mz_values_neg = np.array([])
        self.frag_source_lookup_pos = {}
        self.frag_source_lookup_neg = {}

        self._build_index()

    def _build_index(self):
        for _, row in self.df.iterrows():
            mz_p = row.get('准分子离子（正）', 0)
            mz_n = row.get('准分子离子（负）', 0)
            source = str(row.get('文献来源', ''))
            herb_name = str(row.get(self.herb_col, ''))

            # 【修改】判断是否为对照品：药材名中包含"对照品"的记录，碎片离子m/z<100不参与匹配
            # 其他记录：所有碎片都参与匹配
            is_reference_standard = '对照品' in herb_name
            min_mz = 100 if is_reference_standard else None

            frag_pos, source_map_pos = parse_fragments_with_source(
                row.get('碎片离子（正）', ''), source, self.source_label, min_mz=min_mz
            )
            frag_neg, source_map_neg = parse_fragments_with_source(
                row.get('碎片离子（负）', ''), source, self.source_label, min_mz=min_mz
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

            if pd.notna(mz_p) and str(mz_p).strip():
                try:
                    base = float(mz_p)
                    for shift in [0] + list(ADDUCTS_POSITIVE.values()):
                        mz_val = base + shift
                        self.sorted_idx_pos.append((mz_val, info, frag_pos))
                        for frag_mz, src_set in source_map_pos.items():
                            if frag_mz not in self.frag_source_lookup_pos:
                                self.frag_source_lookup_pos[frag_mz] = set()
                            self.frag_source_lookup_pos[frag_mz].update(src_set)
                except:
                    pass

            if pd.notna(mz_n) and str(mz_n).strip():
                try:
                    base = float(mz_n)
                    for shift in [0] + list(ADDUCTS_NEGATIVE.values()):
                        mz_val = base + shift
                        self.sorted_idx_neg.append((mz_val, info, frag_neg))
                        for frag_mz, src_set in source_map_neg.items():
                            if frag_mz not in self.frag_source_lookup_neg:
                                self.frag_source_lookup_neg[frag_mz] = set()
                            self.frag_source_lookup_neg[frag_mz].update(src_set)
                except:
                    pass

        self.sorted_idx_pos.sort(key=lambda x: x[0])
        self.sorted_idx_neg.sort(key=lambda x: x[0])
        self.mz_values_pos = np.array([x[0] for x in self.sorted_idx_pos])
        self.mz_values_neg = np.array([x[0] for x in self.sorted_idx_neg])

    def binary_search_range(self, mz_array, mz, tolerance_ppm):
        if len(mz_array) == 0:
            return []
        tolerance = mz * tolerance_ppm / 1e6
        mz_min = mz - tolerance
        mz_max = mz + tolerance
        left = bisect_left(mz_array, mz_min)
        right = bisect_right(mz_array, mz_max)
        return list(range(left, min(right, len(mz_array))))

    def search_positive(self, prec_mz, tolerance_ppm):
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
# 多优先级数据库
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
        for herb in self.priority_herbs:
            mask = remaining[self.herb_col].str.contains(herb, na=False, case=False)
            herb_df = remaining[mask].copy()
            remaining = remaining[~mask]
            self.priority_indexes[herb] = EfficientDatabaseIndex(
                herb_df, self.herb_col, True, f'[{herb}]来源'
            )
        self.other_index = EfficientDatabaseIndex(remaining, self.herb_col, False, '其他来源')

    def search(self, prec_mz, mode):
        candidates = []
        for herb, idx in self.priority_indexes.items():
            if mode == 'positive':
                candidates.extend(idx.search_positive(prec_mz, PRIMARY_PPM_TOLERANCE))
            else:
                candidates.extend(idx.search_negative(prec_mz, PRIMARY_PPM_TOLERANCE))
        if mode == 'positive':
            candidates.extend(self.other_index.search_positive(prec_mz, PRIMARY_PPM_TOLERANCE))
        else:
            candidates.extend(self.other_index.search_negative(prec_mz, PRIMARY_PPM_TOLERANCE))
        return candidates

# ============================================================================
# 主鉴定器
# ============================================================================
class UniversalIdentifier:
    def __init__(self, ms_pos_file, ms_neg_file, db_file, priority_herbs, progress_callback=None):
        self.priority_herbs = priority_herbs
        self.progress_callback = progress_callback

        self.ms_pos = load_ms(ms_pos_file) if ms_pos_file else pd.DataFrame()
        self.ms_neg = load_ms(ms_neg_file) if ms_neg_file else pd.DataFrame()
        self.db = load_db(db_file)
        if self.db.empty:
            return

        self.db_idx = MultiPriorityDB(self.db, priority_herbs)
        self.results = {}

    def match_fragments_v12(self, obs_frags, ref_frags, prec_mz, frag_source_lookup):
        if not obs_frags or not ref_frags:
            return [], [], {}

        matched_obs = []
        matched_ref = set()
        matched_sources = {}

        for obs_mz in obs_frags:
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
                    break

        return list(set(matched_obs)), list(matched_ref), matched_sources

    def identify(self):
        pos_data = self._extract_precursors(self.ms_pos, 'positive')
        neg_data = self._extract_precursors(self.ms_neg, 'negative')

        total = len(pos_data) + len(neg_data)
        processed = 0

        for i, prec in enumerate(pos_data):
            self._process_precursor(prec, 'positive')
            processed += 1
            if self.progress_callback:
                self.progress_callback(processed, total)

        for i, prec in enumerate(neg_data):
            self._process_precursor(prec, 'negative')
            processed += 1
            if self.progress_callback:
                self.progress_callback(processed, total)

        self._merge_and_score()
        self._build_fragment_sources()

        return self.results

    def _extract_precursors(self, df, mode):
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
        prec_mz = prec['precursor_mz']
        obs_frags = prec['fragments']
        rt = prec.get('rt')

        candidates = self.db_idx.search(prec_mz, mode)

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

            if mode == 'positive':
                if obs_frags:
                    matched_obs, matched_ref, matched_src = self.match_fragments_v12(
                        obs_frags, cand['ref_frags'], prec_mz, frag_source_lookup
                    )
                    res['matched_frags_pos'].extend(matched_obs)
                    res['matched_sources_pos'].update(matched_src)
                res['obs_frags_pos'].extend(obs_frags)
                res['rt_pos'] = rt
                res['ppm_pos'] = min(res['ppm_pos'], cand['ppm'])
                res['m/z_pos'] = prec_mz
            else:
                if obs_frags:
                    matched_obs, matched_ref, matched_src = self.match_fragments_v12(
                        obs_frags, cand['ref_frags'], prec_mz, frag_source_lookup
                    )
                    res['matched_frags_neg'].extend(matched_obs)
                    res['matched_sources_neg'].update(matched_src)
                res['obs_frags_neg'].extend(obs_frags)
                res['rt_neg'] = rt
                res['ppm_neg'] = min(res['ppm_neg'], cand['ppm'])
                res['m/z_neg'] = prec_mz

    def _merge_and_score(self):
        for key, res in self.results.items():
            all_matched_obs = list(set(res['matched_frags_pos'] + res['matched_frags_neg']))
            all_obs = list(set(res['obs_frags_pos'] + res['obs_frags_neg']))

            all_matched_sources = {}
            all_matched_sources.update(res.get('matched_sources_pos', {}))
            all_matched_sources.update(res.get('matched_sources_neg', {}))

            has_pos = res['m/z_pos'] is not None
            has_neg = res['m/z_neg'] is not None
            has_pos_neg = has_pos and has_neg

            rt_match = False
            if res['rt_pos'] and res['rt_neg']:
                if abs(res['rt_pos'] - res['rt_neg']) <= RT_TOLERANCE:
                    rt_match = True

            if has_pos:
                ppm = res['ppm_pos']
                actual_mz = res['m/z_pos']
                main_mode = 'positive'
            elif has_neg:
                ppm = res['ppm_neg']
                actual_mz = res['m/z_neg']
                main_mode = 'negative'
            else:
                min_ppm = min(res['ppm_pos'], res['ppm_neg'])
                ppm = min_ppm if min_ppm < 999 else 50
                actual_mz = res['m/z_pos'] or res['m/z_neg']
                main_mode = 'positive' if res['m/z_pos'] else 'negative'

            matched_count = len(all_matched_obs)
            all_ref_frags = res.get('ref_frags_pos', []) + res.get('ref_frags_neg', [])
            total_ref_frags = len(all_ref_frags)

            if total_ref_frags > 0:
                coverage_ratio = min(matched_count / total_ref_frags, 1.0)
            else:
                coverage_ratio = 0

            compound_herb = res.get('herb', '')
            is_source_herb_match = False
            if compound_herb and self.priority_herbs:
                for priority_herb in self.priority_herbs:
                    if priority_herb.lower() in compound_herb.lower():
                        is_source_herb_match = True
                        break

            res['matched_frag_count'] = matched_count
            res['total_ref_frags'] = total_ref_frags
            res['coverage_ratio'] = coverage_ratio
            res['has_pos_neg'] = has_pos_neg
            res['rt_match'] = rt_match
            res['main_mode'] = main_mode
            res['ppm'] = ppm
            res['actual_mz'] = actual_mz
            res['has_fragment_data'] = len(all_obs) > 0
            res['is_source_herb_match'] = is_source_herb_match

            confidence = calculate_confidence_v2(
                matched_count, total_ref_frags, res['total_lit_count'],
                has_pos_neg, ppm,
                is_priority=res.get('is_priority', False)
            )
            res['confidence'] = confidence

            level_code, level_name, suggestion = get_confidence_level_v2(
                matched_count, ppm, res['total_lit_count'], is_source_herb_match
            )
            res['level_code'] = level_code
            res['level_name'] = level_name
            res['suggestion'] = suggestion

            res['all_matched_obs'] = all_matched_obs
            res['all_matched_sources'] = all_matched_sources
            res['all_obs_frags'] = all_obs

    def _build_fragment_sources(self):
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

    def export_reports(self, output_prefix='compound', max_confirmed=None):
        sorted_results = sorted(self.results.items(),
                               key=lambda x: (0 if x[1]['is_priority'] else 1, -x[1]['confidence']))

        if max_confirmed is not None:
            confirmed_count = 0
            for key, res in sorted_results:
                if res['level_code'] == 'I':
                    if confirmed_count < max_confirmed:
                        confirmed_count += 1
                    else:
                        res['level_code'] = 'II'
                        res['level_name'] = '高置信级'
                        res['suggestion'] = '超出确证级数量限制，降为高置信级'

        records = []
        for i, (key, res) in enumerate(sorted_results, 1):

            matched_obs_with_sources = []
            for frag_mz in res.get('all_matched_obs', []):
                sources = res.get('matched_sources_pos', {}).get(frag_mz, set()) | \
                         res.get('matched_sources_neg', {}).get(frag_mz, set())
                sources_str = '; '.join(sorted(sources)) if sources else '无来源'
                matched_obs_with_sources.append(f"{frag_mz:.4f}({sources_str})")

            all_ref_frags = res.get('ref_frags_pos', []) + res.get('ref_frags_neg', [])
            all_ref_str = '; '.join([f'{x:.4f}' for x in set(all_ref_frags)][:20])
            obs_str = '; '.join([f'{x:.4f}' for x in res.get('all_obs_frags', [])[:30]])

            coverage_pct = res.get('coverage_ratio', 0) * 100
            coverage_str = f"{coverage_pct:.1f}%"
            is_source_match = res.get('is_source_herb_match', False)

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
                '匹配观测碎片数': res.get('matched_frag_count', 0),
                '参考碎片总数': res.get('total_ref_frags', 0),
                '碎片覆盖率': coverage_str,
                '药材匹配': '是' if is_source_match else '否',
                '主要碎片离子': obs_str if obs_str else all_ref_str,
                '参考碎片离子': all_ref_str,
                '匹配观测碎片(来源)': '; '.join(matched_obs_with_sources[:10]) if matched_obs_with_sources else '无',
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
        return df

# ============================================================================
# Streamlit 应用
# ============================================================================
st.set_page_config(
    page_title="中药化合物智能鉴定平台 v2.0",
    page_icon="data:image/svg+xml,<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 100 100'><text y='.9em' font-size='90'>🌿</text></svg>",
    layout="wide",
    initial_sidebar_state="expanded"
)

st.markdown("""
<style>
    .main-header {
        background: linear-gradient(135deg, #059669 0%, #0891b2 50%, #7c3aed 100%);
        padding: 2rem 2rem;
        border-radius: 20px;
        margin-bottom: 2rem;
        color: white;
        box-shadow: 0 20px 40px rgba(5, 150, 105, 0.3);
    }
    .stat-card {
        background: rgba(255, 255, 255, 0.95);
        border-radius: 16px;
        padding: 1.5rem;
        text-align: center;
        box-shadow: 0 4px 20px rgba(0, 0, 0, 0.08);
    }
    .stat-number {
        font-size: 2.5rem;
        font-weight: 700;
        background: linear-gradient(135deg, #059669, #0891b2);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
    }
</style>
""", unsafe_allow_html=True)

st.markdown("""
<div class="main-header">
    <h1>🌿 中药化合物智能鉴定平台 v2.0</h1>
    <p>新评分体系 | 碎片匹配60分 + 文献来源30分 + ppm误差5分 + 正负离子5分</p>
</div>
""", unsafe_allow_html=True)

# 登录检查
if not check_login():
    st.stop()

logout_button()

# 侧边栏
with st.sidebar:
    st.markdown("### 📊 功能导航")
    page = st.radio("选择功能", [
        "🏠 首页",
        "🔬 开始鉴定",
        "📋 结果分析",
        "📖 使用说明"
    ], index=0)

    st.markdown("---")
    st.markdown("""
    ### v2.0 新评分体系说明

    **基础分（100分）**:
    - 碎片匹配：最高60分
    - 文献来源：最高30分
    - ppm误差：最高5分
    - 正负离子确认：5分

    **优先级加成**: 50分

    **评级标准**:
    - 来源药材：碎片≥2 AND ppm≤50 → 确证级
    - 非来源：碎片≥20 AND 文献≥5 AND ppm≤15 → 确证级
    """)

# 首页
if page == "🏠 首页":
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.markdown("""
        <div class="stat-card">
            <div class="stat-number">48,886</div>
            <div>数据库化合物</div>
        </div>
        """, unsafe_allow_html=True)
    with col2:
        st.markdown("""
        <div class="stat-card">
            <div class="stat-number">300+</div>
            <div>支持药材</div>
        </div>
        """, unsafe_allow_html=True)
    with col3:
        st.markdown("""
        <div class="stat-card">
            <div class="stat-number">5</div>
            <div>置信度级别</div>
        </div>
        """, unsafe_allow_html=True)
    with col4:
        st.markdown("""
        <div class="stat-card">
            <div class="stat-number">v2.0</div>
            <div>最新版本</div>
        </div>
        """, unsafe_allow_html=True)

    st.markdown("---")
    st.markdown("## 🎯 核心特性")

    features = [
        ("🚀", "高效索引", "二分查找算法，处理速度提升10倍"),
        ("🔬", "精准匹配", "碎片离子±0.15Da容差匹配"),
        ("📊", "智能评分", "新评分体系，多维度验证"),
        ("📁", "报告导出", "支持CSV/Excel格式导出"),
    ]

    cols = st.columns(2)
    for i, (icon, title, desc) in enumerate(features):
        with cols[i % 2]:
            st.markdown(f"### {icon} {title}")
            st.markdown(desc)

# 鉴定页面
elif page == "🔬 开始鉴定":
    st.markdown("## 📤 上传质谱数据")

    st.markdown("**请上传正离子和负离子模式质谱数据（Excel格式）**")

    col1, col2 = st.columns(2)
    with col1:
        ms_positive_file = st.file_uploader(
            "📗 正离子模式 (MS+)",
            type=['xlsx', 'xls', 'csv'],
            help="上传正离子模式质谱数据"
        )
    with col2:
        ms_negative_file = st.file_uploader(
            "📘 负离子模式 (MS-)",
            type=['xlsx', 'xls', 'csv'],
            help="上传负离子模式质谱数据"
        )

    st.markdown("---")
    st.markdown("## ⚙️ 优先级设置")

    st.markdown("输入优先级药材名称（多个用逗号分隔），如：`栀子,党参`")

    priority_herbs = st.text_input(
        "优先级药材",
        value="栀子",
        help="数据库中这些药材的化合物会获得50分优先级加成"
    )

    st.markdown("---")
    st.markdown("## 📂 数据库文件")

    db_option = st.radio(
        "选择数据库来源",
        ["使用内置数据库", "上传自定义数据库"],
        help="内置数据库包含48886条化合物，可直接使用"
    )

    use_builtin_db = (db_option == "使用内置数据库")
    db_file = None
    db_path = BUILTIN_DB_PATH

    if not use_builtin_db:
        db_file = st.file_uploader(
            "上传TCM数据库文件 (.csv)",
            type=['csv'],
            help="上传TCM-SM-MS DB.CSV数据库文件"
        )
        if db_file:
            st.success(f"已上传数据库: {db_file.name}")
            with tempfile.NamedTemporaryFile(delete=False, suffix='.csv') as db_tmp:
                db_tmp.write(db_file.getvalue())
                db_path = db_tmp.name
    else:
        st.info(f"✅ 使用内置数据库: TCM-SM-MS DB.CSV")

    # 确证级数量限制
    max_confirmed = st.number_input("确证级数量限制", min_value=0, max_value=10000, value=500, help="确证级最大数量限制")

    # 开始鉴定
    if st.button("🚀 开始化合物鉴定", type="primary", use_container_width=True):
        if not ms_positive_file and not ms_negative_file:
            st.error("请至少上传一个质谱数据文件（正离子或负离子）")
        elif not use_builtin_db and not db_file:
            st.error("请上传数据库文件或选择使用内置数据库")
        else:
            with st.spinner("正在加载数据..."):
                pos_path = neg_path = None

                if ms_positive_file:
                    with tempfile.NamedTemporaryFile(delete=False, suffix='.xlsx') as f:
                        f.write(ms_positive_file.getvalue())
                        pos_path = f.name

                if ms_negative_file:
                    with tempfile.NamedTemporaryFile(delete=False, suffix='.xlsx') as f:
                        f.write(ms_negative_file.getvalue())
                        neg_path = f.name

                herbs = [h.strip() for h in priority_herbs.split(',')]

                progress_bar = st.progress(0)
                status_text = st.empty()

                def update_progress(current, total):
                    progress_bar.progress(int(current / total * 100))
                    status_text.text(f"处理中... {current}/{total}")

                try:
                    identifier = UniversalIdentifier(
                        pos_path, neg_path, db_path, herbs,
                        progress_callback=update_progress
                    )

                    status_text.text("正在鉴定化合物...")
                    results = identifier.identify()

                    status_text.text("正在生成报告...")
                    df = identifier.export_reports('report', max_confirmed=max_confirmed)

                    st.session_state['results_df'] = df
                    st.session_state['identifier'] = identifier

                    progress_bar.progress(100)
                    status_text.text("鉴定完成！")

                    st.success(f"✅ 鉴定完成！共识别出 {len(df)} 个化合物")

                    # 显示统计
                    col1, col2, col3, col4, col5 = st.columns(5)
                    level_counts = df['评级名称'].value_counts()

                    with col1:
                        confirmed = len(df[df['评级'] == 'I'])
                        st.metric("确证级", confirmed)
                    with col2:
                        high_conf = len(df[df['评级'] == 'II'])
                        st.metric("高置信级", high_conf)
                    with col3:
                        probable = len(df[df['评级'] == 'III'])
                        st.metric("推定级", probable)
                    with col4:
                        low_conf = len(df[df['评级'] == 'IV'])
                        st.metric("提示级", low_conf)
                    with col5:
                        excluded = len(df[df['评级'] == 'V'])
                        st.metric("排除级", excluded)

                    st.markdown("---")
                    if st.button("📊 查看详细结果"):
                        st.session_state['current_page'] = "📋 结果分析"
                        st.rerun()

                except Exception as e:
                    st.error(f"鉴定出错: {str(e)}")

# 结果分析页面
elif page == "📋 结果分析":
    st.markdown("## 📊 鉴定结果分析")

    if 'results_df' not in st.session_state:
        st.info("请先进行化合物鉴定")
    else:
        df = st.session_state['results_df']

        # 统计卡片
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("总化合物数", len(df))
        with col2:
            high_count = len(df[df['评级名称'].isin(['确证级', '高置信级'])])
            st.metric("高置信化合物", high_count)
        with col3:
            avg_score = df['置信度'].mean()
            st.metric("平均置信度", f"{avg_score:.1f}")

        st.markdown("---")

        # 评级分布
        st.markdown("### 评级分布")
        level_order = ['确证级', '高置信级', '推定级', '提示级', '排除级']
        level_counts = df['评级名称'].value_counts().reindex(level_order).fillna(0)
        st.bar_chart(level_counts)

        st.markdown("---")

        # 筛选器
        st.markdown("### 🔍 筛选结果")

        col1, col2 = st.columns(2)
        with col1:
            level_filter = st.multiselect(
                "选择评级",
                level_order,
                default=['确证级', '高置信级']
            )
        with col2:
            min_frag = st.slider("最小匹配碎片数", 0, 10, 0)

        if level_filter:
            filtered_df = df[df['评级名称'].isin(level_filter)]
        else:
            filtered_df = df

        if min_frag > 0:
            filtered_df = filtered_df[filtered_df['匹配观测碎片数'] >= min_frag]

        st.markdown(f"筛选结果: {len(filtered_df)} 个化合物")

        # 显示列选择
        all_cols = df.columns.tolist()
        default_cols = ['序号', '化合物中文名', '分子式', 'ppm', '匹配观测碎片数', '评级名称', '置信度', '药材来源']
        selected_cols = st.multiselect("选择显示的列", all_cols, default=[c for c in default_cols if c in all_cols])

        if selected_cols:
            display_df = filtered_df[selected_cols]
        else:
            display_df = filtered_df

        st.dataframe(display_df, use_container_width=True, hide_index=True)

        st.markdown("---")
        st.markdown("### 📥 导出报告")

        col1, col2 = st.columns(2)
        with col1:
            csv = df.to_csv(index=False, encoding='utf-8-sig')
            st.download_button(
                "📥 下载CSV报告",
                csv,
                file_name=f"中药化合物鉴定报告_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
                mime="text/csv"
            )
        with col2:
            buffer = BytesIO()
            with pd.ExcelWriter(buffer, engine='openpyxl') as writer:
                df.to_excel(writer, index=False, sheet_name='完整报告')

                for level in ['I', 'II', 'III', 'IV', 'V']:
                    level_df = df[df['评级'] == level]
                    if not level_df.empty:
                        sheet_name = f"级别{level}"
                        level_df.to_excel(writer, index=False, sheet_name=sheet_name)

            st.download_button(
                "📥 下载Excel报告",
                buffer.getvalue(),
                file_name=f"中药化合物鉴定报告_{datetime.now().strftime('%Y%m%d_%H%M%S')}.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )

# 使用说明页面
elif page == "📖 使用说明":
    st.markdown("## 📖 使用说明")

    st.markdown("""
    ### 数据格式要求

    #### 质谱数据格式
    | 列名 | 说明 |
    |------|------|
    | Precursor M/z | 母离子m/z值 |
    | Peak_1_m/z | 二级碎片离子m/z |
    | 出峰时间t/min | 保留时间（分钟） |

    #### 数据库格式 (CSV)
    | 列名 | 说明 |
    |------|------|
    | 药材名称 | 中药名称 |
    | 名称（中文） | 化合物中文名 |
    | 准分子离子（正） | 正离子模式m/z |
    | 准分子离子（负） | 负离子模式m/z |
    | 碎片离子（正） | 正离子模式碎片 |
    | 碎片离子（负） | 负离子模式碎片 |
    | 文献来源 | 文献信息 |

    ### v2.0 新评分体系

    | 评分项目 | 最高分值 | 说明 |
    |------|---------|------|
    | 碎片匹配 | 60分 | 匹配碎片数量加分 |
    | 文献来源 | 30分 | 文献数量加分 |
    | ppm误差 | 5分 | ppm≤10得5分 |
    | 正负离子确认 | 5分 | 同时确认加5分 |
    | 优先级加成 | 50分 | 优先级药材另外加分 |

    ### 评级标准 (v2.0)

    | 来源类型 | 确证级条件 | 高置信级条件 |
    |---------|-----------|-------------|
    | 来源药材 | ppm≤50 AND 碎片≥2 | ppm≤50 AND 碎片≥1 |
    | 非来源药材 | 碎片≥20 AND 文献≥5 AND ppm≤15 | 碎片≥10 AND 文献≥2 AND ppm≤50 |
    """)
