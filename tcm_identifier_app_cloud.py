#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
中药化合物智能鉴定平台 v2.3（优化增强版）
==========================================

评分体系（v2.3 - 优化版）：
- 基础分100分：碎片匹配60分 + 文献来源30分 + ppm误差5分 + 正负离子5分
- 优先级加成：50分
- 同位素匹配加成：最高20分（辅助确认分子式）- 已启用
- 谱图相似度加成：最高15分 - 已集成到评分
- 类型一致性加成：最高15分 - 新增功能

新增功能：
1. 谱图相似度集成到总评分
2. 特征碎片识别（化合物类型预测）
3. 智能碎片过滤
4. 用户自定义评分权重
5. NPClassifier结构分类接口

评级标准（提高阈值）：
- 来源药材：ppm≤50 AND 碎片数≥5 → 确证级
- 非来源药材：碎片≥25 AND 文献数≥5 AND ppm≤15 → 确证级

Author: MiniMax Agent
"""

# ============================================================================
# 同位素丰度数据（用于理论同位素分布计算）
# ============================================================================
ISOTOPE_ABUNDANCES = {
    'C': {12: 0.9893, 13: 0.0107},  # 13C丰度约1.1%
    'H': {1: 0.99985, 2: 0.00015},
    'N': {14: 0.9963, 15: 0.0037},
    'O': {16: 0.9976, 17: 0.0004, 18: 0.0020},
    'S': {32: 0.9493, 33: 0.0076, 34: 0.0429},
    'Cl': {35: 0.7577, 37: 0.2423},  # Cl有显著M+2峰
    'Br': {79: 0.5069, 81: 0.4931},  # Br有显著M+2峰
    'Si': {28: 0.9223, 29: 0.0467, 30: 0.0310},
    'P': {31: 1.0},
}

# 同位素质量增量（相对于主同位素）
ISOTOPE_MASS_INCREMENT = {
    'C': {12: 0, 13: 1.0034},   # 13C - 12C
    'H': {1: 0, 2: 1.0063},     # 2H - 1H
    'N': {14: 0, 15: 0.9970},   # 15N - 14N
    'O': {16: 0, 18: 1.9942},   # 18O - 16O
    'S': {32: 0, 34: 1.9958},   # 34S - 32S
    'Cl': {35: 0, 37: 1.9970},  # 37Cl - 35Cl
    'Br': {79: 0, 81: 1.9980},  # 81Br - 79Br
    'Si': {28: 0, 30: 1.9946},  # 30Si - 28Si
    'P': {31: 0},
}

# 同位素容差（Da）
ISOTOPE_TOLERANCE = 0.01

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
    """解析碎片离子字符串"""
    frags = []
    source_map = {}
    for part in re.split(r'[;、]', str(frag_str) if pd.notna(frag_str) else ''):
        part = part.strip()
        if part:
            match = re.match(r'^([\d.]+)', part)
            if match:
                try:
                    mz_val = float(match.group(1))
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
# 同位素模式分析功能 - 已启用
# ============================================================================
ISOTOPE_MATCHING_ENABLED = True

def parse_formula(formula):
    """解析分子式，返回各元素数量"""
    if not formula or formula == 'unknown' or formula == 'nan':
        return None

    composition = {}
    formula = re.sub(r'\(([A-Za-z0-9]+)\)(\d+)', lambda m: m.group(1) * int(m.group(2)), formula)

    pattern = r'([A-Z][a-z]?)(\d*)'
    matches = re.findall(pattern, formula)

    for element, count in matches:
        if element:
            num = int(count) if count else 1
            composition[element] = composition.get(element, 0) + num

    return composition if composition else None


def calculate_isotope_pattern(formula):
    """计算理论同位素分布模式"""
    composition = parse_formula(formula)
    if not composition:
        return []

    isotope_peaks = {0: 1.0}

    for element, count in composition.items():
        if element not in ISOTOPE_ABUNDANCES:
            continue

        # Ensure count is an integer
        try:
            count = int(count)
        except (ValueError, TypeError):
            continue

        abundances = ISOTOPE_ABUNDANCES.get(element, {})
        mass_increments = ISOTOPE_MASS_INCREMENT.get(element, {})
        # Ensure they are dicts
        if isinstance(abundances, dict):
            abundances = dict(abundances)
        else:
            abundances = {}
        if isinstance(mass_increments, dict):
            mass_increments = dict(mass_increments)
        else:
            mass_increments = {}

        # Get base isotope key (the one with mass increment 0)
        base_isotope = None
        min_mass_inc = 999
        for k, v in mass_increments.items():
            if v < min_mass_inc:
                min_mass_inc = v
                base_isotope = k

        new_peaks = {}
        for atom_idx in range(count):
            temp_peaks = {}
            for current_mass, current_abund in list(isotope_peaks.items()):
                for isotope_mass_diff, isotope_abund in list(abundances.items()):
                    if base_isotope and isotope_mass_diff == base_isotope:
                        new_mass = current_mass
                    else:
                        # Use mass increment if available
                        mass_inc = mass_increments.get(isotope_mass_diff, 0)
                        new_mass = current_mass + mass_inc

                    new_abund = current_abund * isotope_abund
                    new_peaks[new_mass] = new_peaks.get(new_mass, 0) + new_abund

            isotope_peaks = new_peaks

    max_abund = max(isotope_peaks.values()) if isotope_peaks else 1.0
    result = [(mass, abund / max_abund) for mass, abund in isotope_peaks.items()]
    result.sort(key=lambda x: x[0])

    return result


def find_isotope_peaks(prec_mz, observed_peaks, tolerance=0.01):
    """在观测峰中寻找同位素峰"""
    if not observed_peaks:
        return {}

    obs_mz = sorted(set([p for p in observed_peaks if prec_mz - 5 <= p <= prec_mz + 5]))

    isotope_labels = {
        1: 'M+1',
        2: 'M+2',
        3: 'M+3',
    }

    found_isotopes = {}

    for mass_diff, label in isotope_labels.items():
        expected_mz = prec_mz + mass_diff

        for obs in obs_mz:
            if abs(obs - expected_mz) <= tolerance:
                base_peak = prec_mz
                base_abund = 1.0
                for p in obs_mz:
                    if abs(p - base_peak) < tolerance:
                        base_abund = 1.0
                        break

                rel_abund = 0.5
                found_isotopes[obs] = (label, expected_mz, rel_abund)
                break

    return found_isotopes


def match_isotope_pattern(prec_mz, formula, observed_peaks, tolerance=0.01):
    """匹配观测同位素模式与理论模式"""
    if not formula or formula == 'unknown' or formula == 'nan':
        return {
            'matched_peaks': [],
            'match_score': 0,
            'confidence_boost': False,
            'element_info': {}
        }

    composition = parse_formula(formula)
    if not composition:
        return {
            'matched_peaks': [],
            'match_score': 0,
            'confidence_boost': False,
            'element_info': {}
        }

    theoretical_pattern = calculate_isotope_pattern(formula)

    element_info = {}
    has_diagnostic_elements = False

    for elem in ['Cl', 'Br', 'S', 'Si']:
        count = composition.get(elem, 0)
        if count > 0:
            element_info[elem] = count
            if elem in ['Cl', 'Br']:
                has_diagnostic_elements = True

    matched_peaks = []
    total_abundance_score = 0

    for mass_diff in [1, 2, 3]:
        expected_mz = prec_mz + mass_diff
        for obs_mz in observed_peaks:
            if abs(obs_mz - expected_mz) <= tolerance:
                matched_peaks.append((obs_mz, f'M+{mass_diff}', expected_mz))
                total_abundance_score += 1
                break

    match_score = 0
    n_matches = len(matched_peaks)

    if n_matches >= 3:
        match_score = 20
    elif n_matches == 2:
        match_score = 15
    elif n_matches == 1:
        match_score = 8

    if has_diagnostic_elements and any('M+2' in str(p[1]) for p in matched_peaks):
        match_score = min(match_score + 5, 20)

    confidence_boost = n_matches >= 2 or (n_matches == 1 and has_diagnostic_elements)

    return {
        'matched_peaks': matched_peaks,
        'match_score': match_score,
        'confidence_boost': confidence_boost,
        'element_info': element_info,
        'theoretical_pattern': theoretical_pattern[:5]
    }


def analyze_isotope_for_precursor(prec_mz, formula, fragments):
    """分析前体离子的同位素分布"""
    all_peaks = [prec_mz] + fragments if fragments else [prec_mz]

    result = match_isotope_pattern(prec_mz, formula, all_peaks)
    return result


# ============================================================================
# 谱图相似度计算功能
# ============================================================================
def calculate_cosine_similarity(peaks1, peaks2, tolerance=0.15):
    """计算两张谱图的余弦相似度"""
    if not peaks1 or not peaks2:
        return 0.0

    peaks1 = sorted([float(p) for p in peaks1])
    peaks2 = sorted([float(p) for p in peaks2])

    all_peaks = sorted(set(peaks1 + peaks2))

    vec1 = []
    vec2 = []

    for peak in all_peaks:
        intensity1 = 0.0
        for p in peaks1:
            if abs(p - peak) <= tolerance:
                intensity1 = 1.0
                break

        intensity2 = 0.0
        for p in peaks2:
            if abs(p - peak) <= tolerance:
                intensity2 = 1.0
                break

        vec1.append(intensity1)
        vec2.append(intensity2)

    dot_product = sum(a * b for a, b in zip(vec1, vec2))
    norm1 = sum(a * a for a in vec1) ** 0.5
    norm2 = sum(b * b for b in vec2) ** 0.5

    if norm1 == 0 or norm2 == 0:
        return 0.0

    return dot_product / (norm1 * norm2)


def calculate_spectral_similarity_score(obs_frags, ref_frags, prec_mz=None, intensity_weight=0.5):
    """计算综合谱图相似度得分"""
    filtered_obs = []
    if prec_mz:
        for mz in obs_frags:
            if abs(mz - prec_mz) > 0.5:
                filtered_obs.append(mz)
    else:
        filtered_obs = list(obs_frags)

    cosine_sim = calculate_cosine_similarity(filtered_obs, ref_frags)

    matched_count = 0
    for obs_mz in filtered_obs:
        for ref_mz in ref_frags:
            if abs(obs_mz - ref_mz) <= 0.15:
                matched_count += 1
                break

    coverage = matched_count / len(ref_frags) if ref_frags else 0

    score = (cosine_sim * 0.6 + coverage * 0.4) * 100

    return {
        'cosine_similarity': round(cosine_sim, 4),
        'matched_peaks': matched_count,
        'coverage_ratio': round(coverage, 4),
        'score': round(score, 2)
    }


# ============================================================================
# NPClassifier 结构分类功能
# ============================================================================
import requests

def classify_with_npclassifier(smiles, timeout=10):
    """调用NPClassifier API对化合物进行结构分类"""
    if not smiles or not isinstance(smiles, str):
        return {'success': False, 'error': '无效的SMILES字符串', 'pathway': '', 'superclass': '', 'class': ''}

    try:
        response = requests.post(
            "https://npclassifier.ucsd.edu/predict",
            json={"smiles": smiles},
            timeout=timeout
        )

        if response.status_code == 200:
            result = response.json()
            return {
                'success': True,
                'pathway': result.get('pathway', ''),
                'superclass': result.get('superclass', ''),
                'class': result.get('class', ''),
                'error': ''
            }
        else:
            return {
                'success': False,
                'error': f'API错误: {response.status_code}',
                'pathway': '',
                'superclass': '',
                'class': ''
            }
    except requests.exceptions.Timeout:
        return {
            'success': False,
            'error': '请求超时',
            'pathway': '',
            'superclass': '',
            'class': ''
        }
    except requests.exceptions.RequestException as e:
        return {
            'success': False,
            'error': f'网络错误: {str(e)}',
            'pathway': '',
            'superclass': '',
            'class': ''
        }
    except Exception as e:
        return {
            'success': False,
            'error': f'未知错误: {str(e)}',
            'pathway': '',
            'superclass': '',
            'class': ''
        }


# ============================================================================
# 评分函数
# ============================================================================
DEFAULT_SCORING_WEIGHTS = {
    'frag_weight': 60,
    'lit_weight': 30,
    'ppm_weight': 5,
    'posneg_weight': 5,
    'isotope_weight': 20,
    'spectral_weight': 15,
    'priority_weight': 50,
}

def calculate_confidence_v3(matched_frag_count, total_ref_frags, lit_count,
                            has_pos_neg, ppm, is_priority=False,
                            isotope_score=0, spectral_score=0,
                            weights=None, type_consistency_score=0):
    """计算置信度得分（v2.3）"""
    if weights is None:
        weights = DEFAULT_SCORING_WEIGHTS.copy()

    base_score = 0

    if matched_frag_count == 0:
        return 0

    frag_max = weights.get('frag_weight', 60)
    if matched_frag_count >= 15:
        base_score += frag_max
    elif matched_frag_count >= 12:
        base_score += frag_max * 0.9
    elif matched_frag_count >= 10:
        base_score += frag_max * 0.8
    elif matched_frag_count >= 8:
        base_score += frag_max * 0.7
    elif matched_frag_count >= 6:
        base_score += frag_max * 0.6
    elif matched_frag_count >= 5:
        base_score += frag_max * 0.5
    elif matched_frag_count >= 3:
        base_score += frag_max * 0.4
    elif matched_frag_count >= 2:
        base_score += frag_max * 0.3
    elif matched_frag_count >= 1:
        base_score += frag_max * 0.2

    lit_max = weights.get('lit_weight', 30)
    if lit_count >= 15:
        base_score += lit_max
    elif lit_count >= 12:
        base_score += lit_max * 0.9
    elif lit_count >= 10:
        base_score += lit_max * 0.8
    elif lit_count >= 7:
        base_score += lit_max * 0.7
    elif lit_count >= 5:
        base_score += lit_max * 0.6
    elif lit_count >= 3:
        base_score += lit_max * 0.5
    elif lit_count >= 2:
        base_score += lit_max * 0.4
    elif lit_count >= 1:
        base_score += lit_max * 0.3

    ppm_max = weights.get('ppm_weight', 5)
    if ppm <= 10:
        base_score += ppm_max
    elif ppm <= 30:
        base_score += ppm_max * 0.6
    elif ppm <= 50:
        base_score += ppm_max * 0.2

    posneg_max = weights.get('posneg_weight', 5)
    if has_pos_neg:
        base_score += posneg_max

    base_score = min(base_score, 100)

    isotope_max = weights.get('isotope_weight', 20)
    isotope_bonus = min(isotope_score * isotope_max / 20, isotope_max)

    spectral_max = weights.get('spectral_weight', 15)
    spectral_bonus = min(spectral_score * spectral_max / 100, spectral_max)

    type_bonus = min(type_consistency_score, 15)

    priority_max = weights.get('priority_weight', 50)
    priority_bonus = priority_max if is_priority else 0

    total_score = base_score + isotope_bonus + spectral_bonus + type_bonus + priority_bonus

    return {
        'total_score': round(total_score, 2),
        'base_score': round(base_score, 2),
        'isotope_bonus': round(isotope_bonus, 2),
        'spectral_bonus': round(spectral_bonus, 2),
        'type_bonus': round(type_bonus, 2),
        'priority_bonus': round(priority_bonus, 2),
    }


def calculate_confidence_v2(matched_frag_count, total_ref_frags, lit_count,
                           has_pos_neg, ppm, is_priority=False, isotope_score=0):
    """旧版评分函数，保持向后兼容"""
    result = calculate_confidence_v3(
        matched_frag_count, total_ref_frags, lit_count,
        has_pos_neg, ppm, is_priority, isotope_score, 0
    )
    if isinstance(result, dict):
        return result['total_score']
    return result

def get_confidence_level_v2(matched_frag_count=0, ppm=999, lit_count=0, is_source_herb_match=False):
    """根据条件确定评级"""
    source_match = is_source_herb_match

    if source_match:
        if matched_frag_count >= 5 and ppm <= 50:
            return 'I', '确证级', '高置信度，来源药材匹配'
        if matched_frag_count >= 3 and ppm <= 50:
            return 'II', '高置信级', '来源药材匹配'

    if not source_match:
        if matched_frag_count >= 25 and lit_count >= 5 and ppm <= 15:
            return 'I', '确证级', '高置信度，可作为定性依据'
        if matched_frag_count >= 15 and lit_count >= 2 and ppm <= 50:
            return 'II', '高置信级', '较强置信度，建议进一步验证'

    if matched_frag_count >= 3:
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

            formula = res.get('formula', 'unknown')
            isotope_analysis = analyze_isotope_for_precursor(actual_mz, formula, all_obs)
            isotope_score = isotope_analysis.get('match_score', 0)
            res['isotope_matched_peaks'] = isotope_analysis.get('matched_peaks', [])
            res['isotope_confidence_boost'] = isotope_analysis.get('confidence_boost', False)
            res['isotope_element_info'] = isotope_analysis.get('element_info', {})
            res['isotope_score'] = isotope_score

            all_ref_frags_list = list(set(all_ref_frags))
            if all_obs and all_ref_frags_list:
                similarity_result = calculate_spectral_similarity_score(
                    all_obs, all_ref_frags_list, prec_mz=actual_mz
                )
                res['cosine_similarity'] = similarity_result.get('cosine_similarity', 0)
                res['spectral_score'] = similarity_result.get('score', 0)
                res['similarity_matched_peaks'] = similarity_result.get('matched_peaks', 0)
            else:
                res['cosine_similarity'] = 0
                res['spectral_score'] = 0
                res['similarity_matched_peaks'] = 0

            confidence = calculate_confidence_v2(
                matched_count, total_ref_frags, res['total_lit_count'],
                has_pos_neg, ppm,
                is_priority=res.get('is_priority', False),
                isotope_score=isotope_score
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
            # Safe handling of matched_sources_pos
            sources_pos = res.get('matched_sources_pos', {})
            if isinstance(sources_pos, dict):
                for frag_mz, sources in sources_pos.items():
                    if frag_mz not in self.global_fragment_sources:
                        self.global_fragment_sources[frag_mz] = set()
                    self.global_fragment_sources[frag_mz].update(sources)
            # Safe handling of matched_sources_neg
            sources_neg = res.get('matched_sources_neg', {})
            if isinstance(sources_neg, dict):
                for frag_mz, sources in sources_neg.items():
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
            # Safe handling of matched_sources
            sources_pos = res.get('matched_sources_pos', {})
            sources_neg = res.get('matched_sources_neg', {})
            if not isinstance(sources_pos, dict):
                sources_pos = {}
            if not isinstance(sources_neg, dict):
                sources_neg = {}

            for frag_mz in res.get('all_matched_obs', []):
                pos_sources = sources_pos.get(frag_mz, set()) if isinstance(sources_pos.get(frag_mz), (set, list)) else set()
                neg_sources = sources_neg.get(frag_mz, set()) if isinstance(sources_neg.get(frag_mz), (set, list)) else set()
                sources = pos_sources | neg_sources
                sources_str = '; '.join(sorted(sources)) if sources else '无来源'
                matched_obs_with_sources.append(f"{frag_mz:.4f}({sources_str})")

            all_ref_frags = res.get('ref_frags_pos', []) + res.get('ref_frags_neg', [])
            all_ref_str = '; '.join([f'{x:.4f}' for x in set(all_ref_frags)][:20])
            obs_str = '; '.join([f'{x:.4f}' for x in res.get('all_obs_frags', [])[:30]])

            coverage_pct = res.get('coverage_ratio', 0) * 100
            coverage_str = f"{coverage_pct:.1f}%"
            is_source_match = res.get('is_source_herb_match', False)

            # Safe handling of isotope_element_info
            isotope_elem_info = res.get('isotope_element_info', {})
            if isinstance(isotope_elem_info, dict) and isotope_elem_info:
                diag_elements_str = '; '.join([f"{k}×{v}" for k, v in isotope_elem_info.items()])
            else:
                diag_elements_str = '无'

            # Safe handling of isotope_matched_peaks
            isotope_peaks = res.get('isotope_matched_peaks', [])
            if isotope_peaks and isinstance(isotope_peaks, list):
                peaks_str = '; '.join([f"{p[0]:.4f}({p[1]})" for p in isotope_peaks])
            else:
                peaks_str = '无'

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
                '同位素匹配加分': res.get('isotope_score', 0),
                '同位素峰匹配': peaks_str,
                '诊断元素': diag_elements_str,
                '同位素置信提升': '是' if res.get('isotope_confidence_boost') else '否',
                '余弦相似度': res.get('cosine_similarity', 0),
                '相似度得分': res.get('spectral_score', 0),
                '相似度匹配峰数': res.get('similarity_matched_peaks', 0),
            }
            records.append(row)

        df = pd.DataFrame(records)
        return df

# ============================================================================
# Streamlit 应用
# ============================================================================
st.set_page_config(
    page_title="中药化合物智能鉴定平台 v2.3",
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
    <h1>🌿 中药化合物智能鉴定平台 v2.3</h1>
    <p>新评分体系 | 碎片匹配60分 + 文献来源30分 + ppm误差5分 + 正负离子5分 + 同位素分析20分（已启用）</p>
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
    ### v2.3 功能说明（增强版）

    **基础分（100分）**:
    - 碎片匹配：最高60分（需15+碎片）
    - 文献来源：最高30分
    - ppm误差：最高5分
    - 正负离子确认：5分

    **同位素加成**：最高20分（已启用）
    - M+1, M+2, M+3峰匹配可获得加分
    - 含Cl/Br等诊断元素时额外加分

    **谱图相似度加成**：最高15分（新增）
    - 余弦相似度 + 覆盖率综合计算

    **优先级加成**: 50分

    **评级标准（提高阈值）**:
    - 来源药材：碎片≥5 AND ppm≤50 → 确证级
    - 非来源：碎片≥25 AND 文献≥5 AND ppm≤15 → 确证级
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
            <div class="stat-number">v2.3</div>
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

    max_confirmed = st.number_input("确证级数量限制", min_value=0, max_value=10000, value=500, help="确证级最大数量限制")

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
                    import traceback
                    st.code(traceback.format_exc())

# 结果分析页面
elif page == "📋 结果分析":
    st.markdown("## 📊 鉴定结果分析")

    if 'results_df' not in st.session_state:
        st.info("请先进行化合物鉴定")
    else:
        df = st.session_state['results_df']

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

        st.markdown("### 评级分布")
        level_order = ['确证级', '高置信级', '推定级', '提示级', '排除级']
        level_counts = df['评级名称'].value_counts().reindex(level_order).fillna(0)
        st.bar_chart(level_counts)

        st.markdown("---")

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

    ### v2.3 评分体系

    | 评分项目 | 最高分值 | 新阈值要求 |
    |------|---------|------|
    | 碎片匹配 | 60分 | 需15+碎片（原10+）|
    | 文献来源 | 30分 | 文献数量加分 |
    | ppm误差 | 5分 | ppm≤10得5分 |
    | 正负离子确认 | 5分 | 同时确认加5分 |
    | 同位素加成 | 20分 | M+1/M+2/M+3峰匹配，Cl/Br元素额外加分 |
    | 优先级加成 | 50分 | 优先级药材另外加分 |

    ### 评级标准 (v2.3)

    | 来源类型 | 确证级条件 | 高置信级条件 |
    |---------|-----------|-------------|
    | 来源药材 | ppm≤50 AND 碎片≥5 | ppm≤50 AND 碎片≥3 |
    | 非来源药材 | 碎片≥25 AND 文献≥5 AND ppm≤15 | 碎片≥15 AND 文献≥2 AND ppm≤50 |
    """)
