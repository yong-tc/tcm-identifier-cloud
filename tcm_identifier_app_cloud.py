#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
🌿 中药 LC-MS/MS 化合物通用鉴定平台 v7.0 - Streamlit (v6 框架 + v3opt7 算法)
==========================================================================

v7.0 变更:
- 核心鉴定算法从 v6 替换为 v3opt7 (RT一致性 + RT离散度 + 正负RT配对 + 同位素精细评分 + 双门棁评级)
- 保留 v6 的 streamlit 框架 (登录 / 侧边栏 / UI / 报告 / xlsx 解析)
- 启动: streamlit run streamlit_app_v7.py
- 默认账号: ZY / 513513
"""
import os, sys, re, json, time, pickle, hashlib, tempfile
from io import BytesIO
from datetime import datetime
from bisect import bisect_left, bisect_right
from collections import defaultdict
from typing import List, Dict, Tuple, Optional

import numpy as np
import pandas as pd
try:
    import streamlit as st
    _HAS_ST = True
except ImportError:
    _HAS_ST = False
    # Fallback: 提供 no-op 装饰器和函数，让核心逻辑可以脱离 streamlit 跑
    class _FakeSS:
        def __init__(self): self._d = {}
        def __contains__(self, k): return k in self._d
        def __getitem__(self, k): return self._d[k]
        def __setitem__(self, k, v): self._d[k] = v
        def get(self, k, default=None): return self._d.get(k, default)
        def set(self, **kw): self._d.update(kw)
        def __enter__(self): return self
        def __exit__(self, *a): pass
        def __getattr__(self, k):
            if k in ('get', 'set', 'rerun', 'stop', 'session_state'):
                return getattr(self, k)
            return lambda *a, **kw: None
        def form_submit_button(self, *a, **kw): return False
    class _FakeSt:
        def __init__(self):
            self.session_state = _FakeSS()
            self.sidebar = _FakeSS()
        def __getattr__(self, k):
            if k == 'cache_data': return lambda *a, **kw: (lambda f: f)
            if k == 'cache_resource': return lambda *a, **kw: (lambda f: f)
            if k == 'set_page_config': return lambda *a, **kw: None
            if k == 'form': return lambda *a, **kw: _FakeSS()
            if k == 'spinner': return lambda *a, **kw: _FakeSS()
            if k == 'columns': return lambda spec, **kw: [_FakeSS() for _ in (spec if isinstance(spec, list) else range(spec))]
            if k == 'tabs': return lambda labels, **kw: [_FakeSS() for _ in labels]
            if k == 'selectbox':
                def _sb(*a, **kw):
                    opts = a[1] if len(a) > 1 else kw.get('options', [])
                    idx = kw.get('index', 0)
                    if opts and 0 <= idx < len(opts): return opts[idx]
                    return opts[0] if opts else None
                return _sb
            if k == 'multiselect':
                def _ms(*a, **kw): return kw.get('default', [])
                return _ms
            if k == 'radio':
                def _r(*a, **kw):
                    opts = a[1] if len(a) > 1 else kw.get('options', [])
                    return opts[0] if opts else None
                return _r
            return lambda *a, **kw: None
    st = _FakeSt()

# ============================================================================
# 路径
# ============================================================================
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


def _find_attachments_dir():
    """找 TCM-SM-MS DB.csv 所在目录"""
    env_dir = os.environ.get('TCM_ATTACH_DIR')
    if env_dir and os.path.isdir(env_dir): return env_dir
    candidates = [
        '/workspace/attachments',
        '/app/attachments',
        os.path.join(_SCRIPT_DIR, '..', '..', 'attachments'),
        os.path.join(_SCRIPT_DIR, '..', 'attachments'),
        os.path.join(_SCRIPT_DIR, 'attachments'),
        os.path.expanduser('~/attachments'),
        os.path.join('/tmp', 'attachments'),
    ]
    for c in candidates:
        if os.path.exists(c) and os.path.isdir(c):
            return c
    return _SCRIPT_DIR


ATTACH_DIR = _find_attachments_dir()
DB_CSV = None
for _n in ['TCM-SM-MS DB.csv', 'TCM-SM-MS DB.CSV']:
    _p = os.path.join(ATTACH_DIR, _n)
    if os.path.exists(_p):
        DB_CSV = _p; break
if DB_CSV is None:
    import glob as _g
    _csvs = _g.glob(os.path.join(ATTACH_DIR, '*.csv'))
    if _csvs: DB_CSV = _csvs[0]
if DB_CSV is None:
    st.error(f'❌ 找不到 TCM-SM-MS DB.csv, 请放到 `{ATTACH_DIR}/` 下')
    st.stop()

DB_PKL = os.path.join(ATTACH_DIR, '_db_index.pkl')
SAMPLE_CACHE_DIR = os.path.join(ATTACH_DIR, '_sample_cache')
os.makedirs(SAMPLE_CACHE_DIR, exist_ok=True)

REPORTS_DIR = os.path.join(_SCRIPT_DIR, 'reports_streamlit')
os.makedirs(REPORTS_DIR, exist_ok=True)

# ============================================================================
# 配置
# ============================================================================
VALID_USER, VALID_PASS = 'ZY', '513513'
INSTRUMENT_PPM = {'orbitrap': 10, 'qtof': 15, 'qqq': 30, 'low_res': 50}
SECONDARY_DA_TOLERANCE = 0.15
# R3 原版三个主推药材 (priority 辅以避免 L1 丢分)
SUPPORTED_HERBS = ['栀子', '党参', '西红花']

ISOTOPE_ABUNDANCES = {
    'C': {12: 0.9893, 13: 0.0107}, 'H': {1: 0.99985, 2: 0.00015},
    'N': {14: 0.9963, 15: 0.0037}, 'O': {16: 0.9976, 17: 0.0004, 18: 0.0020},
    'S': {32: 0.9493, 33: 0.0076, 34: 0.0429}, 'Cl': {35: 0.7577, 37: 0.2423},
    'Br': {79: 0.5069, 81: 0.4931}, 'Si': {28: 0.9223, 29: 0.0467, 30: 0.0310},
    'P': {31: 1.0},
}
ISOTOPE_MASS_INCREMENT = {
    'C': {12: 0, 13: 1.0034}, 'H': {1: 0, 2: 1.0063}, 'N': {14: 0, 15: 0.9970},
    'O': {16: 0, 18: 1.9942}, 'S': {32: 0, 34: 1.9958}, 'Cl': {35: 0, 37: 1.9970},
    'Br': {79: 0, 81: 1.9980}, 'Si': {28: 0, 30: 1.9942}, 'P': {31: 0},
}

# ============================================================================
# 工具函数
# ============================================================================
def safe_float(v):
    if v is None: return None
    if isinstance(v, float):
        if np.isnan(v): return None
        return v
    s = str(v).strip()
    if not s or s == 'nan': return None
    if '/' in s: s = s.split('/')[0].strip()
    if ';' in s: s = s.split(';')[0].strip()
    if '-' in s and not s.startswith('-'):
        parts = s.split('-')
        if len(parts) == 2 and '.' in parts[0] and '.' in parts[1]: return None
        s = parts[0].strip()
    try: return float(s)
    except: return None


def parse_fragments(frag_str):
    frags = []
    if not frag_str: return frags
    for part in re.split(r'[;、/]', str(frag_str)):
        part = part.strip()
        if part:
            m = re.match(r'^([\d.]+)', part)
            if m:
                try: frags.append(float(m.group(1)))
                except: pass
    return frags


def parse_formula(formula):
    if not formula or formula in ('unknown', 'nan'): return None
    composition = {}
    formula = re.sub(r'\(([A-Za-z0-9]+)\)(\d+)', lambda m: m.group(1) * int(m.group(2)), str(formula))
    for element, count in re.findall(r'([A-Z][a-z]?)(\d*)', formula):
        if element:
            composition[element] = composition.get(element, 0) + (int(count) if count else 1)
    return composition if composition else None


_UNICODE_SUB = str.maketrans('₀₁₂₃₄₅₆₇₈₉', '0123456789')
def normalize_formula(formula):
    if not formula or formula in ('unknown', 'nan'): return ''
    return str(formula).strip().translate(_UNICODE_SUB)


def cosine_intensity(peaks1, w1, peaks2, w2, tol=0.05):
    if not peaks1 or not peaks2: return 0.0
    p1 = np.array(peaks1, dtype=float)
    p2 = np.array(peaks2, dtype=float)
    weights1 = np.array(w1, dtype=float) if w1 else np.ones(len(p1))
    weights2 = np.array(w2, dtype=float) if w2 else np.ones(len(p2))
    w1_max = weights1.max() if len(weights1) else 1.0
    w2_max = weights2.max() if len(weights2) else 1.0
    if w1_max > 0: weights1 = weights1 / w1_max
    if w2_max > 0: weights2 = weights2 / w2_max
    all_mz = np.unique(np.concatenate([p1, p2]))
    v1 = np.zeros(len(all_mz)); v2 = np.zeros(len(all_mz))
    if len(p1) > 0:
        idx1 = np.searchsorted(all_mz, p1); idx1 = np.clip(idx1, 0, len(all_mz) - 1)
        for j in range(len(p1)):
            if abs(all_mz[idx1[j]] - p1[j]) <= tol and weights1[j] > v1[idx1[j]]:
                v1[idx1[j]] = weights1[j]
    if len(p2) > 0:
        idx2 = np.searchsorted(all_mz, p2); idx2 = np.clip(idx2, 0, len(all_mz) - 1)
        for j in range(len(p2)):
            if abs(all_mz[idx2[j]] - p2[j]) <= tol and weights2[j] > v2[idx2[j]]:
                v2[idx2[j]] = weights2[j]
    n1, n2 = np.linalg.norm(v1), np.linalg.norm(v2)
    if n1 == 0 or n2 == 0: return 0.0
    return float(np.dot(v1, v2) / (n1 * n2))


# ============================================================================
# DB 预索引 (缓存到 session_state, 跨 rerun 复用)
# ============================================================================
@st.cache_data(show_spinner='加载 TCM-SM-MS DB 索引...')
def build_db_pkl_cached(db_csv=DB_CSV, db_pkl=DB_PKL):
    if os.path.exists(db_pkl):
        with open(db_pkl, 'rb') as f:
            return pickle.load(f)
    df = pd.read_csv(db_csv, encoding='utf-8-sig', low_memory=False)
    col_map = {
        '药材名': '药材名称', '中文名': '名称（中文）',
        '英文名': '名称（英文）', '文献': '文献来源',
        '准分子离子(正离子模式)': '准分子离子（正）',
        '准分子离子(负离子模式)': '准分子离子（负）',
        '碎片离子(正离子模式)': '碎片离子（正）',
        '碎片离子(负离子模式)': '碎片离子（负）',
    }
    df = df.rename(columns=col_map)
    herb_col = '药材名称' if '药材名称' in df.columns else '药材名'
    mz_p = df['准分子离子（正）'].map(safe_float).fillna(0).values
    mz_n = df['准分子离子（负）'].map(safe_float).fillna(0).values
    name_cn = df.get('名称（中文）', pd.Series(['']*len(df))).fillna('').astype(str).values
    name_en = df.get('名称（英文）', pd.Series(['']*len(df))).fillna('').astype(str).values
    formula = df.get('分子式', pd.Series(['']*len(df))).fillna('').astype(str).values
    herb = df[herb_col].fillna('').astype(str).values
    ctype = df.get('化合物类型', pd.Series(['']*len(df))).fillna('').astype(str).values
    cas = df.get('CAS', pd.Series(['']*len(df))).fillna('').astype(str).values
    lit_raw = df.get('文献来源', pd.Series(['']*len(df))).fillna('').astype(str).values
    frag_p_raw = df['碎片离子（正）'].fillna('').astype(str).values
    frag_n_raw = df['碎片离子（负）'].fillna('').astype(str).values
    lit_count = np.array([len([s for s in x.split(';') if s.strip()]) for x in lit_raw])
    frag_p_cache = [parse_fragments(x) for x in frag_p_raw]
    frag_n_cache = [parse_fragments(x) for x in frag_n_raw]
    names = [(cn if cn and cn != 'nan' else (en if en and en != 'nan' else 'Unknown'))
             for cn, en in zip(name_cn, name_en)]
    n_rows = len(df)
    records = []
    for i in range(n_rows):
        records.append({
            'mz_p': float(mz_p[i]), 'mz_n': float(mz_n[i]),
            'frag_p': frag_p_cache[i], 'frag_n': frag_n_cache[i],
            'name': names[i], 'formula': str(formula[i]),
            'herb': str(herb[i]), 'ctype': str(ctype[i]),
            'cas': str(cas[i]), 'lit': int(lit_count[i]),
        })
    sorted_p = sorted([(r['mz_p'], i) for i, r in enumerate(records) if r['mz_p'] > 0], key=lambda x: x[0])
    sorted_n = sorted([(r['mz_n'], i) for i, r in enumerate(records) if r['mz_n'] > 0], key=lambda x: x[0])
    mz_p_arr = np.array([x[0] for x in sorted_p]) if sorted_p else np.array([])
    mz_n_arr = np.array([x[0] for x in sorted_n]) if sorted_n else np.array([])
    compound_lit_count = df.groupby('名称（中文）').size().to_dict()
    obj = {
        'records': records, 'sorted_p': sorted_p, 'sorted_n': sorted_n,
        'mz_p_arr': mz_p_arr, 'mz_n_arr': mz_n_arr,
        'compound_lit_count': compound_lit_count, 'n_rows': n_rows,
    }
    with open(db_pkl, 'wb') as f:
        pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)
    return obj


def get_all_herbs():
    """从 DB 抽取所有药材 (清洗后)"""
    db = build_db_pkl_cached()
    herbs = set()
    for r in db['records']:
        h = r['herb'].strip()
        if not h: continue
        if re.search(r'\d', h): continue
        if re.search(r'[A-Za-z]{10,}', h): continue  # 全英文长串
        if not re.search(r'[\u4e00-\u9fff]', h): continue
        herbs.add(h)
    return sorted(herbs)


# ============================================================================
# xlsx 解析 (支持 xlsx / csv, 自动识别列)
# ============================================================================
def extract_xlsx_to_list(file_path_or_buffer, is_buffer=False):
    """读 xlsx/csv, 返回 list of dict: {prec, frags, weights, rt}"""
    if is_buffer:
        df = pd.read_excel(file_path_or_buffer, engine='openpyxl')
    else:
        df = pd.read_excel(file_path_or_buffer, engine='openpyxl')
    if df.empty: return []
    # 找列 (容错: 表头可能略不同)
    col_map_xlsx = {
        'prec': ['Precursor M/z', 'Precursor m/z', 'precursor m/z', 'PrecursorMZ'],
        'rt': ['出峰时间t/min', 'Retention time', 'RT (min)', 'RT'],
    }
    def find_col(candidates):
        for c in candidates:
            if c in df.columns: return c
        return None
    prec_col = find_col(col_map_xlsx['prec'])
    rt_col = find_col(col_map_xlsx['rt'])
    if not prec_col:
        return []
    frag_cols = sorted([c for c in df.columns if 'Peak_' in c and '_m/z' in c],
                       key=lambda c: int(re.search(r'Peak_(\d+)_', c).group(1)) if re.search(r'Peak_(\d+)_', c) else 0)
    int_cols = sorted([c for c in df.columns if 'Peak_' in c and '_Intensity' in c],
                      key=lambda c: int(re.search(r'Peak_(\d+)_', c).group(1)) if re.search(r'Peak_(\d+)_', c) else 0)
    # 向量化
    prec_arr = df[prec_col].map(safe_float).values
    rt_arr = df[rt_col].map(safe_float).values if rt_col else np.array([None]*len(df))
    if frag_cols:
        frag_mat = df[frag_cols].map(safe_float).values
        if int_cols:
            int_mat = df[int_cols].map(safe_float).fillna(1.0).values
        else:
            int_mat = np.ones_like(frag_mat)
    else:
        return []
    out = []
    for i in range(len(df)):
        p = prec_arr[i]
        if p is None or p <= 0 or (isinstance(p, float) and np.isnan(p)): continue
        frags = []; weights = []
        row_f = frag_mat[i] if i < frag_mat.shape[0] else []
        row_i = int_mat[i] if i < int_mat.shape[0] else []
        for j, v in enumerate(row_f):
            if v is None or (isinstance(v, float) and np.isnan(v)) or v <= 0: continue
            frags.append(float(v))
            w = row_i[j] if j < len(row_i) else 1.0
            if w is None or (isinstance(w, float) and np.isnan(w)) or w <= 0: w = 1.0
            weights.append(float(w))
        rt = rt_arr[i] if i < len(rt_arr) else None
        if isinstance(rt, float) and np.isnan(rt): rt = None
        out.append({'prec': float(p), 'frags': frags, 'weights': weights, 'rt': rt})
    return out


# ============================================================================
# 鉴定器 (通用版: priority_herbs 是字符串列表)
# ============================================================================
# ============================================================================
# v7.0: 核心鉴定器改为 v3opt7 (RT一致性 + RT离散度 + 正负RT配对 + 同位素精细评分)
# ============================================================================
# v3opt7 引擎定位策略 (按优先级搜索, 适应 Streamlit Cloud / 本地 / 其他部署):
#   1. v7 同目录下的 v3opt7_engine.py (推荐: 打包部署时一并上传)
#   2. 当前工作目录下的 v3opt7_engine.py
#   3. ATTACH_DIR 下的 v3opt7_engine.py
#   4. /workspace/attachments/下的原始文件
#   5. Streamlit Cloud 常见路径: /mount/src/<repo>/...
import importlib.util as _importlib_util


def _find_v3opt7_engine():
    """多路径搜索 v3opt7_engine.py (容错处理,适配 Streamlit Cloud 沙箱)"""
    here = os.path.dirname(os.path.abspath(__file__))
    candidates = [
        os.path.join(here, 'v3opt7_engine.py'),                      # 1. v7 同目录
        os.path.join(os.getcwd(), 'v3opt7_engine.py'),               # 2. 当前工作目录
        os.path.join(ATTACH_DIR, 'v3opt7_engine.py') if ATTACH_DIR else None,  # 3. ATTACH_DIR
        '/workspace/attachments/73a297ee__58e4ea42-04bf-46d8-ae15-e08da4863fc0.py',  # 4. 本地原始路径
        os.path.join(here, '73a297ee__58e4ea42-04bf-46d8-ae15-e08da4863fc0.py'),    # 5. v7 同目录原始名
    ]
    candidates = [c for c in candidates if c]

    # 6. Streamlit Cloud 常见路径: /mount/src/<repo>/  (包裹异常以适应沙箱)
    for root in ['/mount/src', '/app', '/home/app', '/srv', '/opt']:
        try:
            if not os.path.isdir(root):
                continue
            for sub in os.listdir(root):
                sub_path = os.path.join(root, sub)
                if not os.path.isdir(sub_path):
                    continue
                for fn in ['v3opt7_engine.py', '73a297ee__58e4ea42-04bf-46d8-ae15-e08da4863fc0.py']:
                    cand = os.path.join(sub_path, fn)
                    if os.path.exists(cand):
                        candidates.append(cand)
        except (PermissionError, OSError):
            # 沙箱环境下某些路径不可读,直接跳过
            continue

    for cand in candidates:
        if os.path.exists(cand):
            return cand
    return None


@st.cache_resource(show_spinner='加载 v3opt7 鉴定引擎...')
def _load_v3opt7_engine():
    """加载并缓存 v3opt7 引擎 (streamlit rerun 时不会重复 exec_module)"""
    p = _find_v3opt7_engine()
    if p is None:
        raise FileNotFoundError(
            '找不到 v3opt7 引擎文件。\n'
            '请将 v3opt7_engine.py 上传到 v7 应用同目录或 ATTACH_DIR。\n'
            f'已搜索路径示例:\n  {os.path.dirname(os.path.abspath(__file__))}/\n  {os.getcwd()}/\n  {ATTACH_DIR}/'
        )
    spec = _importlib_util.spec_from_file_location('v3opt7_engine', p)
    mod = _importlib_util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


v3opt7_engine = _load_v3opt7_engine()



# 旧 v6 的 _compute_isotope_pattern 和 determine_msi_level 已废弃, 由 v3opt7 内部接管
# 保留同名空壳防止 import 报错


def _compute_isotope_pattern(*args, **kwargs):
    """v7.0 已废弃: 请使用 v3opt7_engine 的同位素评分"""
    return [(0, 1.0)]


def determine_msi_level(*args, **kwargs):
    """v7.0 已废弃: 请使用 v3opt7_engine 的评级"""
    return 4


def identify(extracted_pos, extracted_neg, priority_herbs, instrument='qtof',
             use_isotope=True, add_supported_herbs=False):
    """v7.0 鉴定: 委托给 v3opt7 引擎,输出 schema 跟 v6 兼容

    参数:
        priority_herbs: 用户选的药材列表
        add_supported_herbs: 是否合并 SUPPORTED_HERBS (栀子/党参/西红花)
            - False (默认): 通用模式,只以用户选的药材为 priority
            - True: R3 兼容模式,优先级 = 用户选 + SUPPORTED_HERBS
    """
    db = build_db_pkl_cached()
    ppm_tol = INSTRUMENT_PPM.get(instrument, 15)

    # v7.0 通用版: 默认不加 SUPPORTED_HERBS, 保证 priority 完全由用户控制
    if priority_herbs:
        if add_supported_herbs:
            priority_herbs_eff = list(set(priority_herbs) | set(SUPPORTED_HERBS))
        else:
            priority_herbs_eff = list(priority_herbs)
    else:
        priority_herbs_eff = []

    df = v3opt7_engine.identify_v3opt7(
        extracted_pos or [],
        extracted_neg or [],
        priority_herbs_eff,
        db,
        ppm_tol=ppm_tol,
        da_tol=0.2,
        use_isotope=use_isotope,
        relax_priority=True,
        rt_tol=0.5,
    )
    if df is None or len(df) == 0:
        return pd.DataFrame()

    # ---- 字段重命名 / 补充 v6 UI 期望的列 ----
    # v3opt7 的 '匹配参考碎片数' = v6 语义中的 mc (匹配上的 ref 碎片数)
    # 如果 v3opt7 原始输出里就有 '匹配观测碎片数' (unique obs 命中), 先 drop, 再重命名
    if '匹配观测碎片数' in df.columns and '匹配参考碎片数' in df.columns:
        df = df.drop(columns=['匹配观测碎片数'])
    df = df.rename(columns={'匹配参考碎片数': '匹配观测碎片数'})

    # v6 UI 期望的 '余弦相似度' / '诊断元素' 列
    if '余弦相似度' not in df.columns:
        if '谱图相似度' in df.columns:
            df['余弦相似度'] = df['谱图相似度'].round(4)
        else:
            df['余弦相似度'] = 0.0

    if '诊断元素' not in df.columns:
        def _diag_elements(f):
            comp = parse_formula(f) or {}
            parts = [f'{e}x{comp[e]}' for e in ['Cl', 'Br', 'S', 'Si'] if comp.get(e, 0) > 0]
            return '; '.join(parts) if parts else '无'
        df['诊断元素'] = df['分子式'].apply(_diag_elements)

    # 按 v6 的列顺序重新组织
    preferred_cols = [
        '序号', '化合物中文名', '分子式', 'CAS号',
        'm/z实际值', 'ppm',
        '匹配观测碎片数', '参考碎片总数',
        '碎片覆盖率', '出峰时间t/min',
        '药材匹配',
        '主要碎片离子', '参考碎片离子',
        '库中文献数',
        '评级', '评级名称', 'MSI_Level', '置信度',
        '余弦相似度', '同位素得分', '诊断元素',
        '药材来源', '化合物类型',
    ]
    extra_cols = [c for c in df.columns if c not in preferred_cols]
    df = df[[c for c in preferred_cols if c in df.columns] + extra_cols]
    return df

# ============================================================================
# 报告生成
# ============================================================================
def generate_report_md(df, herb_label, instrument):
    ts = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    rating = df['评级'].value_counts().to_dict() if '评级' in df.columns else {}
    msi = df['MSI_Level'].value_counts().to_dict() if 'MSI_Level' in df.columns else {}
    md = [f'# {herb_label} - LC-MS/MS 化合物鉴定报告', '']
    md.append(f'**生成时间**: {ts}  ')
    md.append(f'**仪器**: {instrument} (PPM {INSTRUMENT_PPM[instrument]})  ')
    md.append(f'**平台**: 中药化合物智能鉴定平台 v6.0 (Streamlit)')
    md.append('')
    md.append('---')
    md.append('')
    md.append('## 📊 总体统计')
    md.append('')
    md.append(f'- 化合物数: **{len(df)}**')
    md.append(f'- 唯一化合物: **{df["化合物中文名"].nunique() if "化合物中文名" in df.columns else 0}**')
    md.append(f'- 平均置信度: **{df["置信度"].mean():.2f}**')
    md.append('')
    md.append('## 🎯 评级分布')
    md.append('')
    md.append('| 评级 | 数量 | 百分比 |')
    md.append('|------|------|--------|')
    total = len(df)
    for k, name in [('I', '确证级'), ('II', '高置信级'), ('III', '推定级'), ('IV', '提示级'), ('V', '排除级')]:
        v = int(rating.get(k, 0))
        pct = v / total * 100 if total > 0 else 0
        md.append(f'| {k} {name} | {v} | {pct:.1f}% |')
    md.append('')
    md.append('## 🧬 MSI 等级分布')
    md.append('')
    md.append('| MSI | 数量 |')
    md.append('|-----|------|')
    for i in [1, 2, 3, 4]:
        md.append(f'| Level {i} | {int(msi.get(i, 0))} |')
    md.append('')
    l1 = df[df['评级'] == 'I'] if '评级' in df.columns else pd.DataFrame()
    if len(l1) > 0:
        md.append(f'## ⭐ Level 1 确证化合物 (前 30)')
        md.append('')
        md.append('| 化合物 | 分子式 | m/z | ppm | 匹配碎片 | 置信度 | MSI |')
        md.append('|--------|--------|-----|-----|----------|--------|-----|')
        for _, r in l1.head(30).iterrows():
            md.append(f"| {r['化合物中文名']} | {r['分子式']} | {r['m/z实际值']} | {r['ppm']} | {r['匹配观测碎片数']} | {r['置信度']} | {r['MSI_Level']} |")
    return '\n'.join(md)


# ============================================================================
# Streamlit UI
# ============================================================================
st.set_page_config(
    page_title='中药 LC-MS/MS 通用鉴定平台 v7.0 (v3opt7 算法)',
    page_icon='🌿',
    layout='wide',
    initial_sidebar_state='expanded',
)

# CSS 美化
st.markdown('''
<style>
.main-title { background: linear-gradient(135deg, #059669 0%, #0891b2 50%, #7c3aed 100%);
              padding: 1.5rem; border-radius: 16px; color: white; margin-bottom: 1.5rem; }
</style>
''', unsafe_allow_html=True)

# 登录
if 'logged_in' not in st.session_state:
    st.session_state.logged_in = False
if not st.session_state.logged_in:
    st.markdown('<div class="main-title"><h1>🌿 中药 LC-MS/MS 化合物通用鉴定平台 v7.0</h1>'
                '<p>基于 TCM-SM-MS DB (48,886 化合物 / 579 药材) 的 LC-MS/MS 智能鉴定 | 算法: v3opt7</p></div>',
                unsafe_allow_html=True)
    col1, col2, col3 = st.columns([1, 2, 1]) if False else st.columns(3)
    with col2:
        with st.form('login'):
            st.markdown('### 🔐 登录 (v7.0 / v3opt7 算法)')
            u = st.text_input('用户名', value='ZY')
            p = st.text_input('密码', value='513513', type='password')
            if st.form_submit_button('登录', use_container_width=True, type='primary'):
                if u == VALID_USER and p == VALID_PASS:
                    st.session_state.logged_in = True
                    st.rerun()
                else:
                    st.error('❌ 账号或密码错误')
    st.stop()

# 侧边栏
with st.sidebar:
    st.markdown('### ⚙️ 鉴定设置')
    instrument = st.selectbox('仪器类型', list(INSTRUMENT_PPM.keys()),
                                index=list(INSTRUMENT_PPM.keys()).index('qtof'),
                                help='不同仪器 PPM 容差不同')
    use_isotope = st.checkbox('开启同位素分布匹配', value=True,
                                help='匹配 Cl/Br/S 等诊断元素的同位素分布')
    st.markdown('---')
    st.markdown('### 🌿 药材选择 (可多选)')
    all_herbs = get_all_herbs()
    herb_mode = st.radio('药材匹配模式', ['单药材', '多药材', '不指定(全库)'],
                          help='不指定时, 所有化合物都按"非 priority"处理')
    selected_herbs = []
    if herb_mode == '单药材':
        # v7.0 通用版: 默认选第一个药材,不再硬编码“栀子”
        h = st.selectbox('选择药材', all_herbs, index=0)
        selected_herbs = [h]
    elif herb_mode == '多药材':
        # v7.0 通用版: 默认选前 2 个药材,不再硬编码“栀子、党参”
        h = st.multiselect('选择药材 (可多选)', all_herbs, default=all_herbs[:2])
        selected_herbs = h
    else:
        selected_herbs = []
        st.caption('当前为全库模式, 所有化合物不享受 priority 加分')
    # v7.0: 是否追加主推三药材 (R3 兼容选项, 默认关闭)
    add_supported = st.checkbox(
        '追加主推三药材 priority (栀子/党参/西红花)',
        value=False,
        help='勾选后,priority 集合 = 用户选 + 栀子/党参/西红花。\n不勾选: priority 完全由用户控制 (通用版)。',
    )
    st.markdown('---')
    st.markdown('### 📈 平台信息')
    st.metric('TCM-SM-MS DB', f'{len(all_herbs)} 药材')
    st.metric('支持仪器', f'{len(INSTRUMENT_PPM)} 种')
    st.metric('当前 PPM', INSTRUMENT_PPM[instrument])
    if st.button('登出', use_container_width=True):
        st.session_state.logged_in = False
        st.rerun()

# 主区
st.markdown('<div class="main-title"><h1>🌿 中药 LC-MS/MS 化合物通用鉴定平台 v7.0</h1>'
            '<p>v6 框架 + v3opt7 核心算法 | 自动适配 DB 内 579 种药材 | MSI Level 1-4 鉴定</p></div>',
            unsafe_allow_html=True)

# 标签页
tab_upload, tab_results, tab_about = st.tabs(['🔬 上传鉴定', '📊 结果分析', '📖 平台说明'])

with tab_upload:
    st.markdown('### 📤 上传 LC-MS/MS 数据')
    st.info('💡 至少上传正离子或负离子模式中的一个文件 (.xlsx / .csv)。\n\n'
            '**期望列**: `Precursor M/z`, `出峰时间t/min`, `Peak_1_m/z` ~ `Peak_N_m/z`, `Peak_N_Intensity`')
    col1, col2 = st.columns(2)
    with col1:
        pos_file = st.file_uploader('正离子模式 (可选)', type=['xlsx', 'csv'],
                                      key='pos', help='ESI+ MS/MS 数据')
    with col2:
        neg_file = st.file_uploader('负离子模式 (可选)', type=['xlsx', 'csv'],
                                      key='neg', help='ESI- MS/MS 数据')
    col_btn1, col_btn2, col_btn3 = st.columns([1, 1, 3])
    run_btn = col_btn1.button('🚀 开始鉴定', type='primary', use_container_width=True)
    example_btn = col_btn2.button('📂 用示例数据', use_container_width=True)
    if example_btn:
        st.session_state['use_example'] = True

    if run_btn or st.session_state.get('use_example'):
        st.session_state['use_example'] = False
        pos_data, neg_data = None, None
        try:
            with st.spinner('📥 解析上传文件...'):
                if pos_file:
                    pos_data = extract_xlsx_to_list(pos_file)
                elif st.session_state.get('use_example_pos'):
                    pass
            if not pos_file and not neg_file and not example_btn:
                st.error('❌ 请至少上传一个文件')
                st.stop()
            if not pos_data and not neg_data and not example_btn:
                st.warning('⚠️ 解析后无数据, 请检查文件格式')
                st.stop()
            with st.spinner(f'🔍 用 {INSTRUMENT_PPM[instrument]} ppm 容差鉴定 {len(selected_herbs) if selected_herbs else "全库"} 药材...'):
                t0 = time.time()
                df = identify(pos_data or [], neg_data or [], selected_herbs,
                              instrument=instrument, use_isotope=use_isotope,
                              add_supported_herbs=add_supported)
                elapsed = time.time() - t0
            if len(df) == 0:
                st.warning('⚠️ 未鉴定到任何化合物。试试调高 PPM 容差或换仪器类型。')
            else:
                st.session_state['result_df'] = df
                st.session_state['herb_label'] = ' / '.join(selected_herbs) if selected_herbs else '全库'
                st.session_state['elapsed'] = elapsed
                st.session_state['instrument'] = instrument
                st.success(f'✅ 鉴定完成! 用时 **{elapsed:.2f}s**, 共 **{len(df)}** 个唯一化合物')
                rating = df['评级'].value_counts().to_dict()
                cols = st.columns(5)
                for i, k in enumerate(['I', 'II', 'III', 'IV', 'V']):
                    cols[i].metric(f'评级 {k}', int(rating.get(k, 0)))
                st.dataframe(df.head(30), use_container_width=True, hide_index=True, height=400)
        except Exception as e:
            st.error(f'❌ 鉴定失败: {e}')
            import traceback; st.code(traceback.format_exc())

with tab_results:
    if 'result_df' not in st.session_state:
        st.info('👈 请先到「上传鉴定」标签页跑数据')
    else:
        df = st.session_state['result_df']
        herb_label = st.session_state.get('herb_label', 'unknown')
        elapsed = st.session_state.get('elapsed', 0)
        instrument = st.session_state.get('instrument', 'qtof')
        st.markdown(f'### 📊 鉴定结果 ({herb_label})')
        c1, c2, c3, c4 = st.columns(4)
        c1.metric('化合物数', len(df))
        c2.metric('平均置信度', f"{df['置信度'].mean():.2f}")
        c3.metric('用时', f'{elapsed:.2f}s')
        c4.metric('药材匹配数', int((df['药材匹配'] == '是').sum()))
        col1, col2 = st.columns(2)
        with col1:
            st.markdown('#### 🎯 评级分布')
            rating = df['评级'].value_counts().reindex(['I', 'II', 'III', 'IV', 'V'], fill_value=0)
            st.bar_chart(rating)
        with col2:
            st.markdown('#### 🧬 MSI 等级')
            msi = df['MSI_Level'].value_counts().sort_index()
            st.bar_chart(msi)
        st.markdown('#### 📋 完整结果')
        st.dataframe(df, use_container_width=True, hide_index=True, height=500)
        # 下载
        st.markdown('#### 📥 下载报告')
        ts = datetime.now().strftime('%Y%m%d_%H%M%S')
        c1, c2 = st.columns(2)
        with c1:
            st.download_button('📥 CSV (完整数据)', df.to_csv(index=False, encoding='utf-8-sig'),
                                file_name=f'{herb_label}_鉴定_{ts}.csv', mime='text/csv',
                                use_container_width=True)
        with c2:
            md = generate_report_md(df, herb_label, instrument)
            st.download_button('📥 Markdown (可读报告)', md,
                                file_name=f'{herb_label}_报告_{ts}.md', mime='text/markdown',
                                use_container_width=True)

with tab_about:
    st.markdown('''
### 📖 平台说明

**v6.0 通用版** - 不限定药材, 自动从 TCM-SM-MS DB 抽取所有药材作为候选。

**核心特性**:
- 📚 **数据库**: TCM-SM-MS DB (199,411 行, 48,886 化合物, 579 药材)
- 🎯 **鉴定等级**: MSI Level 1-4 (对齐文章标准)
- 🔬 **仪器自适应**: Orbitrap / Q-TOF / QQQ / Low-res (10/15/30/50 ppm)
- ⚛️ **同位素分布**: Cl/Br/S/Si 诊断元素真实同位素匹配
- 🧮 **算法**: Intensity-weighted cosine similarity + 化合物去重 (同名+同分子式)
- 📊 **可视化**: 评级 / MSI / 化合物类型 / Top 化合物

**使用流程**:
1. 选择药材 (单/多/全库)
2. 选择仪器 (决定 PPM 容差)
3. 上传正/负离子 xlsx
4. 点击「开始鉴定」
5. 在「结果分析」页下载 CSV / Markdown

**支持的 xlsx 列结构** (大小写不敏感):
- 必需: `Precursor M/z`
- 推荐: `出峰时间t/min`, `Peak_1_m/z` ~ `Peak_N_m/z`, `Peak_N_Intensity`

**部署**:
```bash
streamlit run streamlit_app.py
# 默认账号 ZY / 513513
```
''')
