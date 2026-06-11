#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
🌿 中药 LC-MS/MS 化合物通用鉴定平台 v6.0 - Streamlit
====================================================

通用版: 不限定药材, 自动从 TCM-SM-MS DB 抽取所有药材名作为候选

启动: streamlit run streamlit_app.py
默认账号: ZY / 513513
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
def determine_msi_level(is_herb_match, frag_count, ppm, has_lit,
                        l1_ppm=10, l1_frag=5, l2_ppm=30, l2_frag=3):
    if is_herb_match and frag_count >= l1_frag and ppm <= l1_ppm: return 1
    if frag_count >= l2_frag and ppm <= l2_ppm and has_lit: return 2
    if ppm <= 50: return 3
    return 4


def identify(extracted_pos, extracted_neg, priority_herbs, instrument='qtof', use_isotope=True):
    """通用鉴定: 1 个 pos + 1 个 neg (可空), 返回 DataFrame"""
    db = build_db_pkl_cached()
    records = db['records']
    sorted_p = db['sorted_p']
    sorted_n = db['sorted_n']
    mz_p_arr = db['mz_p_arr']
    mz_n_arr = db['mz_n_arr']
    compound_lit_count = db['compound_lit_count']
    ppm_tol = INSTRUMENT_PPM.get(instrument, 15)

    priority_set = set(priority_herbs) if priority_herbs else set()
    results = {}

    def process_one(extracted, mode):
        arr = mz_p_arr if mode == 'positive' else mz_n_arr
        idx_map = sorted_p if mode == 'positive' else sorted_n
        if len(arr) == 0: return
        for prec_data in extracted:
            prec_mz = prec_data['prec']
            obs_frags = prec_data['frags']
            obs_w = prec_data['weights']
            tol = prec_mz * ppm_tol / 1e6
            lo = bisect_left(arr, prec_mz - tol)
            hi = bisect_right(arr, prec_mz + tol)
            for i in range(lo, min(hi, len(arr))):
                mz_val, rec_i = idx_map[i]
                rec = records[rec_i]
                ppm_err = abs(prec_mz - mz_val) / mz_val * 1e6
                name = rec['name']
                if name == 'nan': name = f"MZ_{prec_mz:.4f}"
                # 关键: priority_herbs 决定是否是 priority
                is_priority = bool(priority_set) and any(p in rec['herb'] for p in priority_set)
                key = f"{name}|{rec['herb']}|{rec['formula']}"
                if key not in results:
                    results[key] = {
                        'name': name, 'formula': rec['formula'],
                        'herb': rec['herb'], 'ctype': rec['ctype'],
                        'cas': rec['cas'], 'lit': rec['lit'],
                        'total_lit': compound_lit_count.get(name, rec['lit']),
                        'is_priority': is_priority,
                        'ref_pos': rec['frag_p'], 'ref_neg': rec['frag_n'],
                        'matched_pos': set(), 'matched_neg': set(),
                        'obs_pos': [], 'obs_neg': [],
                        'w_pos': [], 'w_neg': [],
                        'rt_pos': [], 'rt_neg': [],
                        'ppm_pos': 999, 'ppm_neg': 999,
                        'mz_pos': None, 'mz_neg': None,
                    }
                res = results[key]
                ref = rec['frag_p'] if mode == 'positive' else rec['frag_n']
                matched = set()
                for of in obs_frags:
                    if abs(of - prec_mz) <= SECONDARY_DA_TOLERANCE: continue
                    for rf in ref:
                        if abs(of - rf) <= SECONDARY_DA_TOLERANCE and of not in matched:
                            matched.add(of); break
                if mode == 'positive':
                    res['matched_pos'].update(matched)
                    res['obs_pos'].extend(obs_frags)
                    res['w_pos'].extend(obs_w)
                    res['rt_pos'].append(prec_data['rt'])
                    if ppm_err < res['ppm_pos']: res['ppm_pos'] = ppm_err
                    res['mz_pos'] = prec_mz
                else:
                    res['matched_neg'].update(matched)
                    res['obs_neg'].extend(obs_frags)
                    res['w_neg'].extend(obs_w)
                    res['rt_neg'].append(prec_data['rt'])
                    if ppm_err < res['ppm_neg']: res['ppm_neg'] = ppm_err
                    res['mz_neg'] = prec_mz

    if extracted_pos: process_one(extracted_pos, 'positive')
    if extracted_neg: process_one(extracted_neg, 'negative')
    if not results: return pd.DataFrame()
    # dedup
    groups = defaultdict(list)
    for k, v in results.items():
        norm = normalize_formula(v['formula'])
        groups[(v['name'], norm)].append((k, v))
    merged = {}
    for gk, entries in groups.items():
        if len(entries) == 1:
            merged[entries[0][0]] = entries[0][1]; continue
        first = dict(entries[0][1])
        for k, v in entries[1:]:
            for attr in ['matched_pos', 'matched_neg', 'obs_pos', 'obs_neg',
                         'w_pos', 'w_neg', 'rt_pos', 'rt_neg']:
                if v.get(attr):
                    s = first.setdefault(attr, set() if attr in ['matched_pos', 'matched_neg'] else [])
                    if isinstance(s, set): s.update(v[attr])
                    else: s.extend(v[attr])
            if v.get('ppm_pos', 999) < first['ppm_pos']: first['ppm_pos'] = v['ppm_pos']
            if v.get('ppm_neg', 999) < first['ppm_neg']: first['ppm_neg'] = v['ppm_neg']
            if v.get('mz_pos') is not None: first['mz_pos'] = v['mz_pos']
            if v.get('mz_neg') is not None: first['mz_neg'] = v['mz_neg']
        merged[entries[0][0]] = first
    # 同位素 (可选)
    if use_isotope:
        _iso_cache = {}
        for k, r in merged.items():
            formula = r['formula']
            if not formula or formula in ('unknown', 'nan'):
                r['iso_score'] = 0; r['iso_elements'] = ''; continue
            prec_mz = r.get('mz_pos') or r.get('mz_neg')
            if not prec_mz:
                r['iso_score'] = 0; r['iso_elements'] = ''; continue
            key = (formula, 3)
            if key not in _iso_cache:
                _iso_cache[key] = _compute_isotope_pattern(formula, n_max=3)
            theo = _iso_cache[key]
            obs = list(set(r.get('obs_pos', []) + r.get('obs_neg', [])))
            matched = sum(1 for off, _ in theo if any(abs(o - (prec_mz + off)) <= 0.02 for o in obs))
            comp = parse_formula(formula) or {}
            has_diag = any(comp.get(e, 0) > 0 for e in ['Cl', 'Br', 'S', 'Si'])
            if matched >= 3: r['iso_score'] = 100
            elif matched == 2: r['iso_score'] = 70
            elif matched == 1: r['iso_score'] = 40 if has_diag else 30
            else: r['iso_score'] = 0
            r['iso_elements'] = '; '.join([f'{e}x{comp[e]}' for e in ['Cl', 'Br', 'S', 'Si'] if comp.get(e, 0) > 0]) or '无'
    else:
        for r in merged.values():
            r['iso_score'] = 0; r['iso_elements'] = '无'
    # 评分
    out = []
    for k, r in merged.items():
        matched = r['matched_pos'] | r['matched_neg']
        obs = list(set(r['obs_pos'] + r['obs_neg']))
        obs_w = r['w_pos'] + r['w_neg']
        all_ref = r['ref_pos'] + r['ref_neg']
        total_ref = len(all_ref)
        mc = len(matched)
        coverage = mc / total_ref if total_ref > 0 else 0
        if r['mz_pos']: ppm = r['ppm_pos']
        else: ppm = r['ppm_neg']
        msi = determine_msi_level(r['is_priority'], mc, ppm, r['total_lit'] >= 1)
        c = cosine_intensity(obs, obs_w, all_ref, None) if obs and all_ref else 0
        base = 0
        if mc > 0:
            if mc >= 15: base += 60
            elif mc >= 10: base += 48
            elif mc >= 5: base += 30
            elif mc >= 3: base += 24
            elif mc >= 2: base += 18
            else: base += 12
        lc = r['total_lit']
        if lc >= 15: base += 30
        elif lc >= 5: base += 18
        elif lc >= 1: base += 9
        if ppm <= 10: base += 4
        elif ppm <= 30: base += 1.5
        elif ppm <= 50: base += 0.5
        if r['mz_pos'] and r['mz_neg']: base += 5
        base = min(base, 100)
        iso_bonus = min(r.get('iso_score', 0) * 0.2, 20)
        pri = 50 if r['is_priority'] else 0
        conf = round(base + iso_bonus + min(c * 15, 15) + pri, 2)
        if r['is_priority'] and mc >= 5 and ppm <= 50:
            rating, rating_name = 'I', '确证级'
        elif mc >= 25 and lc >= 5 and ppm <= 15:
            rating, rating_name = 'I', '确证级'
        elif mc >= 15 and lc >= 2 and ppm <= 50:
            rating, rating_name = 'II', '高置信级'
        elif mc >= 3: rating, rating_name = 'III', '推定级'
        elif mc >= 1: rating, rating_name = 'IV', '提示级'
        else: rating, rating_name = 'V', '排除级'
        rts = [x for x in r['rt_pos'] + r['rt_neg'] if x is not None]
        out.append({
            '化合物中文名': r['name'],
            '分子式': r['formula'],
            'CAS号': r['cas'],
            'm/z实际值': round(r.get('mz_pos') or r.get('mz_neg') or 0, 4),
            'ppm': round(r['ppm_pos'] if r['mz_pos'] else r['ppm_neg'], 4),
            '匹配观测碎片数': mc,
            '参考碎片总数': total_ref,
            '碎片覆盖率': f"{coverage * 100:.1f}%",
            '出峰时间t/min': round(np.mean(rts), 3) if rts else '',
            '药材匹配': '是' if r['is_priority'] else '否',
            '主要碎片离子': '; '.join([f'{x:.3f}' for x in obs[:15]]),
            '参考碎片离子': '; '.join([f'{x:.3f}' for x in all_ref[:15]]),
            '库中文献数': lc,
            '评级': rating, '评级名称': rating_name,
            'MSI_Level': msi, '置信度': conf,
            '余弦相似度': round(c, 4),
            '同位素得分': r.get('iso_score', 0),
            '诊断元素': r.get('iso_elements', '无'),
            '药材来源': r['herb'], '化合物类型': r['ctype'],
        })
    df = pd.DataFrame(out)
    if len(df) > 0:
        df = df.sort_values(['药材匹配', '置信度'], ascending=[False, False]).reset_index(drop=True)
        df.insert(0, '序号', range(1, len(df) + 1))
    return df


def _compute_isotope_pattern(formula, n_max=3):
    comp = parse_formula(formula)
    if not comp: return [(0, 1.0)]
    peaks = {0: 1.0}
    for elem, count in comp.items():
        if elem not in ISOTOPE_ABUNDANCES: continue
        try: count = int(count)
        except: continue
        ab = ISOTOPE_ABUNDANCES[elem]
        mi = ISOTOPE_MASS_INCREMENT[elem]
        base = min(mi.keys(), key=lambda k: mi[k])
        for _ in range(count):
            new = {}
            for cm, ca in peaks.items():
                for im, ia in ab.items():
                    nm = cm if im == base else cm + mi.get(im, 0)
                    new[nm] = new.get(nm, 0) + ca * ia
            peaks = new
            if len(peaks) > 200:
                peaks = dict(sorted(peaks.items(), key=lambda x: -x[1])[:100])
    max_a = max(peaks.values()) if peaks else 1.0
    return sorted([(round(m, 4), a / max_a) for m, a in peaks.items()], key=lambda x: x[0])[:n_max + 1]


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
    page_title='中药 LC-MS/MS 通用鉴定平台 v6.0',
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
    st.markdown('<div class="main-title"><h1>🌿 中药 LC-MS/MS 化合物通用鉴定平台</h1>'
                '<p>基于 TCM-SM-MS DB (48,886 化合物 / 579 药材) 的 LC-MS/MS 智能鉴定</p></div>',
                unsafe_allow_html=True)
    col1, col2, col3 = st.columns([1, 2, 1]) if False else st.columns(3)
    with col2:
        with st.form('login'):
            st.markdown('### 🔐 登录')
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
        h = st.selectbox('选择药材', all_herbs, index=all_herbs.index('栀子') if '栀子' in all_herbs else 0)
        selected_herbs = [h]
    elif herb_mode == '多药材':
        h = st.multiselect('选择药材 (可多选)', all_herbs,
                            default=['栀子', '党参'] if '栀子' in all_herbs and '党参' in all_herbs else all_herbs[:2])
        selected_herbs = h
    else:
        selected_herbs = []
        st.caption('当前为全库模式, 所有化合物不享受 priority 加分')
    st.markdown('---')
    st.markdown('### 📈 平台信息')
    st.metric('TCM-SM-MS DB', f'{len(all_herbs)} 药材')
    st.metric('支持仪器', f'{len(INSTRUMENT_PPM)} 种')
    st.metric('当前 PPM', INSTRUMENT_PPM[instrument])
    if st.button('登出', use_container_width=True):
        st.session_state.logged_in = False
        st.rerun()

# 主区
st.markdown('<div class="main-title"><h1>🌿 中药 LC-MS/MS 化合物通用鉴定平台 v6.0</h1>'
            '<p>通用版 | 自动适配 DB 内 579 种药材 | MSI Level 1-4 鉴定</p></div>',
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
                              instrument=instrument, use_isotope=use_isotope)
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
