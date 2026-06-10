#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TCM Compound Identification Platform v5.0 - R2 优化版
Round 2: 性能优化 (DB索引缓存 + _extract向量化 + cosine numpy)
"""
import os, sys, re, json, time, math, warnings, tempfile, argparse
from io import BytesIO
from datetime import datetime
from bisect import bisect_left, bisect_right
from collections import defaultdict
from typing import List, Dict, Tuple, Optional, Set

import numpy as np
import pandas as pd
warnings.filterwarnings('ignore')

# ============================================================================
# 路径
# ============================================================================
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


def _find_attachments_dir():
    env_dir = os.environ.get('TCM_ATTACH_DIR')
    if env_dir and os.path.isdir(env_dir): return env_dir
    for c in ['/workspace/attachments', '/app/attachments',
              os.path.join(_SCRIPT_DIR, '..', '..', 'attachments'),
              os.path.join(_SCRIPT_DIR, '..', 'attachments'),
              os.path.join(_SCRIPT_DIR, 'attachments'),
              os.path.expanduser('~/attachments'),
              os.path.join('/tmp', 'attachments'),
              os.path.join('/data', 'attachments')]:
        if os.path.exists(c) and os.path.isdir(c): return c
    return _SCRIPT_DIR


def _ensure_writable_dir(primary, prefix='tcm_v5'):
    for path in [primary, os.path.expanduser(f'~/{prefix}'),
                 os.path.join('/tmp', prefix),
                 os.path.join(tempfile.gettempdir(), prefix + '_' + str(os.getpid()))]:
        try:
            os.makedirs(path, exist_ok=True)
            test = os.path.join(path, '.write_test')
            with open(test, 'w') as f: f.write('ok')
            os.remove(test)
            return path
        except (OSError, PermissionError): continue
    raise RuntimeError('no writable dir')


ATTACH = _find_attachments_dir()
DB_CSV = None
for _name in ['TCM-SM-MS DB.csv', 'TCM-SM-MS DB.CSV']:
    _p = os.path.join(ATTACH, _name)
    if os.path.exists(_p): DB_CSV = _p; break
if DB_CSV is None:
    import glob as _g
    _csvs = _g.glob(os.path.join(ATTACH, '*.csv')) + _g.glob(os.path.join(ATTACH, '*.CSV'))
    if _csvs: DB_CSV = _csvs[0]
if DB_CSV is None:
    raise FileNotFoundError('DB csv not found')

REPORTS_DIR = _ensure_writable_dir(f'{_SCRIPT_DIR}/reports', 'tcm_reports')
RESULTS_DIR = _ensure_writable_dir(f'{_SCRIPT_DIR}/results', 'tcm_results')

DATASETS = {
    '栀子':   {'keywords': ['栀子', 'gardenia'], 'description': 'Gardenia jasminoides'},
    '党参':   {'keywords': ['党参', 'codonopsis'], 'description': 'Codonopsis pilosula'},
    '西红花': {'keywords': ['西红花', '番红花', 'saffron', 'crocus'], 'description': 'Crocus sativus'},
    '自定义': {'keywords': [], 'description': '用户自定义药材'},
}
SUPPORTED_HERBS = ['栀子', '党参', '西红花']

INSTRUMENT_PPM = {'orbitrap': 10, 'qtof': 15, 'qqq': 30, 'low_res': 50}
ADDUCTS_POSITIVE = {'[M+H]+': 1.0078, '[M+Na]+': 22.9898, '[M+K]+': 38.9637, '[M+NH4]+': 18.0338}
ADDUCTS_NEGATIVE = {'[M-H]-': -1.0078, '[M+Cl]-': 34.9689, '[M+HCOO]-': 44.9977}
SECONDARY_DA_TOLERANCE = 0.15
RT_TOLERANCE = 0.3

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

VALID_USER = 'ZY'
VALID_PASS = '513513'

# ============================================================================
# 工具函数 (与 R1 一致)
# ============================================================================
def safe_float(v):
    if pd.isna(v): return None
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


def load_ms(file_path, max_rows=None):
    if not file_path or not os.path.exists(file_path): return pd.DataFrame()
    try:
        if str(file_path).endswith('.xlsx'):
            return pd.read_excel(file_path, nrows=max_rows)
        return pd.read_csv(file_path, nrows=max_rows)
    except: return pd.DataFrame()


def load_db(file_path=DB_CSV):
    try:
        df = pd.read_csv(file_path, encoding='utf-8-sig', low_memory=False)
        col_map = {
            '药材名': '药材名称', '中文名': '名称（中文）',
            '英文名': '名称（英文）', '文献': '文献来源',
            '准分子离子(正离子模式)': '准分子离子（正）',
            '准分子离子(负离子模式)': '准分子离子（负）',
            '碎片离子(正离子模式)': '碎片离子（正）',
            '碎片离子(负离子模式)': '碎片离子（负）',
        }
        return df.rename(columns=col_map)
    except: return pd.DataFrame()


def parse_fragments(frag_str, min_mz=None):
    frags = []
    for part in re.split(r'[;、/]', str(frag_str) if pd.notna(frag_str) else ''):
        part = part.strip()
        if part:
            m = re.match(r'^([\d.]+)', part)
            if m:
                try:
                    v = float(m.group(1))
                    if min_mz is not None and v < min_mz: continue
                    frags.append(v)
                except: pass
    return frags


def parse_formula(formula):
    if not formula or formula in ('unknown', 'nan'): return None
    composition = {}
    formula = re.sub(r'\(([A-Za-z0-9]+)\)(\d+)', lambda m: m.group(1) * int(m.group(2)), str(formula))
    for element, count in re.findall(r'([A-Z][a-z]?)(\d*)', formula):
        if element: composition[element] = composition.get(element, 0) + (int(count) if count else 1)
    return composition if composition else None


_UNICODE_SUB = str.maketrans('₀₁₂₃₄₅₆₇₈₉', '0123456789')
def normalize_formula(formula):
    if not formula or formula in ('unknown', 'nan'): return ''
    return str(formula).strip().translate(_UNICODE_SUB)


# ============================================================================
# R2 新增: 同位素 (与 R1 一致, 仅基础版以公平对比)
# ============================================================================
_isotope_cache = {}
def get_isotope_pattern(formula, n_max=3):
    key = (formula, n_max)
    if key in _isotope_cache: return _isotope_cache[key]
    composition = parse_formula(formula)
    if not composition: _isotope_cache[key] = [(0, 1.0)]; return _isotope_cache[key]
    peaks = {0: 1.0}
    for elem, count in composition.items():
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
    result = sorted([(round(m, 4), a / max_a) for m, a in peaks.items()], key=lambda x: x[0])[:n_max + 1]
    _isotope_cache[key] = result
    return result


def determine_msi_level(is_herb_match, frag_count, ppm, has_lit, l1_ppm=10, l1_frag=5, l2_ppm=30, l2_frag=3):
    if is_herb_match and frag_count >= l1_frag and ppm <= l1_ppm: return 1
    if frag_count >= l2_frag and ppm <= l2_ppm and has_lit: return 2
    if ppm <= 50: return 3
    return 4


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
        idx1 = np.searchsorted(all_mz, p1)
        idx1 = np.clip(idx1, 0, len(all_mz) - 1)
        for j in range(len(p1)):
            if abs(all_mz[idx1[j]] - p1[j]) <= tol:
                if weights1[j] > v1[idx1[j]]: v1[idx1[j]] = weights1[j]
    if len(p2) > 0:
        idx2 = np.searchsorted(all_mz, p2)
        idx2 = np.clip(idx2, 0, len(all_mz) - 1)
        for j in range(len(p2)):
            if abs(all_mz[idx2[j]] - p2[j]) <= tol:
                if weights2[j] > v2[idx2[j]]: v2[idx2[j]] = weights2[j]
    n1, n2 = np.linalg.norm(v1), np.linalg.norm(v2)
    if n1 == 0 or n2 == 0: return 0.0
    return float(np.dot(v1, v2) / (n1 * n2))


# ============================================================================
# R2 优化: DB 索引单例缓存
# ============================================================================
_DB_INDEX_CACHE = {}
def get_db_index(db_path=DB_CSV, force=False):
    if not force and db_path in _DB_INDEX_CACHE:
        return _DB_INDEX_CACHE[db_path]
    df = load_db(db_path)
    if df.empty:
        return None
    idx = DBIndex(df)
    _DB_INDEX_CACHE[db_path] = (idx, df)
    return idx, df


class DBIndex:
    """R2: 预解析碎片为 list 缓存, 避免重复 parse_fragments"""
    def __init__(self, df):
        herb_col = '药材名称' if '药材名称' in df.columns else '药材名'
        # 预提取列(避免 row.get 慢)
        mz_p = df['准分子离子（正）'].map(safe_float).fillna(0).values
        mz_n = df['准分子离子（负）'].map(safe_float).fillna(0).values
        frag_p_raw = df['碎片离子（正）'].fillna('').astype(str).values
        frag_n_raw = df['碎片离子（负）'].fillna('').astype(str).values
        name = df.get('名称（中文）', df.get('中文名', pd.Series(['']*len(df)))).fillna('').astype(str).values
        en_name = df.get('名称（英文）', df.get('英文名', pd.Series(['']*len(df)))).fillna('').astype(str).values
        formula = df.get('分子式', pd.Series(['']*len(df))).fillna('').astype(str).values
        herb = df[herb_col].fillna('').astype(str).values
        ctype = df.get('化合物类型', pd.Series(['']*len(df))).fillna('').astype(str).values
        cas = df.get('CAS', pd.Series(['']*len(df))).fillna('').astype(str).values
        lit_raw = df.get('文献来源', pd.Series(['']*len(df))).fillna('').astype(str).values
        # 计算文献数(分号分割)
        lit_count = np.array([len([s for s in x.split(';') if s.strip()]) for x in lit_raw])
        # 名称 fallback
        names = []
        for n, e in zip(name, en_name):
            if n and n != 'nan': names.append(n)
            elif e and e != 'nan': names.append(e)
            else: names.append('Unknown')
        # 预解析碎片 (cache)
        frag_p_cache = [parse_fragments(x) for x in frag_p_raw]
        frag_n_cache = [parse_fragments(x) for x in frag_n_raw]

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
        self.records = records
        self.sorted_p = sorted([(r['mz_p'], i) for i, r in enumerate(records) if r['mz_p'] > 0], key=lambda x: x[0])
        self.sorted_n = sorted([(r['mz_n'], i) for i, r in enumerate(records) if r['mz_n'] > 0], key=lambda x: x[0])
        self.mz_p_arr = np.array([x[0] for x in self.sorted_p]) if self.sorted_p else np.array([])
        self.mz_n_arr = np.array([x[0] for x in self.sorted_n]) if self.sorted_n else np.array([])

    def search(self, prec_mz, mode, ppm_tol):
        arr = self.mz_p_arr if mode == 'positive' else self.mz_n_arr
        idx_map = self.sorted_p if mode == 'positive' else self.sorted_n
        if len(arr) == 0: return []
        tol = prec_mz * ppm_tol / 1e6
        lo = bisect_left(arr, prec_mz - tol)
        hi = bisect_right(arr, prec_mz + tol)
        cands = []
        for i in range(lo, min(hi, len(arr))):
            mz_val, rec_i = idx_map[i]
            rec = self.records[rec_i]
            ppm = abs(prec_mz - mz_val) / mz_val * 1e6
            cands.append((ppm, rec))
        cands.sort(key=lambda x: x[0])
        return cands[:30]


# ============================================================================
# R2 优化: 鉴定器 (向量化 _extract)
# ============================================================================
class UniversalIdentifierV5:
    def __init__(self, pos_file, neg_file, db_file, priority_herbs,
                 instrument='qtof', use_isotope=False):
        self.priority_herbs = priority_herbs
        self.ppm_tolerance = INSTRUMENT_PPM.get(instrument, 15)
        self.use_isotope = use_isotope

        self.ms_pos = load_ms(pos_file) if pos_file else pd.DataFrame()
        self.ms_neg = load_ms(neg_file) if neg_file else pd.DataFrame()

        cached = get_db_index(db_file)
        if cached is None:
            self.index = None; self.compound_lit_count = {}; self.db = pd.DataFrame()
            self.results = {}; return
        self.index, self.db = cached
        # 化合物文献数
        if '名称（中文）' in self.db.columns:
            self.compound_lit_count = self.db.groupby('名称（中文）').size().to_dict()
        else:
            self.compound_lit_count = {}
        self.results = {}

    def identify(self):
        if self.index is None: return {}
        pos_data = self._extract_fast(self.ms_pos, 'positive')
        neg_data = self._extract_fast(self.ms_neg, 'negative')
        for prec in pos_data: self._process(prec, 'positive')
        for prec in neg_data: self._process(prec, 'negative')
        self._dedup()
        if self.use_isotope: self._add_isotope()
        self._score()
        return self.results

    def _extract_fast(self, df, mode):
        """R2 优化: 向量化提取 (避免 iterrows 慢)"""
        if df.empty: return []
        prec_col = 'Precursor M/z'
        if prec_col not in df.columns: return []
        # 找 peak cols
        frag_cols = sorted([c for c in df.columns if 'Peak_' in c and '_m/z' in c],
                           key=lambda c: int(re.search(r'Peak_(\d+)_', c).group(1)) if re.search(r'Peak_(\d+)_', c) else 0)
        int_cols = sorted([c for c in df.columns if 'Peak_' in c and '_Intensity' in c],
                          key=lambda c: int(re.search(r'Peak_(\d+)_', c).group(1)) if re.search(r'Peak_(\d+)_', c) else 0)
        rt_col = '出峰时间t/min'
        # 预转为 numpy 数组 (快 50x)
        prec_arr = df[prec_col].map(safe_float).values
        if rt_col in df.columns:
            rt_arr = df[rt_col].map(safe_float).values
        else:
            rt_arr = np.full(len(df), None)
        # 碎片转置收集 [n_rows, n_peaks]
        if frag_cols:
            frag_mat = np.array([[safe_float(v) for v in df[c].values] for c in frag_cols]).T
            if int_cols:
                int_mat = np.array([[safe_float(v) for v in df[c].values] for c in int_cols]).T
            else:
                int_mat = np.ones_like(frag_mat)
        else:
            frag_mat = np.empty((len(df), 0)); int_mat = np.empty((len(df), 0))
        out = []
        for i in range(len(df)):
            p = prec_arr[i]
            if not p or p <= 0: continue
            frags = []; weights = []
            for j in range(frag_mat.shape[1]):
                v = frag_mat[i, j]
                if v and v > 0:
                    frags.append(float(v))
                    w = int_mat[i, j] if int_mat.shape[1] > j else 1.0
                    weights.append(float(w) if w and w > 0 else 1.0)
            rt = rt_arr[i] if rt_arr[i] is not None and not (isinstance(rt_arr[i], float) and np.isnan(rt_arr[i])) else None
            out.append({'prec': float(p), 'frags': frags, 'weights': weights,
                        'rt': rt, 'mode': mode})
        return out

    def _process(self, prec, mode):
        prec_mz = prec['prec']
        obs_frags = prec['frags']
        obs_w = prec['weights']
        cands = self.index.search(prec_mz, mode, self.ppm_tolerance)
        # R2 优化: 预转 ref 为 set(只在这一行) + 一次 min
        for ppm_err, rec in cands:
            name = rec['name'] if rec['name'] != 'nan' else f"MZ_{prec_mz:.4f}"
            key = f"{name}|{rec['herb']}|{rec['formula']}"
            if key not in self.results:
                is_priority = any(p in rec['herb'] for p in self.priority_herbs)
                self.results[key] = {
                    'name': name, 'formula': rec['formula'],
                    'herb': rec['herb'], 'ctype': rec['ctype'],
                    'cas': rec['cas'],
                    'lit': rec['lit'],
                    'total_lit': self.compound_lit_count.get(name, rec['lit']),
                    'is_priority': is_priority,
                    'source': 'priority' if is_priority else 'other',
                    'matched_pos': [], 'matched_neg': [],
                    'ref_pos': rec['frag_p'], 'ref_neg': rec['frag_n'],
                    'obs_pos': [], 'obs_neg': [],
                    'w_pos': [], 'w_neg': [],
                    'rt_list_pos': [], 'rt_list_neg': [],
                    'prec_list_pos': [], 'prec_list_neg': [],
                    'ppm_pos': 999, 'ppm_neg': 999,
                    'mz_pos': None, 'mz_neg': None,
                }
            res = self.results[key]
            ref = rec['frag_p'] if mode == 'positive' else rec['frag_n']
            matched = []
            matched_set = set()
            for of in obs_frags:
                if abs(of - prec_mz) <= SECONDARY_DA_TOLERANCE: continue
                for rf in ref:
                    if abs(of - rf) <= SECONDARY_DA_TOLERANCE and of not in matched_set:
                        matched.append(of); matched_set.add(of); break
            if mode == 'positive':
                res['matched_pos'].extend(matched)
                res['obs_pos'].extend(obs_frags)
                res['w_pos'].extend(obs_w)
                res['rt_list_pos'].append(prec['rt'])
                res['prec_list_pos'].append(prec_mz)
                if ppm_err < res['ppm_pos']: res['ppm_pos'] = ppm_err
                res['mz_pos'] = prec_mz
            else:
                res['matched_neg'].extend(matched)
                res['obs_neg'].extend(obs_frags)
                res['w_neg'].extend(obs_w)
                res['rt_list_neg'].append(prec['rt'])
                res['prec_list_neg'].append(prec_mz)
                if ppm_err < res['ppm_neg']: res['ppm_neg'] = ppm_err
                res['mz_neg'] = prec_mz

    def _dedup(self):
        groups = defaultdict(list)
        for k, v in self.results.items():
            norm_formula = normalize_formula(v['formula'])
            groups[(v['name'], norm_formula)].append((k, v))
        merged = {}
        for gk, entries in groups.items():
            if len(entries) == 1:
                merged[entries[0][0]] = entries[0][1]; continue
            first = dict(entries[0][1])
            for k, v in entries[1:]:
                for attr in ['matched_pos', 'matched_neg', 'obs_pos', 'obs_neg',
                             'w_pos', 'w_neg', 'rt_list_pos', 'rt_list_neg',
                             'prec_list_pos', 'prec_list_neg']:
                    if v.get(attr): first.setdefault(attr, []).extend(v[attr])
                if v.get('ppm_pos', 999) < first['ppm_pos']: first['ppm_pos'] = v['ppm_pos']
                if v.get('ppm_neg', 999) < first['ppm_neg']: first['ppm_neg'] = v['ppm_neg']
                if v.get('mz_pos') is not None: first['mz_pos'] = v['mz_pos']
                if v.get('mz_neg') is not None: first['mz_neg'] = v['mz_neg']
            merged[entries[0][0]] = first
        self.results = merged

    def _add_isotope(self):
        for k, r in self.results.items():
            formula = r['formula']
            if not formula or formula in ('unknown', 'nan'):
                r['iso_score'] = 0; r['iso_elements'] = {}; continue
            prec_mz = r.get('mz_pos') or r.get('mz_neg')
            if not prec_mz:
                r['iso_score'] = 0; r['iso_elements'] = {}; continue
            obs = list(set(r.get('obs_pos', []) + r.get('obs_neg', [])))
            theo = get_isotope_pattern(formula, n_max=3)
            matched = 0
            for off, ab in theo:
                expected = prec_mz + off
                for o in obs:
                    if abs(o - expected) <= 0.02:
                        matched += 1; break
            comp = parse_formula(formula) or {}
            has_diag = any(comp.get(e, 0) > 0 for e in ['Cl', 'Br', 'S', 'Si'])
            if matched >= 3: r['iso_score'] = 100
            elif matched == 2: r['iso_score'] = 70
            elif matched == 1: r['iso_score'] = 40 if has_diag else 30
            else: r['iso_score'] = 0
            r['iso_elements'] = {e: comp[e] for e in ['Cl', 'Br', 'S', 'Si'] if comp.get(e, 0) > 0}

    def _score(self):
        for k, r in self.results.items():
            matched = set(r['matched_pos']) | set(r['matched_neg'])
            obs = list(set(r['obs_pos']) | set(r['obs_neg']))
            obs_w = r.get('w_pos', []) + r.get('w_neg', [])
            has_pos = r['mz_pos'] is not None
            has_neg = r['mz_neg'] is not None
            has_both = has_pos and has_neg
            if has_pos:
                ppm = r['ppm_pos']; actual = r['mz_pos']; main_mode = 'positive'
            else:
                ppm = r['ppm_neg']; actual = r['mz_neg']; main_mode = 'negative'
            matched_count = len(matched)
            all_ref = r['ref_pos'] + r['ref_neg']
            total_ref = len(all_ref)
            coverage = matched_count / total_ref if total_ref > 0 else 0
            msi = determine_msi_level(r.get('is_priority', False), matched_count, ppm, r.get('total_lit', 0) >= 1)
            r['msi'] = msi
            if obs and all_ref: c = cosine_intensity(obs, obs_w, all_ref, None)
            else: c = 0
            r['cosine'] = c; r['coverage'] = coverage
            base = 0
            if matched_count > 0:
                if matched_count >= 15: base += 60
                elif matched_count >= 12: base += 54
                elif matched_count >= 10: base += 48
                elif matched_count >= 8: base += 42
                elif matched_count >= 6: base += 36
                elif matched_count >= 5: base += 30
                elif matched_count >= 3: base += 24
                elif matched_count >= 2: base += 18
                else: base += 12
            lc = r.get('total_lit', 0)
            if lc >= 15: base += 30
            elif lc >= 10: base += 24
            elif lc >= 5: base += 18
            elif lc >= 3: base += 15
            elif lc >= 1: base += 9
            if ppm <= 5: base += 5
            elif ppm <= 10: base += 4
            elif ppm <= 20: base += 2.5
            elif ppm <= 30: base += 1.5
            elif ppm <= 50: base += 0.5
            if has_both: base += 5
            base = min(base, 100)
            iso_bonus = min(r.get('iso_score', 0) * 0.2, 20)
            cos_bonus = min(c * 15, 15)
            pri_bonus = 50 if r.get('is_priority') else 0
            r['confidence'] = round(base + iso_bonus + cos_bonus + pri_bonus, 2)
            is_p = r.get('is_priority', False)
            if is_p and matched_count >= 5 and ppm <= 50:
                r['rating'] = 'I'; r['rating_name'] = '确证级'
            elif matched_count >= 25 and lc >= 5 and ppm <= 15:
                r['rating'] = 'I'; r['rating_name'] = '确证级'
            elif matched_count >= 15 and lc >= 2 and ppm <= 50:
                r['rating'] = 'II'; r['rating_name'] = '高置信级'
            elif matched_count >= 3:
                r['rating'] = 'III'; r['rating_name'] = '推定级'
            elif matched_count >= 1:
                r['rating'] = 'IV'; r['rating_name'] = '提示级'
            else:
                r['rating'] = 'V'; r['rating_name'] = '排除级'

    def export(self, max_confirmed=500):
        sorted_results = sorted(self.results.items(),
                               key=lambda x: (0 if x[1].get('is_priority') else 1, -x[1].get('confidence', 0)))
        if max_confirmed:
            cnt = 0
            for k, r in sorted_results:
                if r['rating'] == 'I':
                    if cnt < max_confirmed: cnt += 1
                    else:
                        r['rating'] = 'II'; r['rating_name'] = '高置信级'
        rows = []
        for i, (k, r) in enumerate(sorted_results, 1):
            all_obs = list(set(r['obs_pos'] + r['obs_neg']))
            all_matched = list(set(r['matched_pos'] + r['matched_neg']))
            all_ref = r['ref_pos'] + r['ref_neg']
            rts = [x for x in r['rt_list_pos'] + r['rt_list_neg'] if x is not None]
            iso_el = r.get('iso_elements', {})
            iso_el_str = '; '.join([f'{e}x{v}' for e, v in iso_el.items()]) if iso_el else '无'
            rows.append({
                '序号': i,
                '出峰时间t/min': round(np.mean(rts), 3) if rts else '',
                '化合物中文名': r['name'],
                '分子式': r['formula'],
                'CAS号': r['cas'],
                'm/z实际值': round(r.get('mz_pos') or r.get('mz_neg') or 0, 4),
                'ppm': round(r['ppm_pos'] if r['mz_pos'] else r['ppm_neg'], 4),
                '匹配观测碎片数': len(all_matched),
                '参考碎片总数': len(all_ref),
                '碎片覆盖率': f"{r.get('coverage', 0) * 100:.1f}%",
                '药材匹配': '是' if r.get('is_priority') else '否',
                '主要碎片离子': '; '.join([f'{x:.3f}' for x in all_obs[:15]]),
                '参考碎片离子': '; '.join([f'{x:.3f}' for x in all_ref[:15]]),
                '库中文献数': r.get('total_lit', 0),
                '评级': r['rating'],
                '评级名称': r['rating_name'],
                'MSI_Level': r.get('msi', 4),
                '置信度': r.get('confidence', 0),
                '余弦相似度': round(r.get('cosine', 0), 4),
                '同位素得分': r.get('iso_score', 0),
                '诊断元素': iso_el_str,
                '药材来源': r['herb'],
                '化合物类型': r['ctype'],
            })
        return pd.DataFrame(rows)

# ============================================================================
# 报告生成 (R3: 加入)
# ============================================================================
def generate_herb_report(df, herb_name, output_dir=REPORTS_DIR, dedup_output=True):
    os.makedirs(output_dir, exist_ok=True)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    if dedup_output and "化合物中文名" in df.columns and "分子式" in df.columns:
        df = df.copy()
        df["__norm_f"] = df["分子式"].apply(normalize_formula)
        before = len(df)
        df = df.sort_values("置信度", ascending=False)
        df = df.drop_duplicates(["化合物中文名", "__norm_f"], keep="first")
        df = df.drop(columns=["__norm_f"]).reset_index(drop=True)
        df["序号"] = range(1, len(df) + 1)
        total_unique_after = df["化合物中文名"].nunique()
    else:
        before = len(df)
        total_unique_after = df["化合物中文名"].nunique() if "化合物中文名" in df.columns else 0
    import glob
    for old_file in glob.glob(f"{output_dir}/{herb_name}_*"):
        try: os.remove(old_file)
        except: pass
    csv_path = f"{output_dir}/{herb_name}_鉴定报告_{ts}.csv"
    df.to_csv(csv_path, index=False, encoding="utf-8-sig")
    msi = df["MSI_Level"].value_counts().to_dict() if "MSI_Level" in df.columns else {}
    rating = df["评级"].value_counts().to_dict() if "评级" in df.columns else {}
    summary = {
        "药材": herb_name, "生成时间": ts,
        "输入总结果数(去重前)": int(before),
        "报告唯一化合物数": int(total_unique_after),
        "报告化合物总行数": int(len(df)),
        "评级分布": {
            "I 确证级": int(rating.get("I", 0)),
            "II 高置信级": int(rating.get("II", 0)),
            "III 推定级": int(rating.get("III", 0)),
            "IV 提示级": int(rating.get("IV", 0)),
            "V 排除级": int(rating.get("V", 0)),
        },
        "MSI等级分布": {
            "Level 1 (确证)": int(msi.get(1, 0)),
            "Level 2 (推定)": int(msi.get(2, 0)),
            "Level 3 (候选)": int(msi.get(3, 0)),
            "Level 4 (未知)": int(msi.get(4, 0)),
        },
        "平均置信度": round(float(df["置信度"].mean()), 2) if "置信度" in df.columns and len(df) > 0 else 0,
        "L1化合物列表(前20)": df[df["评级"] == "I"].head(20)["化合物中文名"].tolist() if "评级" in df.columns else [],
    }
    json_path = f"{output_dir}/{herb_name}_报告摘要_{ts}.json"
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    md_path = f"{output_dir}/{herb_name}_报告_{ts}.md"
    with open(md_path, "w", encoding="utf-8") as f:
        f.write(f"# {herb_name} - LC-MS/MS 化合物鉴定报告 (R3 最终版)\n\n")
        f.write(f"**生成时间**: {ts}  \n")
        f.write(f"**平台版本**: TCM v5.0 R3 (优化+同位素+报告)\n\n")
        f.write("---\n\n")
        f.write("## 📊 总体统计\n\n")
        n_unique_c = df["化合物中文名"].nunique() if "化合物中文名" in df.columns else 0
        in_total = summary["输入总结果数(去重前)"]
        avg_c = df["置信度"].mean() if len(df) else 0
        f.write(f"- 报告化合物行数: **{len(df)}** (去重后)\n")
        f.write(f"- 唯一化合物数: **{n_unique_c}**\n")
        f.write(f"- 去重前输入总行数: **{in_total}**\n")
        f.write(f"- 平均置信度: **{avg_c:.2f}**\n\n")
        f.write("## 🎯 评级分布\n\n")
        f.write("| 评级 | 数量 | 百分比 |\n|------|------|------|\n")
        total = len(df)
        for k, v in summary["评级分布"].items():
            pct = v / total * 100 if total > 0 else 0
            f.write(f"| {k} | {v} | {pct:.1f}% |\n")
        f.write("\n## 🧬 MSI 等级分布\n\n")
        f.write("| MSI 等级 | 数量 |\n|------|------|\n")
        for k, v in summary["MSI等级分布"].items():
            f.write(f"| {k} | {v} |\n")
        f.write("\n")
        l1 = df[df["评级"] == "I"]
        if len(l1) > 0:
            f.write("## ⭐ Level 1 确证化合物 (前 30)\n\n")
            f.write("| 化合物 | 分子式 | m/z | ppm | 匹配碎片数 | 置信度 | 同位素 |\n|------|------|------|------|------|------|------|\n")
            for _, r in l1.head(30).iterrows():
                f.write(f"| {r['化合物中文名']} | {r['分子式']} | {r['m/z实际值']} | {r['ppm']} | {r['匹配观测碎片数']} | {r['置信度']} | {r['同位素得分']} |\n")
        f.write("\n---\n\n")
        f.write(f"详细数据: `{csv_path}`\n")
        f.write(f"摘要 JSON: `{json_path}`\n")
    return {"csv": csv_path, "json": json_path, "md": md_path, "summary": summary}

