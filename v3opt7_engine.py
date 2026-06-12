#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
v3opt7: v3opt6 + RT 铁证升 L1
- 新升级: mc ≥ 1 + RT 一致性 ≥ 0.99 + obs 谱图数 ≥ 3 → 可从 L3 升 L1 (跳过 L2)
- 解决 v3opt6 "iso=0 的化合物永远升不了 L1" 的过严问题
- 保留 v3opt6 的三大升级: 同位素比例验证 + RT 补救 + 单点精降级
"""
import re
import numpy as np
from bisect import bisect_left, bisect_right
from collections import defaultdict

_UNICODE_SUB = str.maketrans('₀₁₂₃₄₅₆₇₈₉', '0123456789')
def normalize_formula(s):
    if not s: return ''
    return str(s).strip().translate(_UNICODE_SUB)


ADDUCTS = {
    'positive': [
        ('[M+H]+',   1.007276), ('[M+Na]+',  22.989218),
        ('[M+NH4]+', 17.026547), ('[M+K]+',   38.963158),
    ],
    'negative': [
        ('[M-H]-',  -1.007276), ('[M+Cl]-',  34.969402), ('[M+FA-H]-', 44.998),
    ],
}
AW = {'H': 1.007825, 'C': 12.000000, 'N': 14.003074, 'O': 15.994915,
      'S': 31.972071, 'P': 30.973762, 'Cl': 34.968853, 'Br': 78.918336,
      'F': 18.998403, 'I': 126.904468, 'Si': 27.976926, 'B': 10.811}
ISOTOPE_NATURAL = {
    'C': {13: 0.0107}, 'H': {2: 0.00015}, 'N': {15: 0.0037},
    'O': {17: 0.0004, 18: 0.0020}, 'S': {33: 0.0076, 34: 0.0429},
    'Cl': {37: 0.2423}, 'Br': {81: 0.4931}, 'Si': {29: 0.0467, 30: 0.0310},
}
ISOTOPE_MASS_INC = {
    'C': {13: 1.0034}, 'H': {2: 1.0063}, 'N': {15: 0.9970},
    'O': {17: 1.0042, 18: 1.9978}, 'S': {33: 0.9994, 34: 1.9958},
    'Cl': {37: 1.9970}, 'Br': {81: 1.9980}, 'Si': {29: 0.9986, 30: 1.9942},
}


def get_molecular_weight(formula):
    if not formula: return 0
    f = formula.translate(_UNICODE_SUB)
    total = 0.0
    for elem, count in re.findall(r'([A-Z][a-z]?)(\d*)', f):
        if elem in AW:
            n = int(count) if count else 1
            total += AW[elem] * n
    return total


def get_element_counts(formula):
    if not formula: return {}
    f = formula.translate(_UNICODE_SUB)
    counts = {}
    for elem, count in re.findall(r'([A-Z][a-z]?)(\d*)', f):
        if elem:
            n = int(count) if count else 1
            counts[elem] = counts.get(elem, 0) + n
    return counts


def compute_isotope_pattern_detailed(formula):
    counts = get_element_counts(formula)
    peaks = {0: 1.0}
    for elem, nat in ISOTOPE_NATURAL.items():
        n = counts.get(elem, 0)
        if n == 0: continue
        for iso, ab in nat.items():
            inc = ISOTOPE_MASS_INC.get(elem, {}).get(iso, 0)
            if n == 1:
                peaks[inc] = peaks.get(inc, 0) + ab
            else:
                peaks[inc] = peaks.get(inc, 0) + n * ab
                if len(nat) == 1 and n >= 2:
                    double_inc = 2 * inc
                    peaks[double_inc] = peaks.get(double_inc, 0) + (n * (n-1) / 2) * (ab ** 2)
    return peaks


def compute_isotope_pattern(formula):
    return compute_isotope_pattern_detailed(formula)


def check_isotope_match(formula, prec_mz, obs_frags, obs_w, tol=0.02):
    counts = get_element_counts(formula)
    has_diag = any(counts.get(e, 0) > 0 for e in ['Cl', 'Br', 'S', 'Si'])
    if not has_diag and counts.get('C', 0) < 3:
        return 0
    peaks = compute_isotope_pattern_detailed(formula)
    if not peaks: return 0
    m0_intensity = 0
    for j, of in enumerate(obs_frags):
        if abs(of - prec_mz) <= tol:
            m0_intensity = max(m0_intensity, obs_w[j] if j < len(obs_w) else 1.0)
    if m0_intensity == 0:
        m0_intensity = 1.0
    matched_count = 0
    ratio_score = 0.0
    n_ratios = 0
    for offset, expected_ab in peaks.items():
        if offset == 0: continue
        if expected_ab < 0.01: continue
        target_mz = prec_mz + offset
        found_intensity = 0
        for j, of in enumerate(obs_frags):
            if abs(of - target_mz) <= tol:
                found_intensity = max(found_intensity, obs_w[j] if j < len(obs_w) else 1.0)
                matched_count += 1
                break
        if found_intensity > 0:
            observed_ratio = found_intensity / m0_intensity
            expected_ratio = expected_ab
            if expected_ratio > 0:
                rel_err = abs(observed_ratio - expected_ratio) / expected_ratio
                if rel_err < 0.3:
                    ratio_score += 1.0
                elif rel_err < 0.6:
                    ratio_score += 0.5
                n_ratios += 1
    total_offsets = max(1, sum(1 for k in peaks if k > 0 and peaks[k] >= 0.01))
    presence_score = matched_count / total_offsets
    ratio_final = ratio_score / n_ratios if n_ratios > 0 else 0
    return presence_score * 0.6 + ratio_final * 0.4


def find_rt_partner(mz_target, extracted_partner, ppm_tol=20, rt_tol=0.5):
    if not extracted_partner: return (None, 0)
    tol = mz_target * ppm_tol / 1e6
    matches = []
    for prec, d in extracted_partner:
        if abs(prec - mz_target) <= tol and d['rt'] is not None:
            matches.append(d['rt'])
    if not matches: return (None, 0)
    rt_med = float(np.median(matches))
    return (rt_med, len(matches))


def fuse_records_v3opt7(records, name, formula):
    same = [r for r in records if r['name'] == name and r['formula'] == formula]
    if not same: return []
    ref_p_count = defaultdict(int)
    ref_n_count = defaultdict(int)
    for r in same:
        for x in r['frag_p']:
            ref_p_count[round(x, 4)] += 1
        for x in r['frag_n']:
            ref_n_count[round(x, 4)] += 1
    ref_p_fused = sorted(ref_p_count.keys(), key=lambda x: (-ref_p_count[x], x))
    ref_n_fused = sorted(ref_n_count.keys(), key=lambda x: (-ref_n_count[x], x))
    max_count_p = max(ref_p_count.values()) if ref_p_count else 1
    max_count_n = max(ref_n_count.values()) if ref_n_count else 1
    ref_p_weight = {x: ref_p_count[x] / max_count_p for x in ref_p_fused}
    ref_n_weight = {x: ref_n_count[x] / max_count_n for x in ref_n_fused}
    mw = get_molecular_weight(formula)
    herb_set = []
    for r in same:
        for h in re.split(r'[、;/,，\s]+', r['herb']):
            h = h.strip()
            if h and h not in herb_set: herb_set.append(h)
    herb_merged = '、'.join(herb_set)
    lit = max(r['lit'] for r in same)
    ctype = next((r['ctype'] for r in same if r['ctype'] and r['ctype'] != 'nan'), '')
    cas = next((r['cas'] for r in same if r['cas'] and r['cas'] != 'nan'), '')
    candidates = []
    for mode in ['positive', 'negative']:
        for adduct_name, adduct_delta in ADDUCTS[mode]:
            mz = mw + adduct_delta
            if mz <= 0: continue
            candidates.append({
                'name': name, 'formula': formula,
                'mz': mz, 'mode': mode, 'adduct': adduct_name,
                'ref_p': ref_p_fused if mode == 'positive' else [],
                'ref_n': ref_n_fused if mode == 'negative' else [],
                'ref_p_weight': ref_p_weight if mode == 'positive' else {},
                'ref_n_weight': ref_n_weight if mode == 'negative' else {},
                'herb': herb_merged, 'ctype': ctype, 'cas': cas,
                'lit': lit, 'n_records': len(same), 'mw': mw,
            })
    return candidates


def _match_with_weight(obs_frags, obs_w, ref_frags, ref_w_map, prec_mz, da_tol=0.2, min_frag_mz=80.0):
    matched_ref = set()
    matched_obs = set()
    pairs = []
    weighted_score = 0.0
    if not ref_frags or not obs_frags:
        return matched_ref, matched_obs, pairs, weighted_score
    # ★ v3opt7 补丁:质量门 — 过滤 < min_frag_mz 的低质量噪声峰
    #   ref 碎片与 obs 碎片同时过滤
    valid_ref = [rf for rf in ref_frags
                 if abs(rf - prec_mz) > da_tol and rf >= min_frag_mz]
    valid_obs = [(j, of) for j, of in enumerate(obs_frags)
                 if abs(of - prec_mz) > da_tol and of >= min_frag_mz]
    if not valid_ref or not valid_obs:
        return matched_ref, matched_obs, pairs, weighted_score
    for rf in valid_ref:
        best_j, best_diff = -1, da_tol
        for j, of in valid_obs:
            if j in matched_obs: continue
            d = abs(of - rf)
            if d <= da_tol and d < best_diff:
                best_diff = d
                best_j = j
        if best_j >= 0:
            matched_ref.add(rf)
            matched_obs.add(best_j)
            w = obs_w[best_j] if best_j < len(obs_w) else 1.0
            rw = ref_w_map.get(round(rf, 4), 1.0)
            pairs.append((obs_frags[best_j], w, rf, rw))
            weighted_score += rw
    return matched_ref, matched_obs, pairs, weighted_score


def identify_v3opt7(extracted_pos, extracted_neg, priority_herbs, db,
                    ppm_tol=20, da_tol=0.2, use_isotope=True,
                    relax_priority=True, rt_tol=0.5):
    """
    v3opt7: v3opt6 + RT 铁证升 L1
    新增: mc ≥ 1 + RT 一致性 ≥ 0.99 + obs 谱图数 ≥ 3 → 可从 L3 升 L1
    """
    import pandas as pd
    records = db['records']
    print(f'  v3opt7 融合生成多加合物候选 ... ', end='')
    fused_groups = {}
    for r in records:
        norm = normalize_formula(r['formula'])
        key = (r['name'], norm)
        if key not in fused_groups:
            fused_groups[key] = []
        fused_groups[key].append(r)
    all_candidates = []
    for (name, norm), group in fused_groups.items():
        cands = fuse_records_v3opt7(group, name, group[0]['formula'])
        all_candidates.extend(cands)
    print(f'{len(fused_groups)} 个化合物 × 多加合物 = {len(all_candidates)} 个候选')

    sorted_by_mz = sorted([(c['mz'], i) for i, c in enumerate(all_candidates) if c['mz'] > 0],
                          key=lambda x: x[0])
    mz_arr = np.array([x[0] for x in sorted_by_mz]) if sorted_by_mz else np.array([])
    compound_lit_count = db['compound_lit_count']
    if priority_herbs:
        priority_set = set(priority_herbs)
    else:
        priority_set = set()

    pos_sorted = sorted([(d['prec'], d) for d in extracted_pos if d['rt'] is not None], key=lambda x: x[0])
    neg_sorted = sorted([(d['prec'], d) for d in extracted_neg if d['rt'] is not None], key=lambda x: x[0])

    results = {}

    def process_one(extracted, mode):
        if not extracted or len(mz_arr) == 0: return
        for prec_data in extracted:
            prec_mz = prec_data['prec']
            obs_frags = prec_data['frags']
            obs_w = prec_data['weights']
            rt = prec_data.get('rt')
            tol = prec_mz * ppm_tol / 1e6
            lo = bisect_left(mz_arr, prec_mz - tol)
            hi = bisect_right(mz_arr, prec_mz + tol)
            for i in range(lo, min(hi, len(mz_arr))):
                mz_val, c_i = sorted_by_mz[i]
                c = all_candidates[c_i]
                if c['mode'] != mode: continue
                ppm_err = abs(prec_mz - mz_val) / mz_val * 1e6
                name = c['name']
                if not name: continue
                is_priority = bool(priority_set) and any(p in c['herb'] for p in priority_set)
                key = f"{name}|{normalize_formula(c['formula'])}"
                ref = c['ref_p'] if mode == 'positive' else c['ref_n']
                ref_w = c['ref_p_weight'] if mode == 'positive' else c['ref_n_weight']
                matched_ref, matched_obs, pairs, ws = _match_with_weight(
                    obs_frags, obs_w, ref, ref_w, prec_mz, da_tol=da_tol)
                mc = len(matched_ref)
                if mc == 0: continue
                if key not in results:
                    results[key] = {
                        'name': name, 'formula': c['formula'],
                        'herb': c['herb'], 'ctype': c['ctype'],
                        'cas': c['cas'], 'lit': c['lit'],
                        'total_lit': compound_lit_count.get(name, c['lit']),
                        'is_priority': is_priority,
                        'ref_pos': c['ref_p'], 'ref_neg': c['ref_n'],
                        'ref_pos_w': c['ref_p_weight'], 'ref_neg_w': c['ref_n_weight'],
                        'matched_ref_pos': set(), 'matched_ref_neg': set(),
                        'matched_obs_pos': set(), 'matched_obs_neg': set(),
                        'pairs_pos': [], 'pairs_neg': [],
                        'obs_pos': [], 'obs_neg': [],
                        'w_pos': [], 'w_neg': [],
                        'rt_pos': [], 'rt_neg': [],
                        'ppm_pos': 999, 'ppm_neg': 999,
                        'mz_pos': None, 'mz_neg': None,
                        'n_records': c['n_records'], 'mw': c['mw'],
                    }
                res = results[key]
                if mode == 'positive':
                    res['matched_ref_pos'].update(matched_ref)
                    res['matched_obs_pos'].update(obs_frags[j] for j in matched_obs)
                    res['pairs_pos'].extend(pairs)
                    res['obs_pos'].extend(obs_frags)
                    res['w_pos'].extend(obs_w)
                    res['rt_pos'].append(rt)
                    if ppm_err < res['ppm_pos']: res['ppm_pos'] = ppm_err
                    res['mz_pos'] = prec_mz
                else:
                    res['matched_ref_neg'].update(matched_ref)
                    res['matched_obs_neg'].update(obs_frags[j] for j in matched_obs)
                    res['pairs_neg'].extend(pairs)
                    res['obs_neg'].extend(obs_frags)
                    res['w_neg'].extend(obs_w)
                    res['rt_neg'].append(rt)
                    if ppm_err < res['ppm_neg']: res['ppm_neg'] = ppm_err
                    res['mz_neg'] = prec_mz

    if extracted_pos: process_one(extracted_pos, 'positive')
    if extracted_neg: process_one(extracted_neg, 'negative')
    if not results: return pd.DataFrame()

    out = []
    for k, r in results.items():
        matched_ref_set = r['matched_ref_pos'] | r['matched_ref_neg']
        matched_obs_set = r['matched_obs_pos'] | r['matched_obs_neg']
        obs = list(set(r['obs_pos'] + r['obs_neg']))
        obs_w = r['w_pos'] + r['w_neg']
        all_ref = r['ref_pos'] + r['ref_neg']
        total_ref = len([rf for rf in all_ref if abs(rf - (r.get('mz_pos') or r.get('mz_neg') or 0)) > 0.2])
        total_obs = len(obs)
        mc = len(matched_ref_set)
        obs_mc = len(matched_obs_set)
        coverage_v2 = mc / total_ref if total_ref > 0 else 0
        coverage_opt = mc / min(total_ref, 10) if total_ref > 0 else 0
        # v4: 覆盖率只作为展示指标,不再进入评级门槛与 base 分数
        if r['mz_pos']: ppm = r['ppm_pos']
        else: ppm = r['ppm_neg']

        # 辅助信号
        iso_score = 0.0
        if use_isotope and (r['mz_pos'] or r['mz_neg']):
            prec_mz = r.get('mz_pos') or r.get('mz_neg')
            iso_score = check_isotope_match(r['formula'], prec_mz, obs, obs_w, tol=0.02)
        rt_consistency_score = 0.0
        all_rts = [x for x in r['rt_pos'] + r['rt_neg'] if x is not None]
        n_obs_unique = len(all_rts)
        if len(all_rts) >= 2:
            rt_med = float(np.median(all_rts))
            n_within = sum(1 for x in all_rts if abs(x - rt_med) <= rt_tol)
            rt_consistency_score = n_within / len(all_rts)
        elif len(all_rts) == 1:
            rt_consistency_score = 1.0
        # ★ v3opt7 补丁:RT 离散度(只有 RT 真的在多个时间点出现过,才算"铁证")
        rt_spread = (max(all_rts) - min(all_rts)) if all_rts else 0.0
        rt_pair_score = 0.0
        if r['mz_pos'] and r['mz_neg']:
            mz_neg = r['mz_neg']
            mz_pos_expected = mz_neg + 2 * 1.007276
            partner_rt, partner_n = find_rt_partner(
                mz_pos_expected, pos_sorted, ppm_tol, rt_tol)
            if partner_n > 0 and r['rt_neg']:
                diffs = [abs(p - float(np.median(r['rt_neg']))) for p in [partner_rt] if p is not None]
                if diffs and diffs[0] < rt_tol:
                    rt_pair_score = 1.0 - diffs[0] / rt_tol
                else:
                    rt_pair_score = 0.3
            elif r['rt_neg']:
                rt_pair_score = 0.5

        is_single_obs = (n_obs_unique == 1)

        # ★★★ 双门棁评级 v3opt7 ★★★
        if r['is_priority']:
            # Priority 化合物: 评级不依赖辅助信号
            if relax_priority:
                if mc >= 2 and ppm <= 50:
                    rating, rating_name = 'I', '确证级'
                elif mc >= 1 and ppm <= 30:
                    rating, rating_name = 'II', '高置信级'
                elif mc >= 1:
                    rating, rating_name = 'III', '推定级'
                else:
                    rating, rating_name = 'V', '排除级'
            else:
                if mc >= 3 and ppm <= 50:
                    rating, rating_name = 'I', '确证级'
                elif mc >= 1 and ppm <= 30:
                    rating, rating_name = 'II', '高置信级'
                elif mc >= 1:
                    rating, rating_name = 'III', '推定级'
                else:
                    rating, rating_name = 'V', '排除级'
            aux_pri = iso_score * 3 + rt_consistency_score * 2 + rt_pair_score * 2
        else:
            # 非 Priority 化合物
            base_rating = None
            if mc >= 4 and ppm <= 15:
                base_rating = ('I', '确证级')
            elif mc >= 3 and ppm <= 20:
                base_rating = ('II', '高置信级')
            elif mc >= 1 and ppm <= 50:
                base_rating = ('III', '推定级')
            elif mc >= 1:
                base_rating = ('IV', '提示级')
            else:
                base_rating = ('V', '排除级')

            # v3opt6 L1 升级门槛
            promote_to_l1 = (
                iso_score >= 0.7 and rt_consistency_score >= 0.97 and mc >= 2 and ppm <= 30
            )
            # v3opt6 L2 升级门槛
            promote_to_l2 = (
                iso_score >= 0.6 and rt_consistency_score >= 0.9 and mc >= 2 and ppm <= 30
            )
            # v3opt6 RT 补救通道
            rt_rescue_to_l2 = (
                iso_score == 0 and rt_consistency_score >= 0.95 and mc >= 2 and ppm <= 30
            )
            # ★ v3opt7 新升级: RT 铁证升 L1
            # mc ≥ 1 + RT 一致性 ≥ 0.99 + obs 谱图数 ≥ 3 + RT 离散度 ≥ 0.1 min → 可从 L3 升 L1
            # 补丁:加 rt_spread 条件,避免单 RT 重复凑数升 L1
            rt_strong_to_l1 = (
                mc >= 1 and ppm <= 30
                and rt_consistency_score >= 0.99
                and n_obs_unique >= 3
                and rt_spread >= 0.1
            )

            if base_rating[0] == 'II' and promote_to_l1:
                rating, rating_name = 'I', '确证级'
            elif base_rating[0] == 'III' and promote_to_l1:
                rating, rating_name = 'I', '确证级'
            elif base_rating[0] == 'III' and rt_strong_to_l1:
                # ★ v3opt7 RT 铁证升 L1
                rating, rating_name = 'I', '确证级'
            elif base_rating[0] == 'III' and promote_to_l2:
                rating, rating_name = 'II', '高置信级'
            elif base_rating[0] == 'III' and rt_rescue_to_l2:
                rating, rating_name = 'II', '高置信级'
            else:
                rating, rating_name = base_rating

            # 单点精化合物 L1 降级
            if rating == 'I' and is_single_obs:
                rating, rating_name = 'II', '高置信级'

            aux_pri = iso_score * 6 + rt_consistency_score * 4 + rt_pair_score * 4

        # 加权相似度
        all_pairs = r['pairs_pos'] + r['pairs_neg']
        n_match = len(all_pairs)
        weighted_sum = sum(p[3] for p in all_pairs)
        cov = n_match / total_ref if total_ref else 0
        rev = n_match / total_obs if total_obs else 0
        sim_composite = (cov * 0.6 + rev * 0.2 + min(cov, 1.0) * 0.2) * (1 + weighted_sum * 0.1)

        # 置信度 (v4: 移除 coverage 加分; 提高 mc / lc 权重,让匹配碎片数 & 文献数主导评分)
        base = 0
        # 匹配碎片数 mc: 阶梯分 0-55(原 0-50,小幅放大)
        if mc >= 15: base += 55
        elif mc >= 10: base += 45
        elif mc >= 6: base += 32
        elif mc >= 4: base += 22
        elif mc >= 2: base += 14
        else: base += 8
        # 文献数 lc: 阶梯分 0-25(原 0-20,小幅放大)
        lc = r['total_lit']
        if lc >= 15: base += 25
        elif lc >= 8: base += 18
        elif lc >= 3: base += 10
        elif lc >= 1: base += 5
        # ppm 精度
        if ppm <= 5: base += 5
        elif ppm <= 15: base += 3
        elif ppm <= 30: base += 1
        # 正负双模式命中
        if r['mz_pos'] and r['mz_neg']: base += 5
        base = min(base, 100)
        base_with_aux = min(base + aux_pri, 100)
        pri_bonus = 30 if r['is_priority'] else 0
        conf = round(base_with_aux + sim_composite * 20 + pri_bonus, 2)

        rts = [x for x in r['rt_pos'] + r['rt_neg'] if x is not None]
        # v4: 覆盖率仅作为参考展示,不再参与评分
        cov_pct = min(coverage_v2 * 100, 100.0)
        cov_opt_pct = min(coverage_opt * 100, 100.0)
        herb_set = [h.strip() for h in re.split(r'[、;/,，\s]+', r['herb']) if h.strip()]
        priority_herbs_in = [h for h in herb_set if any(p in h for p in priority_set)] if priority_set else []
        rt_med_str = f"{float(np.median(rts)):.2f}" if rts else ''
        out.append({
            '化合物中文名': r['name'],
            '分子式': r['formula'],
            'CAS号': r['cas'],
            'm/z实际值': round(r.get('mz_pos') or r.get('mz_neg') or 0, 4),
            'ppm': round(r['ppm_pos'] if r['mz_pos'] else r['ppm_neg'], 4),
            '匹配观测碎片数': obs_mc,
            '匹配参考碎片数': mc,
            '参考碎片总数': total_ref,
            '碎片覆盖率': f"{cov_pct:.1f}%",
            '碎片覆盖率(优化)': f"{cov_opt_pct:.1f}%",
            '谱图相似度': round(sim_composite, 4),
            'ref加权得分': round(weighted_sum, 3),
            '同位素得分': round(iso_score, 3),
            'RT一致性': round(rt_consistency_score, 3),
            '正负RT配对得分': round(rt_pair_score, 3),
            '辅助信号加分': round(aux_pri, 2),
            'obs谱图数': n_obs_unique,
            '单点精化合物': '是' if is_single_obs else '否',
            'RT中位数/min': rt_med_str,
            '出峰时间t/min': round(np.mean(rts), 3) if rts else '',
            '药材匹配': '是' if r['is_priority'] else '否',
            '主要碎片离子': '; '.join([f'{x:.3f}' for x in sorted(obs)[:20]]),
            '参考碎片离子': '; '.join([f'{x:.3f}' for x in sorted(all_ref)[:20]]),
            '库中文献数': lc,
            '评级': rating, '评级名称': rating_name,
            'MSI_Level': 1 if rating == 'I' else (2 if rating == 'II' else 3),
            '置信度': conf,
            '药材来源': r['herb'],
            '药材来源拆分': '; '.join(herb_set) if herb_set else r['herb'],
            'priority_药材命中': '; '.join(priority_herbs_in) if priority_herbs_in else '',
            '化合物类型': r['ctype'],
            '融合记录数': r.get('n_records', 1),
        })
    df = pd.DataFrame(out)
    if len(df) > 0:
        df = df.sort_values(['药材匹配', '置信度'], ascending=[False, False]).reset_index(drop=True)
        df.insert(0, '序号', range(1, len(df) + 1))
    return df
