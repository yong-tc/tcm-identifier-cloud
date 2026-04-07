# -*- coding: utf-8 -*-
"""
中药化合物智能鉴定平台 - 云端部署增强版 v5.15（跨行合并版）
===========================================================
新增功能（v5.15）：
- 自动合并同一化合物在不同行中的碎片离子和文献信息
- 鉴定时使用所有碎片并集进行匹配，文献支持数统计更准确
- 保留 v5.14 所有功能（碎片-文献映射、评分系统等）
"""

import streamlit as st
import pandas as pd
import numpy as np
import os
import time
import pickle
import re
from datetime import datetime
from io import BytesIO
import tempfile
from bisect import bisect_left, bisect_right
import warnings
warnings.filterwarnings('ignore')


# ============================================================================
# 全局路径常量
# ============================================================================
DB_PATHS = [
    "TCM-SM-MS DB.xlsx",
    "TCM-SM-MS DB（中药小分子化学成分高分辨质谱数据库）.xlsx",
    "data/TCM-SM-MS DB.xlsx",
    "data/TCM-SM-MS DB（中药小分子化学成分高分辨质谱数据库）.xlsx",
    "user_input_files/TCM-SM-MS DB.xlsx",
    "user_input_files/TCM-SM-MS DB（中药小分子化学成分高分辨质谱数据库）.xlsx"
]

DIAGNOSTIC_PATHS = [
    "诊断离子.xlsx",
    "data/诊断离子.xlsx",
    "user_input_files/诊断离子.xlsx"
]

VALID_USERNAME = "ZY"
VALID_PASSWORD = "513513"


# ============================================================================
# 辅助函数：标准化分子式
# ============================================================================
def normalize_formula(formula):
    """标准化分子式（去除Unicode下标）"""
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


# ============================================================================
# 鉴定程序核心代码 v5.15（支持跨行合并碎片-文献）
# ============================================================================

class UltimateGardeniaIdentifier:
    """
    中药化合物鉴定终极版程序 v5.15
    新增：自动合并同一化合物的多行碎片和文献信息
    """

    def __init__(self, database_path, ms_positive_path, ms_negative_path,
                 herb_name=None, config=None, use_parallel=True,
                 rt_tolerance=0.3, loss_tolerance=0.02,
                 external_diagnostic_file=None,
                 rt_fusion_tolerance=0.2,
                 intensity_relative_threshold=1.0,
                 tolerance_type='Da',
                 use_rt_score=True,
                 custom_db_path=None,
                 cache_index=True):
        """初始化鉴定程序"""
        self.config = {
            'gradient_time': 30.0,
            'cid_min': 0.5,
            'cid_max': 1.0,
            'fragment_tolerance': 0.05,
            'fragment_tolerance_ppm': 20,
            'tolerance_ppm': 50,
            'max_candidates': 3,
            'min_fragment_count': 1,
            'min_intensity': 100,
            'ppm_tier1': 10,
            'ppm_tier2': 20,
            'max_ppm': 100,
        }

        if config:
            self.config.update(config)

        self.rt_tolerance = rt_tolerance
        self.loss_tolerance = loss_tolerance
        self.rt_fusion_tolerance = rt_fusion_tolerance
        self.intensity_relative_threshold = intensity_relative_threshold / 100.0
        self.tolerance_type = tolerance_type
        self.use_rt_score = use_rt_score
        self.custom_db_path = custom_db_path
        self.cache_index = cache_index
        self.herb_name = herb_name

        self.use_parallel = use_parallel and os.cpu_count() > 1
        self.num_workers = min(os.cpu_count(), 8)

        print("="*80)
        print("中药化合物鉴定程序 v5.15（跨行合并版）")
        print("="*80)
        print("\n【1/8】正在加载数据库...")
        self.full_database = self._load_data(database_path)
        if custom_db_path and os.path.exists(custom_db_path):
            custom_db = self._load_data(custom_db_path)
            if not custom_db.empty:
                print(f"  加载自定义数据库: {len(custom_db)} 条记录")
                self.full_database = pd.concat([self.full_database, custom_db], ignore_index=True)

        if herb_name:
            print(f"【2/8】正在筛选 {herb_name} 相关数据...")
            self.database = self._filter_by_herb(herb_name)
        else:
            self.database = self.full_database.copy()
            print("【2/8】使用全部数据库进行化合物鉴定")
        print(f"  原始数据库记录数: {len(self.database)}")

        # 新增：合并同一化合物的多行记录（碎片、文献等）
        print("【2.5/8】正在合并同一化合物的多行记录（碎片并集、文献并集）...")
        self.database = self._merge_compound_records(self.database)
        print(f"  合并后化合物数: {len(self.database)}")

        print("【3/8】正在加载质谱数据...")
        self.ms_positive = self._load_data(ms_positive_path) if ms_positive_path else pd.DataFrame()
        self.ms_negative = self._load_data(ms_negative_path) if ms_negative_path else pd.DataFrame()
        print(f"  正离子数据: {len(self.ms_positive)} 条记录")
        print(f"  负离子数据: {len(self.ms_negative)} 条记录")

        print("【4/8】正在构建索引...")
        self._build_or_load_index()

        print("【5/8】正在加载诊断离子库...")
        if external_diagnostic_file is None:
            default_diag_path = find_diagnostic_ion_path()
            if default_diag_path:
                print(f"  自动找到默认诊断离子文件: {default_diag_path}")
                external_diagnostic_file = default_diag_path
        self._build_diagnostic_ion_library(external_diagnostic_file)

        print("【6/8】正在加载辅助数据...")
        self._load_auxiliary_data()

        self.stats = {
            'total_precursors': 0,
            'identified_compounds': 0,
            'scored_candidates': 0
        }

        self._print_initialization_info()
        print("【7/8】初始化完成，准备鉴定")

    # ----------------------------- 新增方法：合并同一化合物的多行记录 ----------------------------------
    def _merge_compound_records(self, df):
        """
        合并同一化合物在不同行中的信息（碎片、文献、药材等）
        返回合并后的DataFrame，每个化合物仅保留一行。
        """
        if df.empty:
            return df

        # 确定合并键：优先使用CAS，若无则用“中文名_英文名_分子式”
        df = df.copy()
        df['_merge_key'] = None
        for idx, row in df.iterrows():
            cas = row.get('CAS', '')
            if pd.notna(cas) and str(cas).strip() and str(cas) != 'nan':
                key = f"CAS_{cas}"
            else:
                name_cn = str(row.get('名称（中文）', ''))
                name_en = str(row.get('名称（英文）', ''))
                formula = normalize_formula(row.get('分子式', ''))
                key = f"NAME_{name_cn}_{name_en}_{formula}"
            df.at[idx, '_merge_key'] = key

        # 按合并键分组
        grouped = df.groupby('_merge_key')
        merged_rows = []

        for key, group in grouped:
            # 基础信息取第一条（大部分字段相同）
            first = group.iloc[0].to_dict()

            # 药材名称：合并去重
            herbs = set()
            for _, r in group.iterrows():
                herb = r.get('药材名称') or r.get('药材名')
                if pd.notna(herb) and str(herb).strip():
                    herbs.add(str(herb).strip())
            first['药材名称'] = '; '.join(sorted(herbs)) if herbs else ''

            # 文献来源：合并去重，并建立全局文献列表
            all_sources = []
            for _, r in group.iterrows():
                src = r.get('文献来源') or r.get('文献')
                if pd.notna(src) and str(src).strip():
                    parts = [s.strip() for s in re.split(r'[;,；，]+', str(src)) if s.strip()]
                    all_sources.extend(parts)
            unique_sources = []
            [unique_sources.append(x) for x in all_sources if x not in unique_sources]  # 保持顺序去重
            first['文献来源'] = '; '.join(unique_sources)

            # 辅助：获取全局文献索引映射（旧索引 -> 新索引）
            # 但碎片解析需要根据每行原有的文献列表重新映射，所以这里先不构建映射，
            # 直接对每个碎片重新解析并合并。

            # 合并正离子碎片（带文献索引）
            pos_fragments_dict = {}  # mz -> set of new lit indices
            for _, r in group.iterrows():
                frag_str = r.get('碎片离子（正）', '')
                if pd.isna(frag_str) or not str(frag_str).strip():
                    continue
                # 获取该行对应的原始文献列表（用于解析索引）
                src_str = r.get('文献来源') or r.get('文献')
                if pd.notna(src_str) and str(src_str).strip():
                    row_sources = [s.strip() for s in re.split(r'[;,；，]+', str(src_str)) if s.strip()]
                else:
                    row_sources = []
                # 解析碎片
                parsed = self._parse_fragments_with_literature(frag_str, row_sources)
                for mz, lit_set in parsed:
                    if mz not in pos_fragments_dict:
                        pos_fragments_dict[mz] = set()
                    # lit_set 中的索引是相对于 row_sources 的，需要映射到全局 unique_sources 中的位置
                    for old_idx in lit_set:
                        if old_idx < len(row_sources):
                            old_lit = row_sources[old_idx]
                            # 在 unique_sources 中查找新索引
                            if old_lit in unique_sources:
                                new_idx = unique_sources.index(old_lit)
                                pos_fragments_dict[mz].add(new_idx)
            # 生成新的碎片离子（正）字符串
            if pos_fragments_dict:
                pos_items = []
                for mz in sorted(pos_fragments_dict.keys()):
                    idxs = sorted(pos_fragments_dict[mz])
                    idx_str = ','.join(str(i) for i in idxs)
                    pos_items.append(f"{mz:.3f}:{idx_str}")
                first['碎片离子（正）'] = '; '.join(pos_items)
            else:
                first['碎片离子（正）'] = ''

            # 合并负离子碎片（同理）
            neg_fragments_dict = {}
            for _, r in group.iterrows():
                frag_str = r.get('碎片离子（负）', '')
                if pd.isna(frag_str) or not str(frag_str).strip():
                    continue
                src_str = r.get('文献来源') or r.get('文献')
                if pd.notna(src_str) and str(src_str).strip():
                    row_sources = [s.strip() for s in re.split(r'[;,；，]+', str(src_str)) if s.strip()]
                else:
                    row_sources = []
                parsed = self._parse_fragments_with_literature(frag_str, row_sources)
                for mz, lit_set in parsed:
                    if mz not in neg_fragments_dict:
                        neg_fragments_dict[mz] = set()
                    for old_idx in lit_set:
                        if old_idx < len(row_sources):
                            old_lit = row_sources[old_idx]
                            if old_lit in unique_sources:
                                new_idx = unique_sources.index(old_lit)
                                neg_fragments_dict[mz].add(new_idx)
            if neg_fragments_dict:
                neg_items = []
                for mz in sorted(neg_fragments_dict.keys()):
                    idxs = sorted(neg_fragments_dict[mz])
                    idx_str = ','.join(str(i) for i in idxs)
                    neg_items.append(f"{mz:.3f}:{idx_str}")
                first['碎片离子（负）'] = '; '.join(neg_items)
            else:
                first['碎片离子（负）'] = ''

            # 合并中性丢失（去重）
            losses_set = set()
            for _, r in group.iterrows():
                loss_str = r.get('中性丢失', '')
                if pd.notna(loss_str) and str(loss_str).strip():
                    for loss in str(loss_str).split(','):
                        loss = loss.strip()
                        try:
                            losses_set.add(float(loss))
                        except:
                            pass
            first['中性丢失'] = ','.join(str(l) for l in sorted(losses_set)) if losses_set else ''

            # 保留时间：取第一个非空
            rt_val = None
            for _, r in group.iterrows():
                rt = r.get('保留时间(min)')
                if pd.notna(rt):
                    rt_val = rt
                    break
            first['保留时间(min)'] = rt_val

            # 准分子离子（正/负）：取第一个非空（通常相同）
            for col in ['准分子离子（正）', '准分子离子（负）']:
                val = None
                for _, r in group.iterrows():
                    v = r.get(col)
                    if pd.notna(v):
                        val = v
                        break
                first[col] = val

            # 其他字段保持第一条的值（如化合物类型、CAS等）
            merged_rows.append(first)

        result_df = pd.DataFrame(merged_rows)
        # 删除辅助列
        if '_merge_key' in result_df.columns:
            result_df.drop(columns=['_merge_key'], inplace=True)
        return result_df

    # ----------------------------- 辅助方法 ----------------------------------
    def _normalize_columns(self, df):
        """标准化列名"""
        column_mapping = {
            '药材名': '药材名称',
            '文献': '文献来源',
            '保留时间(min)': '保留时间',
            '中丢失': '中性丢失'
        }
        for old_name, new_name in column_mapping.items():
            if old_name in df.columns and new_name not in df.columns:
                df.rename(columns={old_name: new_name}, inplace=True)
        return df

    def _load_data(self, filepath):
        if filepath and os.path.exists(filepath):
            try:
                if filepath.endswith('.xlsx'):
                    df = pd.read_excel(filepath)
                elif filepath.endswith('.csv'):
                    df = pd.read_csv(filepath)
                else:
                    return pd.DataFrame()
                df = self._normalize_columns(df)
                return df
            except Exception as e:
                print(f"警告: 无法加载文件 {filepath}: {e}")
                return pd.DataFrame()
        print(f"警告: 文件不存在或路径无效: {filepath}")
        return pd.DataFrame()

    def _filter_by_herb(self, herb_name):
        if herb_name is None:
            return self.full_database.copy()
        herb_col = '药材名称' if '药材名称' in self.full_database.columns else '药材名'
        mask = self.full_database[herb_col].str.contains(herb_name, na=False, case=False)
        filtered_db = self.full_database[mask].copy()
        if len(filtered_db) == 0:
            print(f"警告: 未找到 '{herb_name}' 相关数据，将使用全部数据")
            return self.full_database.copy()
        return filtered_db

    def _get_db_cache_key(self, db_path):
        if not db_path or not os.path.exists(db_path):
            return None
        try:
            mtime = os.path.getmtime(db_path)
            size = os.path.getsize(db_path)
            return f"{db_path}_{mtime}_{size}"
        except:
            return None

    def _build_or_load_index(self):
        cache_file = "index_cache.pkl"
        db_key = self._get_db_cache_key(find_database_path())
        if db_key is None:
            db_key = "unknown"

        if self.cache_index and os.path.exists(cache_file):
            try:
                with open(cache_file, 'rb') as f:
                    cached = pickle.load(f)
                if cached.get('db_key') == db_key:
                    print("  从缓存加载索引...")
                    self.sorted_idx_pos = cached['sorted_idx_pos']
                    self.sorted_idx_neg = cached['sorted_idx_neg']
                    self.mz_values_pos = cached['mz_values_pos']
                    self.mz_values_neg = cached['mz_values_neg']
                    self.db_frag_pos = cached['db_frag_pos']
                    self.db_frag_neg = cached['db_frag_neg']
                    self.compound_info = cached['compound_info']
                    self.rt_values = cached.get('rt_values', {})
                    self.neutral_losses = cached.get('neutral_losses', {})
                    print("  缓存加载成功")
                    return
            except Exception as e:
                print(f"  缓存加载失败，重新构建: {e}")
                st.warning(f"缓存文件损坏，正在重建索引。若频繁出现此提示，可手动删除 '{cache_file}' 文件。")

        self._build_optimized_index()

        if self.cache_index:
            try:
                cached = {
                    'db_key': db_key,
                    'sorted_idx_pos': self.sorted_idx_pos,
                    'sorted_idx_neg': self.sorted_idx_neg,
                    'mz_values_pos': self.mz_values_pos,
                    'mz_values_neg': self.mz_values_neg,
                    'db_frag_pos': self.db_frag_pos,
                    'db_frag_neg': self.db_frag_neg,
                    'compound_info': self.compound_info,
                    'rt_values': self.rt_values,
                    'neutral_losses': self.neutral_losses
                }
                with open(cache_file, 'wb') as f:
                    pickle.dump(cached, f)
                print("  索引已缓存")
            except Exception as e:
                print(f"  缓存保存失败: {e}")

    def _parse_fragments_with_literature(self, fragment_string, sources):
        """
        解析碎片离子字符串，支持文献索引后缀。
        格式示例: "151.003:0,1; 137.024:1,2"
        返回: list of (mz, set(lit_indices))
        """
        if pd.isna(fragment_string) or not fragment_string or str(fragment_string).strip() == '':
            return []
        s = str(fragment_string)
        parts = re.split(r'[、;，,\s]+', s)
        result = []
        for part in parts:
            part = part.strip()
            if not part:
                continue
            if ':' in part:
                mz_str, lit_str = part.split(':', 1)
                try:
                    mz = float(mz_str)
                    lit_indices = set()
                    for token in lit_str.split(','):
                        token = token.strip()
                        if token.isdigit():
                            lit_indices.add(int(token))
                    if lit_indices:
                        result.append((mz, lit_indices))
                    else:
                        # 解析失败则默认属于所有文献
                        result.append((mz, set(range(len(sources)))))
                except ValueError:
                    continue
            else:
                try:
                    mz = float(part)
                    result.append((mz, set(range(len(sources)))))
                except ValueError:
                    continue
        return result

    def _build_optimized_index(self):
        self.sorted_idx_pos = []
        self.sorted_idx_neg = []
        self.mz_values_pos = []
        self.mz_values_neg = []
        self.db_frag_pos = []
        self.db_frag_neg = []
        self.compound_info = {}
        self.rt_values = {}
        self.neutral_losses = {}

        for idx, row in self.database.iterrows():
            mz_pos = row.get('准分子离子（正）', 0)
            mz_neg = row.get('准分子离子（负）', 0)

            cn_name = str(row.get('名称（中文）', ''))
            en_name = str(row.get('名称（英文）', ''))
            formula = str(row.get('分子式', ''))
            if formula and formula != 'nan':
                formula = normalize_formula(formula)
            herb = str(row.get('药材名称') or row.get('药材名') or '')
            compound_type = str(row.get('化合物类型', ''))
            source_str = str(row.get('文献来源') or row.get('文献') or '')
            # 解析文献列表
            if pd.notna(source_str) and source_str:
                sources = [s.strip() for s in re.split(r'[;,；，]+', source_str) if s.strip()]
            else:
                sources = []
            cas = str(row.get('CAS', ''))

            self.compound_info[idx] = {
                'name_cn': cn_name,
                'name_en': en_name,
                'formula': formula,
                'herb': herb,
                'compound_type': compound_type,
                'source': source_str,
                'sources': sources,
                'cas': cas if cas and cas != 'nan' else '',
                'adduct_pos': str(row.get('加合物（正）', '')),
                'adduct_neg': str(row.get('加合物（负）', ''))
            }

            if '保留时间(min)' in row and pd.notna(row['保留时间(min)']):
                try:
                    self.rt_values[idx] = float(row['保留时间(min)'])
                except:
                    pass

            if '中性丢失' in row and pd.notna(row['中性丢失']):
                losses = self._parse_losses(row['中性丢失'])
                if losses:
                    self.neutral_losses[idx] = losses

            # 正离子碎片（带文献）
            if pd.notna(mz_pos) and str(mz_pos).strip() != '':
                try:
                    mz_val = float(mz_pos)
                    if mz_val > 0:
                        fragments_with_lit = self._parse_fragments_with_literature(
                            row.get('碎片离子（正）', ''), sources
                        )
                        self.sorted_idx_pos.append((mz_val, idx, fragments_with_lit))
                except (ValueError, TypeError):
                    pass

            # 负离子碎片（带文献）
            if pd.notna(mz_neg) and str(mz_neg).strip() != '':
                try:
                    mz_val = float(mz_neg)
                    if mz_val > 0:
                        fragments_with_lit = self._parse_fragments_with_literature(
                            row.get('碎片离子（负）', ''), sources
                        )
                        self.sorted_idx_neg.append((mz_val, idx, fragments_with_lit))
                except (ValueError, TypeError):
                    pass

        self.sorted_idx_pos.sort(key=lambda x: x[0])
        self.sorted_idx_neg.sort(key=lambda x: x[0])

        self.mz_values_pos = np.array([x[0] for x in self.sorted_idx_pos])
        self.mz_values_neg = np.array([x[0] for x in self.sorted_idx_neg])

        self.db_frag_pos = [x[2] for x in self.sorted_idx_pos]
        self.db_frag_neg = [x[2] for x in self.sorted_idx_neg]

    def _parse_losses(self, loss_string):
        if pd.isna(loss_string):
            return []
        losses = []
        for part in str(loss_string).split(','):
            part = part.strip()
            try:
                losses.append(float(part))
            except ValueError:
                continue
        return losses

    def _build_diagnostic_ion_library(self, external_file=None):
        if external_file and os.path.exists(external_file):
            try:
                df = pd.read_excel(external_file)
                required_cols = ['化合物类型', '诊断碎片离子m/z']
                if not all(col in df.columns for col in required_cols):
                    raise ValueError(f"外部诊断离子文件缺少必要列：{required_cols}")

                self.diagnostic_ions = {}
                for category, group in df.groupby('化合物类型'):
                    ions_dict = {}
                    for _, row in group.iterrows():
                        mz = row['诊断碎片离子m/z']
                        if pd.notna(mz):
                            try:
                                mz = float(mz)
                            except:
                                continue
                            weight = row.get('权重', 1)
                            try:
                                weight = float(weight)
                            except (ValueError, TypeError):
                                weight = 1
                            ions_dict[mz] = ions_dict.get(mz, 0) + weight
                    description = group['描述'].iloc[0] if '描述' in group.columns else f'来自外部文件，{len(ions_dict)}个离子'
                    self.diagnostic_ions[category] = {
                        'ions': list(ions_dict.keys()),
                        'weights': list(ions_dict.values()),
                        'description': description
                    }
                print(f"  成功加载外部诊断离子库：{len(self.diagnostic_ions)} 类，{len(df)} 条记录（去重后）")
            except Exception as e:
                print(f"  加载外部诊断离子文件失败：{e}，将使用内置库")
                self._build_default_diagnostic_ions()
        else:
            print("  未找到外部诊断离子文件，使用内置诊断离子库")
            self._build_default_diagnostic_ions()

    def _build_default_diagnostic_ions(self):
        self.diagnostic_ions = {
            '环烯醚萜类': {
                'ions': [138.055, 124.039, 110.023, 96.008, 82.029, 67.029, 127.039],
                'weights': [1]*7,
                'description': '环烯醚萜类特征脱水碎片'
            },
            '有机酸类': {
                'ions': [191.056, 179.034, 173.045, 135.045, 93.034, 85.029],
                'weights': [1]*6,
                'description': '咖啡酰奎尼酸系列特征离子'
            },
            '黄酮类': {
                'ions': [151.003, 137.024, 121.029, 107.049, 81.034, 65.039],
                'weights': [1]*6,
                'description': '黄酮苷元特征碎片'
            },
            '萜类': {
                'ions': [127.076, 113.060, 99.044, 85.029, 71.013],
                'weights': [1]*5,
                'description': '萜类特征碎片'
            },
            '栀子特异': {
                'ions': [127.039, 113.024, 101.024, 69.034, 97.028],
                'weights': [1]*5,
                'description': '栀子特有成分特征离子'
            },
            '生物碱类': {
                'ions': [105.070, 91.054, 79.054, 65.039],
                'weights': [1]*4,
                'description': '生物碱特征碎片'
            },
            '酚酸类': {
                'ions': [137.024, 123.044, 109.028, 95.049],
                'weights': [1]*4,
                'description': '酚酸类特征碎片'
            }
        }

    def _load_auxiliary_data(self):
        pass

    def _find_neutral_losses(self, fragments, precursor_mz):
        if len(fragments) == 0:
            return []
        observed_losses = []
        precursor_mz = float(precursor_mz)
        for frag in fragments:
            loss = precursor_mz - float(frag)
            if 10 < loss < 500:
                observed_losses.append(round(loss, 4))
        return observed_losses

    def _classify_compound(self, name, compound_type):
        name = name.lower()
        compound_type = str(compound_type).lower()

        if any(keyword in compound_type for keyword in ['环烯醚', 'iridoid']):
            return '环烯醚萜类'
        if any(keyword in compound_type for keyword in ['有机酸', '酚酸']):
            return '有机酸类'
        if any(keyword in compound_type for keyword in ['黄酮', 'flavon']):
            return '黄酮类'
        if any(keyword in compound_type for keyword in ['萜', 'terpen']):
            return '萜类'
        if any(keyword in compound_type for keyword in ['生物碱', 'alkaloid']):
            return '生物碱类'
        if any(keyword in name for keyword in ['京尼平', '栀子', 'genipos', 'gardenia']):
            return '栀子特异'
        return '其他'

    def _binary_search_range(self, mz_array, mz, tolerance_ppm):
        if len(mz_array) == 0:
            return range(0, 0)
        mz = float(mz)
        tolerance = mz * tolerance_ppm / 1e6
        mz_min = mz - tolerance
        mz_max = mz + tolerance
        left = bisect_left(mz_array, mz_min)
        right = bisect_right(mz_array, mz_max)
        return range(left, right)

    def _match_fragments_fast(self, observed, reference_frags, tolerance_value, tolerance_type='Da', precursor_mz=None):
        """
        reference_frags: list of (mz, lit_indices)
        observed: array of observed mz values
        返回: (matched_mz_list, matched_lit_indices)
        """
        if len(reference_frags) == 0 or len(observed) == 0:
            return [], set()
        observed_sorted = np.sort(observed)
        matched_mz = []
        matched_lit = set()
        for ref_mz, lit_set in reference_frags:
            if pd.isna(ref_mz) or ref_mz <= 0:
                continue
            if tolerance_type == 'Da':
                tol = tolerance_value
                left = np.searchsorted(observed_sorted, ref_mz - tol)
                right = np.searchsorted(observed_sorted, ref_mz + tol)
                for i in range(left, right):
                    obs_val = observed_sorted[i]
                    if precursor_mz is not None and abs(obs_val - precursor_mz) <= tolerance_value:
                        continue
                    matched_mz.append(obs_val)
                    matched_lit.update(lit_set)
                    break
            else:  # ppm
                tol = ref_mz * tolerance_value / 1e6
                left = np.searchsorted(observed_sorted, ref_mz - tol)
                right = np.searchsorted(observed_sorted, ref_mz + tol)
                for i in range(left, right):
                    obs_val = observed_sorted[i]
                    if precursor_mz is not None and abs(obs_val - precursor_mz) <= tolerance_value:
                        continue
                    matched_mz.append(obs_val)
                    matched_lit.update(lit_set)
                    break
        return list(set(matched_mz)), matched_lit

    def _find_diagnostic_ions_fast(self, matched_fragments, category, precursor_mz=None):
        if len(matched_fragments) == 0:
            return [], []
        if category not in self.diagnostic_ions:
            return [], []

        diag_data = self.diagnostic_ions[category]
        diag_ions = np.array(diag_data['ions'])
        diag_weights = diag_data['weights']
        matched_arr = np.asarray(matched_fragments)
        tolerance = self.config['fragment_tolerance']

        diagnostic = []
        weights_used = []
        for diag_val, weight in zip(diag_ions, diag_weights):
            for matched_val in matched_arr:
                if abs(matched_val - diag_val) <= tolerance:
                    if precursor_mz is not None and abs(matched_val - precursor_mz) <= tolerance:
                        continue
                    diagnostic.append(float(matched_val))
                    weights_used.append(weight)
                    break
        return diagnostic, weights_used

    def _determine_confidence_level(self, ppm, matched_count, diagnostic_count, has_fragment_data):
        ppm_tier1 = self.config.get('ppm_tier1', 10)
        ppm_tier2 = self.config.get('ppm_tier2', 20)
        max_ppm = self.config.get('max_ppm', 100)

        if ppm > max_ppm:
            return 6, '排除级', '未识别', f'ppm > {max_ppm}，不符合报告要求'
        if ppm <= 5 and matched_count >= 3:
            return 1, '确证级', '最高', '高精度+丰富碎片，可直接报告'
        if ppm <= ppm_tier1 and matched_count >= 2 and diagnostic_count >= 1:
            return 1, '确证级', '最高', '可直接报告'
        if ppm <= ppm_tier2 and matched_count >= 2:
            return 2, '高置信级', '良好', '建议复核后报告'
        if ppm <= 50 and matched_count >= 1:
            return 3, '推定级', '中等', '需验证后报告'
        if ppm <= 50 and not has_fragment_data:
            return 4, '提示级', '较低', '仅供筛查'
        if ppm <= 50 and has_fragment_data and matched_count == 0:
            return 5, '参考级', '受限', '信息受限'
        if 50 < ppm <= 75:
            return 4, '提示级', '较低', 'ppm稍大，仅供筛查参考'
        if 75 < ppm <= max_ppm:
            return 5, '参考级', '受限', 'ppm较大，需谨慎使用'
        return 6, '排除级', '未识别', '不符合评级标准'

    def _calculate_base_score(self, rating, ppm, matched_frag_count, diag_count,
                              diag_weights=None, neutral_loss_match_count=0, 
                              rt_deviation=None, lit_match_count=0):
        if rating == 1:
            base = 85
        elif rating == 2:
            base = 65
        elif rating == 3:
            base = 45
        elif rating == 4:
            base = 25
        else:
            base = 0

        if ppm <= 5:
            ppm_adj = 5
        elif ppm <= 10:
            ppm_adj = 0
        elif ppm <= 20:
            ppm_adj = -5
        elif ppm <= 30:
            ppm_adj = -10
        else:
            ppm_adj = -15

        frag_adj = min(matched_frag_count * 2, 20)

        if diag_weights:
            diag_score = min(sum(diag_weights) * 5, 15)
        else:
            diag_score = min(diag_count * 5, 15)

        loss_adj = min(neutral_loss_match_count * 2, 10)

        rt_adj = 0
        if rt_deviation is not None and self.use_rt_score:
            if rt_deviation < 0.2:
                rt_adj = 5
            elif rt_deviation < 0.5:
                rt_adj = 2

        lit_adj = min(lit_match_count * 3, 12)   # 每匹配一篇文献 +3 分，上限 12

        total = base + ppm_adj + frag_adj + diag_score + loss_adj + rt_adj + lit_adj
        total = max(0, total)

        if rating == 1 and total < 80:
            total = 80
        if rating == 2 and total < 60:
            total = 60

        total = min(total, 100)
        return round(total, 2)

    def extract_precursor_ions(self, ms_data, ionization_mode):
        if ms_data.empty:
            return []

        precursor_col = None
        possible_names = ['Precursor M/z', 'Precursor_mz', 'Precursor', 'precursor m/z', 'precursor_mz']
        for col in ms_data.columns:
            if col in possible_names:
                precursor_col = col
                break
        if precursor_col is None:
            for col in ms_data.columns:
                if 'precursor' in col.lower() and 'm/z' in col.lower():
                    precursor_col = col
                    break
        if precursor_col is None:
            st.error("未找到母离子列（Precursor M/z），请检查文件格式或使用自定义列名功能。")
            return []

        mz_columns = []
        for col in ms_data.columns:
            if ('Peak_' in col and '_m/z' in col) or ('m/z' in col.lower() and col != precursor_col):
                mz_columns.append(col)
        if not mz_columns:
            st.warning("未找到碎片离子列，将只进行一级匹配。")

        min_intensity_abs = self.config['min_intensity']
        rel_threshold = self.intensity_relative_threshold

        precursors = []
        for idx, row in ms_data.iterrows():
            precursor_mz = row.get(precursor_col)
            if pd.notna(precursor_mz) and float(precursor_mz) > 0:
                fragments = []
                fragments_dict = {}
                base_peak = 0
                for col in mz_columns:
                    intensity = 0
                    intensity_col = col.replace('_m/z', '_Intensity') if '_m/z' in col else None
                    if intensity_col and intensity_col in row.index:
                        intensity = float(row[intensity_col]) if pd.notna(row[intensity_col]) else 0
                    else:
                        intensity = 1
                    if intensity > base_peak:
                        base_peak = intensity

                for col in mz_columns:
                    fragment_mz = row[col]
                    if pd.notna(fragment_mz) and float(fragment_mz) > 0:
                        intensity = 0
                        intensity_col = col.replace('_m/z', '_Intensity') if '_m/z' in col else None
                        if intensity_col and intensity_col in row.index:
                            intensity = float(row[intensity_col]) if pd.notna(row[intensity_col]) else 0
                        else:
                            intensity = 1

                        if intensity >= min_intensity_abs:
                            if base_peak > 0:
                                rel_intensity = intensity / base_peak
                                if rel_intensity >= rel_threshold:
                                    mz_value = float(fragment_mz)
                                    fragments.append(mz_value)
                                    fragments_dict[mz_value] = intensity
                            else:
                                mz_value = float(fragment_mz)
                                fragments.append(mz_value)
                                fragments_dict[mz_value] = intensity

                rt = row.get('出峰时间t/min', np.nan)
                if pd.isna(rt):
                    rt = row.get('出峰时间', np.nan)

                if pd.isna(rt) and 'CID' in row.index:
                    cid = row['CID']
                    gt = self.config['gradient_time']
                    cid_min = self.config['cid_min']
                    cid_max = self.config['cid_max']
                    if pd.isna(cid):
                        rt = 0.5
                    else:
                        if cid >= cid_max:
                            rt = gt
                        else:
                            rt = round(gt * (cid - cid_min) / (cid_max - cid_min), 2)
                            rt = max(0, min(gt, rt))

                precursors.append({
                    'precursor_mz': float(precursor_mz),
                    'retention_time': rt if not pd.isna(rt) else None,
                    'fragments': np.array(sorted(fragments, reverse=True)),
                    'fragments_dict': fragments_dict,
                    'mode': ionization_mode,
                    'total_fragments': len(fragments),
                    'base_peak': base_peak
                })

        return precursors

    def identify_compound(self, precursor, return_top=1, score_threshold=0):
        mode = precursor['mode']
        if mode == '正离子':
            sorted_idx = self.sorted_idx_pos
            mz_array = self.mz_values_pos
            db_frags = self.db_frag_pos
        else:
            sorted_idx = self.sorted_idx_neg
            mz_array = self.mz_values_neg
            db_frags = self.db_frag_neg

        candidate_range = self._binary_search_range(mz_array, precursor['precursor_mz'], self.config['tolerance_ppm'])

        scored_candidates = []
        for idx in candidate_range:
            mz_val, db_idx, ref_frags = sorted_idx[idx]
            if db_idx not in self.compound_info:
                continue

            info = self.compound_info[db_idx]
            category = self._classify_compound(info['name_cn'] + ' ' + info['name_en'], info['compound_type'])

            if self.tolerance_type == 'Da':
                tolerance_value = self.config['fragment_tolerance']
            else:
                tolerance_value = self.config['fragment_tolerance_ppm']

            matched_fragments, matched_lit_indices = self._match_fragments_fast(
                precursor['fragments'],
                ref_frags,
                tolerance_value,
                self.tolerance_type,
                precursor['precursor_mz']
            )

            diagnostic_ions, diag_weights = self._find_diagnostic_ions_fast(
                matched_fragments,
                category,
                precursor['precursor_mz']
            )

            neutral_loss_matches = 0
            if db_idx in self.neutral_losses:
                expected_losses = self.neutral_losses[db_idx]
                observed_losses = self._find_neutral_losses(precursor['fragments'], precursor['precursor_mz'])
                for exp in expected_losses:
                    for obs in observed_losses:
                        if abs(exp - obs) <= self.loss_tolerance:
                            neutral_loss_matches += 1
                            break

            theoretical_mz = mz_val
            ppm = abs(float(precursor['precursor_mz']) - theoretical_mz) / theoretical_mz * 1e6

            rt_deviation = None
            if self.use_rt_score and db_idx in self.rt_values and precursor['retention_time'] is not None:
                rt_deviation = abs(precursor['retention_time'] - self.rt_values[db_idx])

            total_sources = len(info.get('sources', []))
            matched_lit_count = len(matched_lit_indices) if total_sources > 0 else 0

            candidate = {
                'name_cn': info['name_cn'],
                'name_en': info['name_en'],
                'formula': info['formula'],
                'cas': info['cas'],
                'herb': info['herb'],
                'compound_type': info['compound_type'],
                'observed_mz': precursor['precursor_mz'],
                'theoretical_mz': theoretical_mz,
                'ppm': ppm,
                'adduct': info['adduct_pos'] if mode == '正离子' else info['adduct_neg'],
                'matched_fragments': matched_fragments,
                'matched_lit_indices': matched_lit_indices,
                'matched_lit_count': matched_lit_count,
                'diagnostic_ions': diagnostic_ions,
                'diag_weights': diag_weights,
                'mode': mode,
                'source': info['source'],
                'sources': info.get('sources', []),
                'category': category,
                'db_index': db_idx,
                'neutral_loss_matches': neutral_loss_matches,
                'rt_deviation': rt_deviation
            }

            candidate['temp_score'] = -ppm
            scored_candidates.append(candidate)

        scored_candidates.sort(key=lambda x: (-len(x['matched_fragments']), -x['temp_score']))

        if return_top == 1:
            return scored_candidates[0] if scored_candidates else None
        else:
            return scored_candidates[:return_top]

    def generate_report(self, herb_name=None):
        if herb_name is None:
            herb_name = self.herb_name if self.herb_name else '中药'

        best_records = {}
        herbs_collection = {}

        pos_precursors = self.extract_precursor_ions(self.ms_positive, '正离子') if not self.ms_positive.empty else []
        neg_precursors = self.extract_precursor_ions(self.ms_negative, '负离子') if not self.ms_negative.empty else []

        all_precursors = pos_precursors + neg_precursors
        total = len(all_precursors)

        if total == 0:
            print("警告: 没有提取到任何母离子，请检查质谱数据文件格式。")
            return pd.DataFrame()

        for i, precursor in enumerate(all_precursors):
            self.stats['total_precursors'] += 1
            if (i + 1) % 500 == 0:
                print(f"    处理进度: {i+1}/{total}")

            best_candidate = self.identify_compound(precursor, return_top=1)
            if best_candidate is None:
                continue

            self.stats['identified_compounds'] += 1
            self.stats['scored_candidates'] += 1

            lvl_id, lvl_name, confidence, suggestion = self._determine_confidence_level(
                best_candidate['ppm'],
                len(best_candidate['matched_fragments']),
                len(best_candidate['diagnostic_ions']),
                len(precursor['fragments']) > 0
            )

            if lvl_id == 6:
                continue

            formula = best_candidate['formula'] if best_candidate['formula'] and best_candidate['formula'] != 'nan' else '待确定'
            en_name = best_candidate['name_en'] if best_candidate['name_en'] and best_candidate['name_en'] != 'nan' else ''

            source_str = best_candidate.get('source', '')
            if pd.notna(source_str) and source_str:
                source_list = [s.strip() for s in re.split(r'[;,；，]+', str(source_str)) if s.strip()]
                source_count = len(source_list)
                source_display = '; '.join(source_list)
            else:
                source_count = 0
                source_display = ''

            matched_frags = best_candidate['matched_fragments']
            diag_ions = best_candidate['diagnostic_ions']
            diag_weights = best_candidate['diag_weights']

            adduct = best_candidate['adduct']
            if pd.isna(adduct) or adduct == 'nan' or adduct == '':
                adduct = ''
            else:
                for sep in [',', ';']:
                    if sep in adduct:
                        adduct = adduct.split(sep)[0].strip()
                        break

            base_score = self._calculate_base_score(
                lvl_id,
                best_candidate['ppm'],
                len(matched_frags),
                len(diag_ions),
                diag_weights,
                best_candidate['neutral_loss_matches'],
                best_candidate['rt_deviation'],
                best_candidate['matched_lit_count']   # 新增文献匹配数
            )

            cas = best_candidate.get('cas', '')
            if cas and cas != 'nan' and cas != '':
                merge_key = cas
            elif en_name and en_name != '':
                merge_key = f"{en_name}_{formula}"
            else:
                merge_key = f"{best_candidate['name_cn']}_{formula}"

            record = {
                '序号': 0,
                '出峰时间t/min': precursor['retention_time'],
                '化合物中文名': best_candidate['name_cn'],
                '化合物英文名': en_name,
                '分子式': formula,
                'CAS号': best_candidate['cas'],
                '药材名称': best_candidate['herb'],
                '化合物类型': best_candidate['compound_type'],
                '离子化方式': best_candidate['mode'],
                '加和离子': adduct,
                'm/z实际值': round(best_candidate['observed_mz'], 4),
                'm/z理论值': round(best_candidate['theoretical_mz'], 4),
                'ppm': round(best_candidate['ppm'], 4),
                '是否有碎片数据': '是' if precursor['fragments'].size > 0 else '否',
                '主要碎片离子': '; '.join([f'{f:.4f}' for f in matched_frags]) if matched_frags else '',
                '匹配碎片数': len(matched_frags),
                '匹配文献数': best_candidate['matched_lit_count'],   # 新增字段
                '诊断性离子个数': len(diag_ions),
                '诊断性离子': '; '.join([f'{f:.4f}' for f in diag_ions]) if diag_ions else '',
                '文献来源数': source_count,
                '文献来源': source_display,
                '评级': lvl_id,
                '评级名称': lvl_name,
                '置信度': confidence,
                '报告建议': suggestion,
                '基础得分': base_score,
                '综合得分': base_score,
                '_merge_key': merge_key,
                '_fragments_set': set(matched_frags)
            }

            if merge_key not in herbs_collection:
                herbs_collection[merge_key] = set()
            herbs_collection[merge_key].add(best_candidate['herb'])

            current_score = base_score
            if merge_key not in best_records or current_score > best_records[merge_key][1]:
                best_records[merge_key] = (record, current_score)

        records_list = []
        for merge_key, (record, score) in best_records.items():
            herbs = herbs_collection[merge_key]
            record['药材名称'] = '; '.join(sorted(herbs))
            records_list.append(record)

        if not records_list:
            return pd.DataFrame()

        records_list.sort(key=lambda x: x['出峰时间t/min'] if x['出峰时间t/min'] is not None else 0)

        rt_tol = self.rt_fusion_tolerance
        fused_records = []
        i = 0
        n = len(records_list)
        while i < n:
            group = [records_list[i]]
            j = i + 1
            while j < n:
                rt1 = group[0]['出峰时间t/min']
                rt2 = records_list[j]['出峰时间t/min']
                if rt1 is None or rt2 is None:
                    break
                if abs(rt1 - rt2) <= rt_tol:
                    group.append(records_list[j])
                    j += 1
                else:
                    break

            compound_groups = {}
            for rec in group:
                key = rec['_merge_key']
                if key not in compound_groups:
                    compound_groups[key] = []
                compound_groups[key].append(rec)

            for comp_key, comp_recs in compound_groups.items():
                if len(comp_recs) == 1:
                    rec = comp_recs[0].copy()
                    del rec['_merge_key']
                    del rec['_fragments_set']
                    if '基础得分' in rec:
                        del rec['基础得分']
                    fused_records.append(rec)
                else:
                    best_rec = max(comp_recs, key=lambda x: x['基础得分']).copy()
                    adducts = set()
                    modes = set()
                    for rec in comp_recs:
                        if rec['加和离子']:
                            adducts.add(rec['加和离子'])
                        if rec['离子化方式']:
                            modes.add(rec['离子化方式'])
                    best_rec['加和离子'] = '; '.join(sorted(adducts))
                    best_rec['离子化方式'] = '/'.join(sorted(modes))

                    all_frag_sets = [rec['_fragments_set'] for rec in comp_recs]
                    similarities = []
                    for a in range(len(all_frag_sets)):
                        for b in range(a+1, len(all_frag_sets)):
                            set_a = all_frag_sets[a]
                            set_b = all_frag_sets[b]
                            if len(set_a) == 0 and len(set_b) == 0:
                                sim = 1.0
                            elif len(set_a) == 0 or len(set_b) == 0:
                                sim = 0.0
                            else:
                                intersection = len(set_a & set_b)
                                union = len(set_a | set_b)
                                sim = intersection / union if union > 0 else 0
                            similarities.append(sim)
                    avg_similarity = np.mean(similarities) if similarities else 1.0

                    extra_count = len(adducts) - 1
                    if extra_count > 0:
                        base_bonus = extra_count * 5
                        if avg_similarity < 0.5:
                            fusion_bonus = base_bonus * 0.5
                        else:
                            fusion_bonus = base_bonus
                        fusion_bonus = min(fusion_bonus, 15)
                    else:
                        fusion_bonus = 0

                    final_score = min(best_rec['基础得分'] + fusion_bonus, 100)

                    if best_rec['评级'] == 1 and final_score < 80:
                        final_score = 80
                    if best_rec['评级'] == 2 and final_score < 60:
                        final_score = 60

                    best_rec['综合得分'] = final_score

                    if len(comp_recs) >= 2 and all(rec['评级'] <= 2 for rec in comp_recs):
                        best_rec['评级'] = 1
                        best_rec['评级名称'] = '确证级'
                        best_rec['置信度'] = '最高'
                        best_rec['报告建议'] = '可直接报告'

                    del best_rec['_merge_key']
                    del best_rec['_fragments_set']
                    if '基础得分' in best_rec:
                        del best_rec['基础得分']
                    fused_records.append(best_rec)

            i = j

        report_df = pd.DataFrame(fused_records)

        if not report_df.empty:
            report_df = report_df.sort_values(by=['评级', '综合得分', 'ppm'], ascending=[True, False, True])
            report_df = report_df.reset_index(drop=True)
            report_df['序号'] = range(1, len(report_df) + 1)

        return report_df

    def save_report(self, report_df, output_path):
        report_df.to_excel(output_path, index=False)
        print(f"\n报告已保存至: {output_path}")

    def print_summary(self, report_df, herb_name):
        print("\n" + "="*100)
        print(f"【{herb_name}化合物鉴定报告】")
        print("="*100)
        print(f"鉴定药材: {herb_name}")
        print(f"鉴定时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"共鉴定出 {len(report_df)} 个化合物")

        print(f"\n【处理统计】")
        print(f"  - 总母离子数: {self.stats['total_precursors']}")
        print(f"  - 鉴定化合物数: {self.stats['identified_compounds']}")
        print(f"  - 数据库记录数: {len(self.database)}")

        if not report_df.empty:
            print(f"\n【评级分布】")
            for level, count in report_df['评级名称'].value_counts().items():
                print(f"  - {level}: {count} 个")

            print(f"\n【ppm误差分布】")
            ppm_10 = (report_df['ppm'] <= 10).sum()
            ppm_20 = ((report_df['ppm'] > 10) & (report_df['ppm'] <= 20)).sum()
            ppm_50 = ((report_df['ppm'] > 20) & (report_df['ppm'] <= 50)).sum()
            print(f"  - ≤10ppm: {ppm_10} 个")
            print(f"  - 10-20ppm: {ppm_20} 个")
            print(f"  - 20-50ppm: {ppm_50} 个")

            print(f"\n【综合得分分布】")
            score_90_100 = (report_df['综合得分'] >= 90).sum()
            score_80_89 = ((report_df['综合得分'] >= 80) & (report_df['综合得分'] < 90)).sum()
            score_70_79 = ((report_df['综合得分'] >= 70) & (report_df['综合得分'] < 80)).sum()
            score_60_69 = ((report_df['综合得分'] >= 60) & (report_df['综合得分'] < 70)).sum()
            score_below_60 = (report_df['综合得分'] < 60).sum()
            print(f"  - 90-100分: {score_90_100} 个")
            print(f"  - 80-89分: {score_80_89} 个")
            print(f"  - 70-79分: {score_70_79} 个")
            print(f"  - 60-69分: {score_60_69} 个")
            print(f"  - <60分: {score_below_60} 个")

            print(f"\n【药材来源分布（前10）】")
            for herb, count in report_df['药材名称'].value_counts().head(10).items():
                print(f"  - {herb}: {count} 个")

        print("\n" + "="*100)

    def _print_initialization_info(self):
        print("\n" + "="*80)
        print("程序初始化完成（跨行合并版）")
        print("="*80)
        print(f"  - 数据库记录数: {len(self.database)} 条")
        print(f"  - 正离子索引: {len(self.mz_values_pos)} 条")
        print(f"  - 负离子索引: {len(self.mz_values_neg)} 条")
        print(f"  - 诊断性离子库: {len(self.diagnostic_ions)} 类")
        print(f"  - 保留时间数据: {len(self.rt_values)} 条")
        print(f"  - 中性丢失数据: {len(self.neutral_losses)} 条")
        print(f"  - 容差类型: {self.tolerance_type}")
        print(f"  - 强度相对阈值: {self.intensity_relative_threshold*100:.1f}%")
        print(f"  - RT得分: {'启用' if self.use_rt_score else '禁用'}")
        print(f"  - 并行处理: {'启用' if self.use_parallel else '禁用'}")
        print("="*80)


# ============================================================================
# 数据库加载函数（带缓存）
# ============================================================================

@st.cache_data
def load_database_cached(db_filename=None):
    if db_filename and os.path.exists(db_filename):
        try:
            df = pd.read_excel(db_filename)
            return df
        except:
            pass
    for path in DB_PATHS:
        if os.path.exists(path):
            try:
                df = pd.read_excel(path)
                return df
            except Exception as e:
                continue
    return pd.DataFrame()


def find_database_path():
    for path in DB_PATHS:
        if os.path.exists(path):
            return path
    return None


@st.cache_data
def load_diagnostic_ions_cached():
    for path in DIAGNOSTIC_PATHS:
        if os.path.exists(path):
            try:
                df = pd.read_excel(path)
                if '诊断碎片离子m/z' in df.columns:
                    df['诊断碎片离子m/z'] = pd.to_numeric(df['诊断碎片离子m/z'], errors='coerce')
                return df
            except Exception as e:
                st.warning(f"加载诊断离子文件时出错: {e}")
                continue
    return pd.DataFrame()


def find_diagnostic_ion_path():
    for path in DIAGNOSTIC_PATHS:
        if os.path.exists(path):
            return path
    return None


def match_diagnostic_ions(user_mz_values, diagnostic_df, tolerance_ppm=10, ion_mode=None):
    if diagnostic_df.empty or not user_mz_values:
        return pd.DataFrame()
    if ion_mode and ion_mode != "全部":
        filtered_df = diagnostic_df[diagnostic_df['离子模式'] == ion_mode].copy()
    else:
        filtered_df = diagnostic_df.copy()
    if '诊断碎片离子m/z' not in filtered_df.columns:
        return pd.DataFrame()
    filtered_df = filtered_df.dropna(subset=['诊断碎片离子m/z'])
    if filtered_df.empty:
        return pd.DataFrame()
    results = []
    for user_mz in user_mz_values:
        user_mz = float(user_mz)
        tolerance = user_mz * tolerance_ppm / 1e6
        mz_min = user_mz - tolerance
        mz_max = user_mz + tolerance
        matches = filtered_df[
            (filtered_df['诊断碎片离子m/z'] >= mz_min) &
            (filtered_df['诊断碎片离子m/z'] <= mz_max)
        ]
        for _, row in matches.iterrows():
            ref_mz = row['诊断碎片离子m/z']
            ppm_error = abs(user_mz - ref_mz) / ref_mz * 1e6
            results.append({
                '输入m/z': user_mz,
                '匹配诊断离子m/z': ref_mz,
                '误差(ppm)': round(ppm_error, 4),
                '化合物类型': row.get('化合物类型', ''),
                '离子模式': row.get('离子模式', ''),
                '中文名称': row.get('中文名称', ''),
                '英文名称': row.get('英文名称', ''),
                '分子式': row.get('分子式', ''),
                '药材名': row.get('药材名', ''),
                '准分子离子m/z': row.get('准分子离子m/z', ''),
                '加合形式': row.get('加合形式', ''),
                '相对丰度': row.get('相对丰度', ''),
                '类特征性离子': row.get('类特征性离子', False)
            })
    return pd.DataFrame(results)


# ============================================================================
# Streamlit 网页应用部分
# ============================================================================

st.set_page_config(
    page_title="中药化合物智能鉴定平台 v5.15",
    page_icon="🌿",
    layout="wide",
    initial_sidebar_state="expanded"
)


def load_optimized_css():
    st.markdown("""
    <style>
        .stApp { background: linear-gradient(135deg, #f8fafc 0%, #e2e8f0 50%, #f1f5f9 100%); }
        .main-header { background: linear-gradient(135deg, #059669 0%, #0891b2 50%, #7c3aed 100%); padding: 2rem; border-radius: 20px; margin-bottom: 2rem; }
        .main-header h1, .main-header p { color: white; }
        .stat-card { background: rgba(255,255,255,0.95); border-radius: 16px; padding: 1.5rem; text-align: center; }
        .feature-card { background: white; border-radius: 16px; padding: 1.5rem; margin-bottom: 1rem; }
        [data-testid="stSidebar"] { background: #ffffff; }
        .stButton > button { background: linear-gradient(135deg, #059669, #0891b2); color: white; border-radius: 12px; }
    </style>
    """, unsafe_allow_html=True)


def login_page():
    st.markdown("""
    <div style="max-width: 450px; margin: 0 auto; padding: 3rem; background: white; border-radius: 24px; box-shadow: 0 20px 60px rgba(0,0,0,0.15); text-align: center;">
        <h2 style="color: #059669;">中药化合物智能鉴定平台</h2>
    </div>
    """, unsafe_allow_html=True)
    with st.form("login_form"):
        username = st.text_input("用户名", placeholder="请输入用户名")
        password = st.text_input("密码", type="password", placeholder="请输入密码")
        submitted = st.form_submit_button("登录", use_container_width=True)
        if submitted:
            if username == VALID_USERNAME and password == VALID_PASSWORD:
                st.session_state.logged_in = True
                st.session_state.username = username
                st.success("登录成功！")
                time.sleep(0.5)
                st.rerun()
            else:
                st.error("用户名或密码错误")


def create_sidebar():
    st.sidebar.markdown("""
    <div style="text-align: center; padding: 1rem;">
        <h3 style="color: #1e293b;">TCM Identifier</h3>
        <p style="color: #64748b;">中药化合物鉴定系统</p>
    </div>
    """, unsafe_allow_html=True)
    if st.session_state.get('logged_in'):
        st.sidebar.markdown(f"👤 {st.session_state.username}")
    page = st.sidebar.radio("导航", ["🏠 首页", "🔬 开始鉴定", "🧪 诊断离子筛查", "📖 使用指南", "🗃️ 数据库预览", "📊 结果分析"])
    page = page.split(' ', 1)[1] if ' ' in page else page
    st.sidebar.markdown("---")
    st.sidebar.info("💡 **性能提示**：大数据量时请启用并行处理和索引缓存（默认已启用）。")
    if st.sidebar.button("🚪 登出", use_container_width=True):
        st.session_state.logged_in = False
        st.rerun()
    return page


def create_header():
    username = st.session_state.get('username', 'Guest')
    st.markdown(f"""
    <div class="main-header">
        <h1>🌿 中药化合物智能鉴定平台</h1>
        <p>v5.15 跨行合并版 | 欢迎回来，{username}</p>
    </div>
    """, unsafe_allow_html=True)


def show_home_page():
    create_header()
    st.markdown("## 📊 系统概览")
    col1, col2, col3, col4 = st.columns(4)
    with col1: st.metric("数据库化合物数", "35,828")
    with col2: st.metric("支持药材种类", "291")
    with col3: st.metric("鉴定评级级别", "6")
    with col4: st.metric("化合物类型", "400+")
    st.markdown("---")
    st.markdown("## ✨ 核心功能")
    st.markdown("""
    - **智能化合物鉴定**：基于高分辨质谱数据精准匹配
    - **六级评级标准**：科学评估鉴定结果可靠性
    - **碎片-文献映射**：支持碎片关联文献，统计匹配文献数并加分
    - **跨行碎片合并**：同一化合物的多行碎片和文献自动合并，提高匹配覆盖度
    - **综合评分系统**：多维度加权计算（含文献匹配加分）
    """)
    if st.button("🚀 立即开始鉴定", type="primary"):
        st.session_state['page'] = '开始鉴定'
        st.rerun()


def show_analysis_page():
    create_header()
    st.markdown("## 📁 上传质谱数据")
    col1, col2 = st.columns(2)
    with col1:
        ms_positive_file = st.file_uploader("正离子模式数据 (.xlsx)", type=['xlsx'], key='ms_pos')
    with col2:
        ms_negative_file = st.file_uploader("负离子模式数据 (.xlsx)", type=['xlsx'], key='ms_neg')
    
    with st.expander("📥 下载示例数据模板"):
        st.markdown("点击下方按钮下载CSV模板，按照模板格式准备数据。")
        sample_data = pd.DataFrame({
            'Precursor M/z': [500.1234, 600.5678],
            '出峰时间t/min': [5.23, 8.45],
            'Peak_1_m/z': [150.0, 200.0],
            'Peak_1_Intensity': [1000, 800],
            'Peak_2_m/z': [250.0, 300.0],
            'Peak_2_Intensity': [500, 400]
        })
        csv_data = sample_data.to_csv(index=False).encode('utf-8')
        st.download_button("下载 CSV 模板", data=csv_data, file_name="ms_data_template.csv", mime="text/csv")
    
    st.markdown("---")
    st.markdown("## 📚 自定义数据库（可选）")
    custom_db_file = st.file_uploader("上传自定义数据库 (.xlsx，需与主数据库列一致)", type=['xlsx'], key='custom_db')
    if custom_db_file:
        st.success(f"✅ 已上传: {custom_db_file.name}")
    
    st.markdown("---")
    st.markdown("## 🧪 外部诊断离子（可选）")
    diagnostic_file = st.file_uploader("上传自定义诊断离子文件 (.xlsx，可包含权重列)", type=['xlsx'], key='diagnostic')
    if diagnostic_file:
        st.success(f"✅ 已上传: {diagnostic_file.name}")
    
    st.markdown("---")
    st.markdown("## ⚙️ 鉴定参数配置")
    
    with st.expander("📖 参数调优指南（点击展开）"):
        st.markdown("""
        - **ppm误差容限**：仪器精度高（<5 ppm）建议设为10，普通设为50。
        - **碎片匹配容差**：高分辨率质谱建议用ppm模式（5-10），低分辨率用Da（0.05）。
        - **最小绝对强度**：基质复杂样品可提高至200-500，减少噪音干扰。
        - **相对强度阈值**：仅保留强度占比≥该值的碎片，推荐1-5%。
        - **保留时间容差**：若色谱稳定可设为0.2 min，不稳定则放宽至0.5-1.0 min。
        - **启用RT得分**：仅在数据库有可靠保留时间时启用。
        """)
    
    # 参数控件
    if 'tolerance_ppm' not in st.session_state:
        st.session_state.tolerance_ppm = 50
    if 'fragment_tolerance' not in st.session_state:
        st.session_state.fragment_tolerance = 0.05
    if 'fragment_tolerance_ppm' not in st.session_state:
        st.session_state.fragment_tolerance_ppm = 20
    if 'tolerance_type' not in st.session_state:
        st.session_state.tolerance_type = 'Da'
    if 'min_intensity_abs' not in st.session_state:
        st.session_state.min_intensity_abs = 100
    if 'intensity_rel_threshold' not in st.session_state:
        st.session_state.intensity_rel_threshold = 1.0
    if 'rt_tolerance' not in st.session_state:
        st.session_state.rt_tolerance = 0.3
    if 'loss_tolerance' not in st.session_state:
        st.session_state.loss_tolerance = 0.02
    if 'use_rt_score' not in st.session_state:
        st.session_state.use_rt_score = True
    if 'use_parallel' not in st.session_state:
        st.session_state.use_parallel = True
    if 'cache_index' not in st.session_state:
        st.session_state.cache_index = True
    if 'max_candidates' not in st.session_state:
        st.session_state.max_candidates = 3
    
    colA, colB, colC, colD = st.columns(4)
    with colA: st.number_input("ppm误差容限", 10, 100, key='tolerance_ppm')
    with colB: st.number_input("最大候选数", 1, 10, key='max_candidates')
    with colC: st.number_input("最小绝对强度", 0, 1000, key='min_intensity_abs')
    with colD: st.checkbox("并行处理", key='use_parallel')
    
    colE, colF, colG = st.columns(3)
    with colE:
        st.selectbox("碎片匹配容差类型", ['Da', 'ppm'], key='tolerance_type')
        if st.session_state.tolerance_type == 'Da':
            st.number_input("碎片容差(Da)", 0.01, 1.0, step=0.01, key='fragment_tolerance')
        else:
            st.number_input("碎片容差(ppm)", 1, 100, key='fragment_tolerance_ppm')
    with colF:
        st.number_input("保留时间容差(min)", 0.1, 2.0, step=0.1, key='rt_tolerance')
    with colG:
        st.number_input("中性丢失容差(Da)", 0.01, 0.5, step=0.01, key='loss_tolerance')
    
    colH, colI, colJ = st.columns(3)
    with colH: st.slider("相对强度阈值(%)", 0.0, 100.0, key='intensity_rel_threshold')
    with colI: st.checkbox("启用RT得分", key='use_rt_score')
    with colJ: st.checkbox("索引缓存", key='cache_index')
    
    herb_filter = st.selectbox("药材筛选", ["使用全部数据库", "筛选特定药材"])
    herb_name = None
    if herb_filter == "筛选特定药材":
        herb_name = st.text_input("药材名称", placeholder="如：栀子")
    
    if ms_positive_file or ms_negative_file:
        if st.button("🚀 开始鉴定", type="primary"):
            temp_files = []
            try:
                temp_dir = tempfile.gettempdir()
                pos_path = None
                neg_path = None
                diag_path = None
                custom_db_path = None
                if ms_positive_file:
                    pos_path = os.path.join(temp_dir, ms_positive_file.name)
                    with open(pos_path, 'wb') as f:
                        f.write(ms_positive_file.getbuffer())
                    temp_files.append(pos_path)
                if ms_negative_file:
                    neg_path = os.path.join(temp_dir, ms_negative_file.name)
                    with open(neg_path, 'wb') as f:
                        f.write(ms_negative_file.getbuffer())
                    temp_files.append(neg_path)
                if diagnostic_file:
                    diag_path = os.path.join(temp_dir, diagnostic_file.name)
                    with open(diag_path, 'wb') as f:
                        f.write(diagnostic_file.getbuffer())
                    temp_files.append(diag_path)
                if custom_db_file:
                    custom_db_path = os.path.join(temp_dir, custom_db_file.name)
                    with open(custom_db_path, 'wb') as f:
                        f.write(custom_db_file.getbuffer())
                    temp_files.append(custom_db_path)
                
                db_path = find_database_path()
                if not db_path:
                    st.error("未找到主数据库文件！请将 TCM-SM-MS DB.xlsx 放在项目目录下。")
                    return
                
                config = {
                    'min_intensity': st.session_state.min_intensity_abs,
                    'fragment_tolerance': st.session_state.fragment_tolerance,
                    'fragment_tolerance_ppm': st.session_state.fragment_tolerance_ppm,
                    'tolerance_ppm': st.session_state.tolerance_ppm,
                    'max_candidates': st.session_state.max_candidates,
                }
                
                with st.spinner("鉴定中，请稍候..."):
                    identifier = UltimateGardeniaIdentifier(
                        database_path=db_path,
                        ms_positive_path=pos_path,
                        ms_negative_path=neg_path,
                        herb_name=herb_name,
                        config=config,
                        use_parallel=st.session_state.use_parallel,
                        rt_tolerance=st.session_state.rt_tolerance,
                        loss_tolerance=st.session_state.loss_tolerance,
                        external_diagnostic_file=diag_path,
                        rt_fusion_tolerance=st.session_state.rt_tolerance,
                        intensity_relative_threshold=st.session_state.intensity_rel_threshold,
                        tolerance_type=st.session_state.tolerance_type,
                        use_rt_score=st.session_state.use_rt_score,
                        custom_db_path=custom_db_path,
                        cache_index=st.session_state.cache_index
                    )
                    report = identifier.generate_report('样品')
                    st.session_state['analysis_results'] = report
                    st.success(f"鉴定完成！共鉴定出 {len(report)} 个化合物。")
                    if st.button("📊 查看结果"):
                        st.session_state['page'] = '结果分析'
                        st.rerun()
            except Exception as e:
                st.error(f"鉴定出错: {str(e)}")
                st.exception(e)
            finally:
                for f in temp_files:
                    try: os.remove(f)
                    except: pass
    else:
        st.info("请至少上传一个质谱数据文件")


def show_diagnostic_ion_page():
    create_header()
    st.markdown("## 🔬 诊断离子筛查")
    diagnostic_df = load_diagnostic_ions_cached()
    if diagnostic_df.empty:
        st.warning("未找到诊断离子数据库文件")
        return
    st.success(f"已加载 {len(diagnostic_df)} 条诊断离子记录")
    
    with st.expander("📥 下载诊断离子模板"):
        template = pd.DataFrame({
            '化合物类型': ['黄酮类', '生物碱类'],
            '诊断碎片离子m/z': [151.003, 105.070],
            '权重': [1, 2],
            '描述': ['特征碎片', '特征碎片'],
            '离子模式': ['正离子', '正离子']
        })
        buffer = BytesIO()
        with pd.ExcelWriter(buffer, engine='openpyxl') as writer:
            template.to_excel(writer, index=False, sheet_name='诊断离子模板')
        st.download_button("下载模板.xlsx", data=buffer.getvalue(), file_name="diagnostic_template.xlsx")
    
    mz_input = st.text_area("输入m/z值（每行一个，或逗号分隔）", height=150)
    tolerance_ppm = st.number_input("ppm容差", 1, 100, 10)
    ion_mode = st.selectbox("离子模式", ["全部", "正离子", "负离子"])
    
    if mz_input:
        mz_list = []
        for token in re.split(r'[,\s\n]+', mz_input):
            token = token.strip()
            if token:
                try:
                    mz_list.append(float(token))
                except:
                    pass
        if mz_list:
            results = match_diagnostic_ions(mz_list, diagnostic_df, tolerance_ppm, ion_mode)
            if not results.empty:
                st.dataframe(results)
                csv = results.to_csv(index=False).encode('utf-8')
                st.download_button("导出CSV", csv, "diagnostic_results.csv")
            else:
                st.info("未找到匹配的诊断离子")


def show_guide_page():
    create_header()
    st.markdown("## 📖 使用指南")
    st.markdown("""
    ### 快速入门
    1. 准备质谱数据（Excel或CSV），包含母离子m/z、保留时间、碎片m/z及强度。
    2. 上传文件，设置参数（新手可使用快速模式）。
    3. 点击开始鉴定，等待结果。
    4. 查看结果报告，根据评级和得分判断可靠性。

    ### 评分规则（v5.15 新增跨行合并）
    - 基础分：确证级85，高置信级65，推定级45，提示级25，参考级0
    - 加分项：
      - 匹配碎片：每个+2（上限20）
      - 诊断离子：按权重+5/单位（上限15）
      - 中性丢失：每个+2（上限10）
      - RT匹配：偏差<0.2 +5，<0.5 +2
      - **文献匹配：每篇+3（上限12）**
    - 扣分项：ppm>10开始扣分
    - 保底：1级不低于80，2级不低于60

    ### 碎片-文献映射格式
    在数据库的“碎片离子（正/负）”列中，可以为每个碎片指定文献索引：
    151.003:0,1; 137.024:1,2
    表示碎片151.003出现在第0、1篇文献，碎片137.024出现在第1、2篇文献。
文献索引对应“文献来源”列中分号分隔的顺序（从0开始）。
若不指定索引（如`151.003`），则默认属于所有文献。

### 常见问题
- **鉴定结果为空**：检查ppm容差是否过小，或药材筛选是否正确。
- **缓存加载失败**：删除 `index_cache.pkl` 文件后重试。
- **文献匹配数为0**：检查数据库是否配置了带索引的碎片字段。
- **同一化合物多行信息未合并**：v5.15已自动合并，无需手动处理。
""")


def show_database_page():
    create_header()
    st.markdown("## 🗃️ 数据库预览")
    db_path = find_database_path()
    if not db_path:
        st.warning("未找到数据库文件")
        return
    df = load_database_cached()
    st.dataframe(df.head(10))
    st.info(f"共 {len(df)} 条记录")


def show_results_page():
    create_header()
    if 'analysis_results' not in st.session_state:
        st.warning("暂无鉴定结果，请先进行鉴定")
        if st.button("前往鉴定"):
            st.session_state['page'] = '开始鉴定'
            st.rerun()
        return
    report = st.session_state['analysis_results']
    if report.empty:
        st.warning("鉴定结果为空")
        return

    st.markdown("## 📊 鉴定结果")
    col1, col2, col3, col4 = st.columns(4)
    with col1: st.metric("鉴定总数", len(report))
    with col2: st.metric("确证级", (report['评级名称']=='确证级').sum())
    with col3: st.metric("90分以上", (report['综合得分']>=90).sum())
    with col4: st.metric("平均得分", f"{report['综合得分'].mean():.1f}")

    with st.expander("📖 结果解读示例（基于当前最高分化合物）"):
        top = report.iloc[0]
        st.markdown(f"""
        **化合物**: {top['化合物中文名']} (m/z {top['m/z实际值']:.4f})
        - ppm误差 = {top['ppm']:.2f} → {"高精度匹配" if top['ppm']<10 else "可接受"}
        - 匹配碎片数 = {top['匹配碎片数']} → {"二级质谱验证充分" if top['匹配碎片数']>=3 else "碎片信息较少"}
        - **匹配文献数 = {top['匹配文献数']}** → {"文献支持充分" if top['匹配文献数']>=2 else "文献支持较少"}
        - 诊断离子个数 = {top['诊断性离子个数']} → {"符合类别特征" if top['诊断性离子个数']>0 else "无类别特异性"}
        - 综合得分 = {top['综合得分']:.1f} → {top['评级名称']}
        - 报告建议: {top['报告建议']}
        """)

    st.dataframe(report, use_container_width=True)

    csv = report.to_csv(index=False).encode('utf-8')
    st.download_button("导出CSV", csv, "identification_report.csv")


def main():
    load_optimized_css()
    if 'logged_in' not in st.session_state:
        st.session_state.logged_in = False
    if not st.session_state.logged_in:
        login_page()
        return
    if 'page' not in st.session_state:
        st.session_state['page'] = '首页'
    page = create_sidebar()
    st.session_state['page'] = page
    if page == "首页":
        show_home_page()
    elif page == "开始鉴定":
        show_analysis_page()
    elif page == "诊断离子筛查":
        show_diagnostic_ion_page()
    elif page == "使用指南":
        show_guide_page()
    elif page == "数据库预览":
        show_database_page()
    elif page == "结果分析":
        show_results_page()


if __name__ == "__main__":
    main()
    
