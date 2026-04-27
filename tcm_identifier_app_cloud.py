# -*- coding: utf-8 -*-
"""
中药化合物智能鉴定平台 - 云端部署增强版 v5.13（文献-碎片详细映射）
==========================================================
更新说明：
- 一级匹配（准分子离子）：使用ppm容差控制
- 二级匹配（碎片离子）：使用绝对Da容差控制（可调整）
- 文献输出：支持详细映射格式 139.05(文献1,文献2); 134.03(文献1,文献3)
- 支持碎片对应多篇文献的映射关系
- 修复 min_intensity_abs 未定义错误
"""

import streamlit as st
import pandas as pd
import numpy as np
import os
import time
import pickle
import hashlib
import re
from datetime import datetime
from io import BytesIO
import tempfile
from bisect import bisect_left, bisect_right
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')


# ============================================================================
# 登录验证常量
# ============================================================================
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
# 碎片离子解析函数（支持多值映射，支持文献来源绑定）
# ============================================================================

def parse_fragments_with_sources(fragment_string, source_string):
    """
    解析碎片离子字符串，并绑定文献来源
    返回: [(mz, [source1, source2, ...]), ...] 每个碎片对应的文献列表
    """
    if pd.isna(fragment_string) or not fragment_string or str(fragment_string).strip() == '':
        return []
    
    s = str(fragment_string).strip()
    fragments = []
    
    # 尝试多种分隔符
    if '、' in s:
        parts = s.split('、')
    elif ';' in s:
        parts = s.split(';')
    elif ',' in s:
        parts = s.split(',')
    elif ' ' in s:
        parts = s.split()
    else:
        parts = [s]
    
    # 解析文献来源列表
    sources = []
    if not pd.isna(source_string) and source_string:
        source_str = str(source_string).strip()
        if ';' in source_str:
            sources = [s.strip() for s in source_str.split(';') if s.strip()]
        elif ',' in source_str:
            sources = [s.strip() for s in source_str.split(',') if s.strip()]
        else:
            sources = [source_str]
    
    for part in parts:
        part = part.strip()
        if part:
            try:
                match = re.search(r'[-+]?\d+\.?\d*', part)
                if match:
                    frag_value = float(match.group())
                    fragments.append((frag_value, sources.copy()))
                else:
                    fragments.append((float(part), sources.copy()))
            except ValueError:
                continue
    
    return fragments


def parse_fragments_simple(fragment_string):
    """简单解析碎片离子，仅返回mz列表"""
    if pd.isna(fragment_string) or not fragment_string or str(fragment_string).strip() == '':
        return []
    
    s = str(fragment_string).strip()
    fragments = []
    
    if '、' in s:
        parts = s.split('、')
    elif ';' in s:
        parts = s.split(';')
    elif ',' in s:
        parts = s.split(',')
    elif ' ' in s:
        parts = s.split()
    else:
        parts = [s]
    
    for part in parts:
        part = part.strip()
        if part:
            try:
                match = re.search(r'[-+]?\d+\.?\d*', part)
                if match:
                    fragments.append(float(match.group()))
                else:
                    fragments.append(float(part))
            except ValueError:
                continue
    
    return fragments


def format_fragment_sources(fragment_sources_map):
    """
    格式化碎片-文献映射关系
    输入: {139.05: ['文献1', '文献2'], 134.03: ['文献1', '文献3']}
    输出: "139.05(文献1,文献2); 134.03(文献1,文献3)"
    """
    if not fragment_sources_map:
        return ''
    
    items = []
    for mz in sorted(fragment_sources_map.keys()):
        sources = fragment_sources_map[mz]
        if sources:
            unique_sources = []
            for s in sources:
                if s not in unique_sources:
                    unique_sources.append(s)
            sources_str = ','.join(unique_sources)
            items.append(f"{mz:.4f}({sources_str})")
        else:
            items.append(f"{mz:.4f}")
    
    return '; '.join(items)


# ============================================================================
# 数据库加载函数
# ============================================================================

@st.cache_data
def load_database_cached(db_filename="TCM-SM-MS DB.xlsx"):
    """加载主数据库"""
    db_filenames = [
        "TCM-SM-MS DB.xlsx",
        "TCM-SM-MS DB（中药小分子化学成分高分辨质谱数据库）.xlsx",
        "data/TCM-SM-MS DB.xlsx",
        "data/TCM-SM-MS DB（中药小分子化学成分高分辨质谱数据库）.xlsx",
        "user_input_files/TCM-SM-MS DB.xlsx",
        "user_input_files/TCM-SM-MS DB（中药小分子化学成分高分辨质谱数据库）.xlsx"
    ]
    if db_filename:
        db_filenames = [db_filename] + [f for f in db_filenames if f != db_filename]
    for path in db_filenames:
        if os.path.exists(path):
            try:
                df = pd.read_excel(path)
                print(f"主数据库加载成功: {path}, 共 {len(df)} 条记录")
                df = standardize_db_columns(df)
                return df
            except Exception as e:
                print(f"主数据库加载失败 {path}: {e}")
                continue
    print("警告: 未找到主数据库文件")
    return pd.DataFrame()


@st.cache_data
def load_english_database_cached(db_filename="数据库（英文）.xlsx"):
    """加载英文数据库"""
    english_paths = [
        "数据库（英文）.xlsx",
        "data/数据库（英文）.xlsx",
        "user_input_files/数据库（英文）.xlsx"
    ]
    if db_filename:
        english_paths = [db_filename] + [f for f in english_paths if f != db_filename]
    for path in english_paths:
        if os.path.exists(path):
            try:
                df = pd.read_excel(path)
                print(f"英文数据库加载成功: {path}, 共 {len(df)} 条记录")
                df = standardize_db_columns(df)
                return df
            except Exception as e:
                print(f"英文数据库加载失败 {path}: {e}")
                continue
    print("提示: 未找到英文数据库文件，将仅使用主数据库")
    return pd.DataFrame()


def standardize_db_columns(df):
    """标准化数据库列名"""
    column_mapping = {
        '药材名': '药材名称',
        '文献': '文献来源',
        '保留时间(min)': '保留时间',
        '中丢失': '中性丢失',
        'CAS': 'CAS号'
    }
    for old_name, new_name in column_mapping.items():
        if old_name in df.columns and new_name not in df.columns:
            df.rename(columns={old_name: new_name}, inplace=True)
    
    if '名称（中文）' not in df.columns:
        df['名称（中文）'] = ''
    if '名称（英文）' not in df.columns:
        df['名称（英文）'] = ''
    if '碎片离子（正）' not in df.columns:
        df['碎片离子（正）'] = ''
    if '碎片离子（负）' not in df.columns:
        df['碎片离子（负）'] = ''
    if '准分子离子（正）' not in df.columns:
        df['准分子离子（正）'] = ''
    if '准分子离子（负）' not in df.columns:
        df['准分子离子（负）'] = ''
    if '加合物（正）' not in df.columns:
        df['加合物（正）'] = ''
    if '加合物（负）' not in df.columns:
        df['加合物（负）'] = ''
    if '化合物类型' not in df.columns:
        df['化合物类型'] = ''
    if '文献来源' not in df.columns:
        df['文献来源'] = ''
    
    return df


def find_database_path():
    """查找主数据库路径"""
    db_paths = [
        "TCM-SM-MS DB.xlsx",
        "TCM-SM-MS DB（中药小分子化学成分高分辨质谱数据库）.xlsx",
        "data/TCM-SM-MS DB.xlsx",
        "data/TCM-SM-MS DB（中药小分子化学成分高分辨质谱数据库）.xlsx",
        "user_input_files/TCM-SM-MS DB.xlsx",
        "user_input_files/TCM-SM-MS DB（中药小分子化学成分高分辨质谱数据库）.xlsx"
    ]
    for path in db_paths:
        if os.path.exists(path):
            return path
    return None


def find_english_database_path():
    """查找英文数据库路径"""
    english_paths = [
        "数据库（英文）.xlsx",
        "data/数据库（英文）.xlsx",
        "user_input_files/数据库（英文）.xlsx"
    ]
    for path in english_paths:
        if os.path.exists(path):
            return path
    return None


@st.cache_data
def load_diagnostic_ions_cached():
    """加载诊断离子库"""
    diagnostic_ion_paths = [
        "诊断离子.xlsx",
        "data/诊断离子.xlsx",
        "user_input_files/诊断离子.xlsx"
    ]
    for path in diagnostic_ion_paths:
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
    """查找诊断离子库路径"""
    diagnostic_ion_paths = [
        "诊断离子.xlsx",
        "data/诊断离子.xlsx",
        "user_input_files/诊断离子.xlsx"
    ]
    for path in diagnostic_ion_paths:
        if os.path.exists(path):
            return path
    return None


def match_diagnostic_ions(user_mz_values, diagnostic_df, tolerance_ppm=10, ion_mode=None):
    """匹配诊断离子"""
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
# 鉴定程序核心代码
# ============================================================================

class UltimateGardeniaIdentifier:
    """
    中药化合物鉴定终极版程序 v5.13（文献-碎片详细映射）
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
                 cache_index=True,
                 english_db_path=None,
                 fragment_tolerance_no_ppm=0.1):
        """初始化鉴定程序"""
        self.config = {
            'gradient_time': 30.0,
            'cid_min': 0.5,
            'cid_max': 1.0,
            'fragment_tolerance': 0.05,
            'fragment_tolerance_ppm': 20,
            'fragment_tolerance_no_ppm': fragment_tolerance_no_ppm,
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
        print("中药化合物鉴定程序 v5.13（文献-碎片详细映射）")
        print("="*80)
        
        self.full_database = pd.DataFrame()
        self.database = pd.DataFrame()
        
        print("\n【1/9】正在加载主数据库...")
        if database_path and os.path.exists(database_path):
            self.full_database = self._load_data(database_path)
            print(f"  主数据库加载成功: {len(self.full_database)} 条记录")
        else:
            print("  警告: 未找到主数据库文件")
        
        print("【2/9】正在加载英文数据库...")
        if english_db_path is None:
            english_db_path = find_english_database_path()
        if english_db_path and os.path.exists(english_db_path):
            english_db = self._load_data(english_db_path)
            if not english_db.empty:
                print(f"  英文数据库加载成功: {len(english_db)} 条记录")
                english_db['_db_source'] = 'english'
                self.full_database = pd.concat([self.full_database, english_db], ignore_index=True)
            else:
                print("  英文数据库为空或加载失败")
        else:
            print("  未找到英文数据库文件，将仅使用主数据库")
        
        if custom_db_path and os.path.exists(custom_db_path):
            custom_db = self._load_data(custom_db_path)
            if not custom_db.empty:
                print(f"  加载自定义数据库: {len(custom_db)} 条记录")
                custom_db['_db_source'] = 'custom'
                self.full_database = pd.concat([self.full_database, custom_db], ignore_index=True)

        if herb_name:
            print(f"【3/9】正在筛选 {herb_name} 相关数据...")
            self.database = self._filter_by_herb(herb_name)
        else:
            self.database = self.full_database.copy()
            print("【3/9】使用全部数据库进行化合物鉴定")
        print(f"  数据库记录数: {len(self.database)}")

        print("【4/9】正在加载质谱数据...")
        self.ms_positive = self._load_data(ms_positive_path) if ms_positive_path else pd.DataFrame()
        self.ms_negative = self._load_data(ms_negative_path) if ms_negative_path else pd.DataFrame()
        print(f"  正离子数据: {len(self.ms_positive)} 条记录")
        print(f"  负离子数据: {len(self.ms_negative)} 条记录")

        print("【5/9】正在构建索引...")
        self._build_or_load_index()

        print("【6/9】正在加载诊断离子库...")
        if external_diagnostic_file is None:
            default_diag_path = find_diagnostic_ion_path()
            if default_diag_path:
                print(f"  自动找到默认诊断离子文件: {default_diag_path}")
                external_diagnostic_file = default_diag_path
        self._build_diagnostic_ion_library(external_diagnostic_file)

        print("【7/9】正在加载辅助数据...")
        self._load_auxiliary_data()

        self.stats = {
            'total_precursors': 0,
            'identified_compounds': 0,
            'scored_candidates': 0
        }

        self._print_initialization_info()
        print("【8/9】初始化完成，准备鉴定")

    def _load_data(self, filepath):
        """加载数据文件"""
        if filepath and os.path.exists(filepath):
            try:
                if filepath.endswith('.xlsx'):
                    df = pd.read_excel(filepath)
                elif filepath.endswith('.csv'):
                    return pd.read_csv(filepath)

                column_mapping = {
                    '药材名': '药材名称',
                    '文献': '文献来源',
                    '保留时间(min)': '保留时间',
                    '中丢失': '中性丢失',
                    'CAS': 'CAS号'
                }
                for old_name, new_name in column_mapping.items():
                    if old_name in df.columns and new_name not in df.columns:
                        df.rename(columns={old_name: new_name}, inplace=True)

                if '名称（中文）' not in df.columns:
                    df['名称（中文）'] = ''
                if '名称（英文）' not in df.columns:
                    df['名称（英文）'] = ''
                if '碎片离子（正）' not in df.columns:
                    df['碎片离子（正）'] = ''
                if '碎片离子（负）' not in df.columns:
                    df['碎片离子（负）'] = ''
                if '准分子离子（正）' not in df.columns:
                    df['准分子离子（正）'] = ''
                if '准分子离子（负）' not in df.columns:
                    df['准分子离子（负）'] = ''
                if '加合物（正）' not in df.columns:
                    df['加合物（正）'] = ''
                if '加合物（负）' not in df.columns:
                    df['加合物（负）'] = ''
                if '化合物类型' not in df.columns:
                    df['化合物类型'] = ''
                if '文献来源' not in df.columns:
                    df['文献来源'] = ''

                return df
            except Exception as e:
                print(f"警告: 无法加载文件 {filepath}: {e}")
                return pd.DataFrame()
        print(f"警告: 文件不存在或路径无效: {filepath}")
        return pd.DataFrame()

    def _filter_by_herb(self, herb_name):
        """按药材筛选数据库"""
        if herb_name is None:
            return self.full_database.copy()

        herb_col = '药材名称' if '药材名称' in self.full_database.columns else '药材名'
        mask = self.full_database[herb_col].str.contains(herb_name, na=False, case=False)
        filtered_db = self.full_database[mask].copy()
        if len(filtered_db) == 0:
            print(f"警告: 未找到 '{herb_name}' 相关数据，将使用全部数据")
            return self.full_database.copy()
        return filtered_db

    def _build_or_load_index(self):
        """构建或加载索引缓存"""
        cache_file = "index_cache_v513.pkl"
        db_hash = hashlib.md5(pd.util.hash_pandas_object(self.database).values).hexdigest()

        if self.cache_index and os.path.exists(cache_file):
            try:
                with open(cache_file, 'rb') as f:
                    cached = pickle.load(f)
                if cached.get('db_hash') == db_hash:
                    print("  从缓存加载索引...")
                    self.sorted_idx_pos = cached['sorted_idx_pos']
                    self.sorted_idx_neg = cached['sorted_idx_neg']
                    self.mz_values_pos = cached['mz_values_pos']
                    self.mz_values_neg = cached['mz_values_neg']
                    self.compound_info = cached['compound_info']
                    self.rt_values = cached.get('rt_values', {})
                    self.neutral_losses = cached.get('neutral_losses', {})
                    print("  缓存加载成功")
                    return
            except Exception as e:
                print(f"  缓存加载失败，重新构建: {e}")

        self._build_optimized_index()

        if self.cache_index:
            try:
                cached = {
                    'db_hash': db_hash,
                    'sorted_idx_pos': self.sorted_idx_pos,
                    'sorted_idx_neg': self.sorted_idx_neg,
                    'mz_values_pos': self.mz_values_pos,
                    'mz_values_neg': self.mz_values_neg,
                    'compound_info': self.compound_info,
                    'rt_values': self.rt_values,
                    'neutral_losses': self.neutral_losses
                }
                with open(cache_file, 'wb') as f:
                    pickle.dump(cached, f)
                print("  索引已缓存")
            except Exception as e:
                print(f"  缓存保存失败: {e}")

    def _build_optimized_index(self):
        """构建优化索引"""
        self.sorted_idx_pos = []
        self.sorted_idx_neg = []
        self.mz_values_pos = []
        self.mz_values_neg = []
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
            source = str(row.get('文献来源') or row.get('文献') or '')
            cas = str(row.get('CAS号') or row.get('CAS', ''))
            db_source = row.get('_db_source', 'main')

            self.compound_info[idx] = {
                'name_cn': cn_name,
                'name_en': en_name,
                'formula': formula,
                'herb': herb,
                'compound_type': compound_type,
                'source': source,
                'cas': cas if cas and cas != 'nan' else '',
                'adduct_pos': str(row.get('加合物（正）', '')),
                'adduct_neg': str(row.get('加合物（负）', '')),
                'db_source': db_source
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

            if pd.notna(mz_pos) and str(mz_pos).strip() != '':
                try:
                    mz_val = float(mz_pos)
                    if mz_val > 0:
                        fragments_with_sources = parse_fragments_with_sources(
                            row.get('碎片离子（正）', ''),
                            source
                        )
                        self.sorted_idx_pos.append((mz_val, idx, fragments_with_sources))
                except (ValueError, TypeError):
                    pass

            if pd.notna(mz_neg) and str(mz_neg).strip() != '':
                try:
                    mz_val = float(mz_neg)
                    if mz_val > 0:
                        fragments_with_sources = parse_fragments_with_sources(
                            row.get('碎片离子（负）', ''),
                            source
                        )
                        self.sorted_idx_neg.append((mz_val, idx, fragments_with_sources))
                except (ValueError, TypeError):
                    pass

        self.sorted_idx_pos.sort(key=lambda x: x[0])
        self.sorted_idx_neg.sort(key=lambda x: x[0])

        self.mz_values_pos = np.array([x[0] for x in self.sorted_idx_pos])
        self.mz_values_neg = np.array([x[0] for x in self.sorted_idx_neg])

    def _match_fragments_with_sources(self, observed_fragments, db_fragments_with_sources, 
                                       tolerance_value, precursor_mz=None):
        """匹配碎片并返回碎片-文献映射"""
        if len(db_fragments_with_sources) == 0:
            return [], {}
        
        observed = np.asarray(observed_fragments)
        matched_fragments = []
        fragment_sources_map = defaultdict(list)
        
        for ref_mz, sources in db_fragments_with_sources:
            if pd.notna(ref_mz) and float(ref_mz) > 0:
                for obs_mz in observed:
                    if precursor_mz is not None and abs(obs_mz - precursor_mz) <= tolerance_value:
                        continue
                    
                    if abs(obs_mz - ref_mz) <= tolerance_value:
                        if obs_mz not in matched_fragments:
                            matched_fragments.append(obs_mz)
                        for s in sources:
                            if s and s not in fragment_sources_map[obs_mz]:
                                fragment_sources_map[obs_mz].append(s)
                        break
        
        return matched_fragments, fragment_sources_map

    def _build_diagnostic_ion_library(self, external_file=None):
        """构建诊断性离子库"""
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
                            mz = float(mz)
                            weight = row.get('权重', 1)
                            try:
                                weight = float(weight)
                            except:
                                weight = 1
                            ions_dict[mz] = ions_dict.get(mz, 0) + weight
                    description = group['描述'].iloc[0] if '描述' in group.columns else f'来自外部文件，{len(ions_dict)}个离子'
                    self.diagnostic_ions[category] = {
                        'ions': list(ions_dict.keys()),
                        'weights': list(ions_dict.values()),
                        'description': description
                    }
                print(f"  成功加载外部诊断离子库：{len(self.diagnostic_ions)} 类，{len(df)} 条记录")
            except Exception as e:
                print(f"  加载外部诊断离子文件失败：{e}，将使用内置库")
                self._build_default_diagnostic_ions()
        else:
            print("  未找到外部诊断离子文件，使用内置诊断离子库")
            self._build_default_diagnostic_ions()

    def _build_default_diagnostic_ions(self):
        """构建默认诊断离子库"""
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

    def _parse_losses(self, loss_string):
        """解析中性丢失字符串"""
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

    def _find_neutral_losses(self, fragments, precursor_mz):
        """查找中性丢失"""
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
        """分类化合物类型"""
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
        """二分查找匹配范围"""
        if len(mz_array) == 0:
            return range(0, 0)
        mz = float(mz)
        tolerance = mz * tolerance_ppm / 1e6
        mz_min = mz - tolerance
        mz_max = mz + tolerance
        left = bisect_left(mz_array, mz_min)
        right = bisect_right(mz_array, mz_max)
        return range(left, right)

    def _find_diagnostic_ions_fast(self, matched_fragments, category, precursor_mz=None):
        """快速查找诊断离子"""
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
        """确定置信等级"""
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
                              diag_weights=None, neutral_loss_match_count=0, rt_deviation=None):
        """计算基础得分"""
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

        total = base + ppm_adj + frag_adj + diag_score + loss_adj + rt_adj
        total = max(0, total)

        if rating == 1 and total < 80:
            total = 80
        if rating == 2 and total < 60:
            total = 60

        total = min(total, 100)
        return round(total, 2)

    def extract_precursor_ions(self, ms_data, ionization_mode):
        """提取母离子"""
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
        """鉴定化合物"""
        mode = precursor['mode']
        if mode == '正离子':
            sorted_idx = self.sorted_idx_pos
            mz_array = self.mz_values_pos
        else:
            sorted_idx = self.sorted_idx_neg
            mz_array = self.mz_values_neg

        # 一级匹配：使用ppm容差
        candidate_range = self._binary_search_range(
            mz_array, 
            precursor['precursor_mz'], 
            self.config['tolerance_ppm']
        )

        scored_candidates = []
        for idx in candidate_range:
            mz_val, db_idx, db_fragments_with_sources = sorted_idx[idx]
            if db_idx not in self.compound_info:
                continue

            info = self.compound_info[db_idx]
            category = self._classify_compound(info['name_cn'] + ' ' + info['name_en'], info['compound_type'])

            # 二级匹配：使用绝对Da容差
            tolerance_value = self.config['fragment_tolerance_no_ppm']
            
            matched_fragments, fragment_sources_map = self._match_fragments_with_sources(
                precursor['fragments'],
                db_fragments_with_sources,
                tolerance_value,
                precursor['precursor_mz']
            )
            
            matched_fragments_str = format_fragment_sources(fragment_sources_map)
            matched_fragments_list = list(fragment_sources_map.keys())

            diagnostic_ions, diag_weights = self._find_diagnostic_ions_fast(
                matched_fragments_list,
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

            all_sources = set()
            for sources in fragment_sources_map.values():
                for s in sources:
                    if s:
                        all_sources.add(s)
            source_display = '; '.join(sorted(all_sources)) if all_sources else info['source']

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
                'matched_fragments': matched_fragments_list,
                'matched_fragments_str': matched_fragments_str,
                'diagnostic_ions': diagnostic_ions,
                'diag_weights': diag_weights,
                'mode': mode,
                'source': source_display,
                'category': category,
                'db_index': db_idx,
                'neutral_loss_matches': neutral_loss_matches,
                'rt_deviation': rt_deviation,
                'db_source': info.get('db_source', 'main'),
                'fragment_sources_map': fragment_sources_map
            }

            candidate['temp_score'] = -ppm
            scored_candidates.append(candidate)

        scored_candidates.sort(key=lambda x: (-len(x['matched_fragments']), -x['temp_score']))

        if return_top == 1:
            return scored_candidates[0] if scored_candidates else None
        else:
            return scored_candidates[:return_top]

    def generate_report(self, herb_name=None):
        """生成鉴定报告"""
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

            source_display = best_candidate['source']
            matched_frags = best_candidate['matched_fragments']
            diag_ions = best_candidate['diagnostic_ions']
            diag_weights = best_candidate['diag_weights']
            
            matched_fragments_display = best_candidate.get('matched_fragments_str', '')
            if not matched_fragments_display:
                matched_fragments_display = '; '.join([f'{f:.4f}' for f in matched_frags]) if matched_frags else ''

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
                best_candidate['rt_deviation']
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
                '主要碎片离子(文献映射)': matched_fragments_display,
                '匹配碎片数': len(matched_frags),
                '诊断性离子个数': len(diag_ions),
                '诊断性离子': '; '.join([f'{f:.4f}' for f in diag_ions]) if diag_ions else '',
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
                    all_frag_sets = []
                    for rec in comp_recs:
                        if rec['加和离子']:
                            adducts.add(rec['加和离子'])
                        if rec['离子化方式']:
                            modes.add(rec['离子化方式'])
                        all_frag_sets.append(rec['_fragments_set'])
                    
                    best_rec['加和离子'] = '; '.join(sorted(adducts))
                    best_rec['离子化方式'] = '/'.join(sorted(modes))

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
        print("\n" + "="*100)

    def _print_initialization_info(self):
        print("\n" + "="*80)
        print("程序初始化完成（v5.13，文献-碎片详细映射）")
        print("="*80)
        print(f"  - 数据库记录数: {len(self.database)} 条")
        print(f"  - 正离子索引: {len(self.mz_values_pos)} 条")
        print(f"  - 负离子索引: {len(self.mz_values_neg)} 条")
        print(f"  - 一级匹配ppm容差: {self.config['tolerance_ppm']} ppm")
        print(f"  - 二级匹配绝对Da容差: {self.config['fragment_tolerance_no_ppm']} Da")
        print("="*80)


# ============================================================================
# Streamlit 网页应用部分
# ============================================================================

st.set_page_config(
    page_title="中药化合物智能鉴定平台 v5.13",
    page_icon="🌿",
    layout="wide",
    initial_sidebar_state="expanded"
)


def load_optimized_css():
    st.markdown("""
    <style>
        .stApp { background: linear-gradient(135deg, #f8fafc 0%, #e2e8f0 50%, #f1f5f9 100%); min-height: 100vh; }
        .main-header { background: linear-gradient(135deg, #059669 0%, #0891b2 50%, #7c3aed 100%); padding: 2.5rem 2rem; border-radius: 20px; margin-bottom: 2rem; color: white; box-shadow: 0 20px 40px rgba(5,150,105,0.3); }
        .main-header h1 { color: white !important; font-size: 2.5rem; font-weight: 700; margin: 0; }
        .main-header p { color: rgba(255,255,255,0.9); font-size: 1.1rem; margin-top: 0.5rem; }
        .stat-card { background: rgba(255,255,255,0.95); border-radius: 16px; padding: 1.5rem; text-align: center; margin: 0.5rem; box-shadow: 0 4px 20px rgba(0,0,0,0.08); transition: all 0.3s ease; }
        .stat-card:hover { transform: translateY(-5px); box-shadow: 0 8px 30px rgba(5,150,105,0.2); }
        .stat-number { font-size: 2.5rem; font-weight: 700; background: linear-gradient(135deg, #059669, #0891b2); -webkit-background-clip: text; -webkit-text-fill-color: transparent; }
        .stat-label { font-size: 0.95rem; color: #64748b; margin-top: 0.5rem; font-weight: 500; }
        .feature-card { background: white; border-radius: 16px; padding: 2rem; box-shadow: 0 4px 20px rgba(0,0,0,0.08); margin-bottom: 1.5rem; transition: all 0.3s ease; }
        .feature-card:hover { transform: translateX(5px); box-shadow: 0 8px 30px rgba(5,150,105,0.15); }
        .feature-card h4 { color: #1e293b; font-size: 1.25rem; font-weight: 600; margin-bottom: 0.75rem; display: flex; align-items: center; gap: 0.5rem; }
        .feature-card p { color: #64748b; font-size: 0.95rem; line-height: 1.6; margin: 0; }
        [data-testid="stSidebar"] { background: linear-gradient(180deg, #ffffff 0%, #f8fafc 100%); border-right: 1px solid #e2e8f0; }
        .stButton > button { background: linear-gradient(135deg, #059669 0%, #0891b2 100%); color: white; border: none; border-radius: 12px; padding: 0.75rem 1.5rem; font-weight: 600; font-size: 1rem; transition: all 0.3s ease; }
        .stButton > button:hover { transform: translateY(-3px); box-shadow: 0 8px 25px rgba(5,150,105,0.4); }
        .stProgress > div > div { background: linear-gradient(90deg, #059669, #0891b2, #7c3aed); border-radius: 10px; }
        [data-testid="stFileUploader"] { border: 2px dashed #cbd5e1; border-radius: 16px; padding: 2rem; background: rgba(5,150,105,0.02); transition: all 0.3s ease; }
        [data-testid="stFileUploader"]:hover { border-color: #059669; background: rgba(5,150,105,0.05); }
        hr { border: none; height: 2px; background: linear-gradient(90deg, transparent, #e2e8f0, transparent); margin: 2rem 0; }
        h1, h2, h3 { color: #1e293b !important; font-weight: 600 !important; }
        h2 { font-size: 1.75rem; margin-bottom: 1rem; }
        @media (max-width: 768px) { .main-header { padding: 1.5rem 1rem; } .main-header h1 { font-size: 1.75rem; } .stat-card { padding: 1rem; } .stat-number { font-size: 2rem; } }
        .empty-state { text-align: center; padding: 3rem; color: #94a3b8; }
    </style>
    """, unsafe_allow_html=True)


def login_page():
    st.markdown("""
    <div style="max-width: 450px; margin: 0 auto; padding: 3rem; background: white; border-radius: 24px; box-shadow: 0 20px 60px rgba(0,0,0,0.15); text-align: center;">
        <h2 style="color: #059669; font-size: 1.75rem; font-weight: 700; margin-bottom: 2rem;">中药化合物智能鉴定平台</h2>
    </div>
    """, unsafe_allow_html=True)

    with st.form("login_form"):
        username = st.text_input("用户名", placeholder="请输入用户名", help="默认用户名: ZY")
        password = st.text_input("密码", type="password", placeholder="请输入密码", help="默认密码: 513513")
        submitted = st.form_submit_button("登录", use_container_width=True)

        if submitted:
            if username == VALID_USERNAME and password == VALID_PASSWORD:
                st.session_state.logged_in = True
                st.session_state.username = username
                st.success("登录成功！正在跳转...")
                time.sleep(0.5)
                st.rerun()
            else:
                st.error("用户名或密码错误，请重试。")

    st.markdown("""
        <div style="margin-top: 2rem; color: #94a3b8; font-size: 0.85rem; text-align: center;">
            <p>中药化合物智能鉴定平台 v5.13</p>
            <p>支持主数据库 + 英文数据库联合鉴定 | 文献-碎片详细映射</p>
        </div>
    </div>
    """, unsafe_allow_html=True)


def logout_button():
    if st.sidebar.button("🚪 登出", use_container_width=True):
        st.session_state.logged_in = False
        st.session_state.pop('username', None)
        st.rerun()


def create_header():
    username = st.session_state.get('username', 'Guest')
    st.markdown(f"""
    <div class="main-header">
        <h1>🌿 中药化合物智能鉴定平台</h1>
        <p>v5.13 文献-碎片详细映射 | 一级ppm/二级Da双容差控制 | 欢迎回来，{username}</p>
    </div>
    """, unsafe_allow_html=True)


def create_sidebar():
    st.sidebar.markdown("""
    <div style="text-align: center; padding: 1.5rem 0; border-bottom: 2px solid #e2e8f0;">
        <h3 style="color: #1e293b; margin: 0.5rem 0; font-weight: 700;">TCM Identifier</h3>
        <p style="color: #64748b; font-size: 0.8rem; margin: 0;">中药化合物鉴定系统</p>
    </div>
    """, unsafe_allow_html=True)

    if st.session_state.get('logged_in'):
        st.sidebar.markdown(f"""
        <div style="background: linear-gradient(135deg, #ecfdf5, #d1fae5); padding: 1rem; border-radius: 12px; margin: 1rem 0; text-align: center;">
            <p style="margin: 0; color: #059669; font-weight: 600;">👤 {st.session_state.username}</p>
        </div>
        """, unsafe_allow_html=True)

    page = st.sidebar.radio("选择功能", ["🏠 首页", "🔬 开始鉴定", "🧪 诊断离子筛查", "📖 使用指南", "🗃️ 数据库预览", "📊 结果分析"], index=0)
    page = page.split(' ', 1)[1] if ' ' in page else page

    st.sidebar.markdown("---")
    st.sidebar.markdown("""
    <div style="background: white; border-radius: 12px; padding: 1rem; box-shadow: 0 2px 10px rgba(0,0,0,0.05);">
        <h4 style="color: #1e293b; margin: 0 0 0.75rem 0; font-size: 0.95rem;">📋 版本信息</h4>
        <ul style="color: #64748b; font-size: 0.85rem; padding-left: 1.2rem; margin: 0;">
            <li>程序版本：v5.13</li>
            <li>主数据库：已加载</li>
            <li>一级匹配：ppm控制</li>
            <li>二级匹配：绝对Da容差</li>
            <li>文献-碎片映射：✅ 已启用</li>
        </ul>
    </div>
    """, unsafe_allow_html=True)

    return page


def show_home_page():
    create_header()
    st.markdown("## 📊 系统概览")
    
    cols = st.columns(4)
    stats_data = [("35,828+", "数据库总记录数", "📚"), ("291", "支持药材种类", "🌱"), ("6", "鉴定评级级别", "⭐"), ("400+", "化合物类型", "🧪")]
    for col, (number, label, icon) in zip(cols, stats_data):
        with col:
            st.markdown(f"""
            <div class="stat-card">
                <div style="font-size: 2rem; margin-bottom: 0.5rem;">{icon}</div>
                <div class="stat-number">{number}</div>
                <div class="stat-label">{label}</div>
            </div>
            """, unsafe_allow_html=True)

    st.markdown("---")
    st.markdown("## ✨ 核心功能")
    
    col1, col2 = st.columns(2)
    with col1:
        st.markdown("""
        <div class="feature-card">
            <h4>🔍 智能化合物鉴定</h4>
            <p>基于高分辨质谱数据，在主数据库和英文数据库中进行精准匹配，支持正负离子模式。</p>
        </div>
        <div class="feature-card">
            <h4>📈 六级评级标准</h4>
            <p>确证级、高置信级、推定级、提示级、参考级、排除级，科学评估鉴定结果可靠性。</p>
        </div>
        <div class="feature-card">
            <h4>🧪 外部诊断离子</h4>
            <p>支持上传自定义诊断离子文件，自动去重合并权重，避免重复。</p>
        </div>
        """, unsafe_allow_html=True)
    with col2:
        st.markdown("""
        <div class="feature-card">
            <h4>⚡ 双层次匹配策略</h4>
            <p>一级匹配（准分子离子）使用ppm容差控制；二级匹配（碎片离子）使用绝对Da容差控制，可独立调整。</p>
        </div>
        <div class="feature-card">
            <h4>📚 文献-碎片详细映射</h4>
            <p>支持碎片对应多篇文献，输出格式：139.05(文献1,文献2); 134.03(文献1,文献3)</p>
        </div>
        <div class="feature-card">
            <h4>🌱 药材来源合并</h4>
            <p>同一化合物的所有药材来源自动合并显示，实现"所有的药材"。</p>
        </div>
        """, unsafe_allow_html=True)

    if st.button("🚀 立即开始鉴定", type="primary", use_container_width=True):
        st.session_state['page'] = '开始鉴定'
        st.rerun()


def show_analysis_page():
    """鉴定分析页面"""
    create_header()

    st.markdown("## 📁 上传质谱数据")
    st.markdown("**至少上传一个文件（正离子或负离子）**")

    col1, col2 = st.columns(2)
    with col1:
        ms_positive_file = st.file_uploader("上传正离子模式质谱数据 (.xlsx，可选)", type=['xlsx'], key='ms_positive')
        if ms_positive_file:
            st.success(f"✅ 已上传: {ms_positive_file.name}")
    with col2:
        ms_negative_file = st.file_uploader("上传负离子模式质谱数据 (.xlsx，可选)", type=['xlsx'], key='ms_negative')
        if ms_negative_file:
            st.success(f"✅ 已上传: {ms_negative_file.name}")

    if not ms_positive_file and not ms_negative_file:
        st.warning("⚠️ 请至少上传一个质谱数据文件")

    st.markdown("---")
    st.markdown("## 🧪 外部诊断离子（可选）")

    diagnostic_file = st.file_uploader("上传自定义诊断离子文件 (.xlsx，可包含权重列)", type=['xlsx'], key='diagnostic')
    if diagnostic_file:
        st.success(f"✅ 已上传诊断离子文件: {diagnostic_file.name}")

    st.markdown("---")
    st.markdown("## 📚 自定义数据库（可选）")

    custom_db_file = st.file_uploader("上传自定义数据库 (.xlsx，需与主数据库列一致)", type=['xlsx'], key='custom_db')
    if custom_db_file:
        st.success(f"✅ 已上传自定义数据库: {custom_db_file.name}")

    st.markdown("---")
    st.markdown("## 🌍 数据库状态")
    
    english_db_path = find_english_database_path()
    if english_db_path:
        try:
            english_df = pd.read_excel(english_db_path)
            st.success(f"✅ 英文数据库已找到: {english_db_path}，共 {len(english_df)} 条记录")
        except Exception as e:
            st.warning(f"⚠️ 英文数据库文件存在但读取失败: {e}")
    else:
        st.info("ℹ️ 未找到英文数据库文件，将仅使用主数据库进行鉴定")
    
    st.info("💡 文献-碎片映射：碎片离子列支持分号(;)、逗号(,)、顿号(、)、空格分隔，输出格式: 139.05(文献1,文献2); 134.03(文献1,文献3)")

    st.markdown("---")
    st.markdown("## ⚙️ 鉴定参数配置")

    preset = st.radio("参数预设", options=["自定义", "快速模式", "高精度模式"], horizontal=True)

    col1, col2, col3, col4 = st.columns(4)
    with col1:
        tolerance_ppm = st.number_input("一级ppm误差容限", min_value=10, max_value=100, value=50, help="准分子离子匹配使用ppm容差")
    with col2:
        fragment_tolerance_no_ppm = st.number_input("二级绝对Da容差", min_value=0.01, max_value=1.0, value=0.1, step=0.01, help="碎片离子匹配使用绝对Da容差")
    with col3:
        max_candidates = st.number_input("最大候选数", min_value=1, max_value=10, value=3)
    with col4:
        use_parallel = st.checkbox("启用多进程并行", value=True)

    col1, col2, col3 = st.columns(3)
    with col1:
        tolerance_type = st.selectbox("诊断离子容差类型", options=["Da", "ppm"], index=0)
        if tolerance_type == "Da":
            fragment_tolerance = st.number_input("诊断离子容差 (Da)", min_value=0.01, max_value=1.0, value=0.05, step=0.01)
            fragment_tolerance_ppm = 20
        else:
            fragment_tolerance_ppm = st.number_input("诊断离子容差 (ppm)", min_value=1, max_value=100, value=20, step=1)
            fragment_tolerance = 0.05
    with col2:
        rt_tolerance = st.number_input("保留时间容差 (min)", min_value=0.1, max_value=2.0, value=0.3, step=0.1)
    with col3:
        loss_tolerance = st.number_input("中性丢失容差 (Da)", min_value=0.01, max_value=0.5, value=0.02, step=0.01)

    col1, col2, col3 = st.columns(3)
    with col1:
        intensity_rel_threshold = st.slider("相对强度阈值 (%)", min_value=0.0, max_value=100.0, value=1.0, step=0.1)
    with col2:
        use_rt_score = st.checkbox("启用保留时间得分", value=True)
    with col3:
        cache_index = st.checkbox("启用索引缓存（加速）", value=True)

    herb_filter = st.selectbox("药材筛选模式", options=["使用全部数据库", "筛选特定药材"], index=0)
    if herb_filter == "筛选特定药材":
        herb_name = st.text_input("输入药材名称", placeholder="如：栀子、黄芩")
    else:
        herb_name = None

    # 初始化变量，确保在预设模式中定义
    min_intensity_abs = 100  # 默认值
    
    if preset == "快速模式":
        tolerance_ppm = 50
        fragment_tolerance_no_ppm = 0.1
        fragment_tolerance = 0.05
        min_intensity_abs = 100
        intensity_rel_threshold = 1.0
    elif preset == "高精度模式":
        tolerance_ppm = 10
        fragment_tolerance_no_ppm = 0.02
        fragment_tolerance = 0.02
        min_intensity_abs = 50
        intensity_rel_threshold = 0.5

    col1, col2 = st.columns(2)
    with col1:
        if st.button("💾 保存当前参数", use_container_width=True):
            st.session_state['saved_params'] = {
                'tolerance_ppm': tolerance_ppm,
                'fragment_tolerance_no_ppm': fragment_tolerance_no_ppm,
                'fragment_tolerance': fragment_tolerance,
                'fragment_tolerance_ppm': fragment_tolerance_ppm,
                'tolerance_type': tolerance_type,
                'min_intensity_abs': min_intensity_abs,
                'intensity_rel_threshold': intensity_rel_threshold,
                'rt_tolerance': rt_tolerance,
                'loss_tolerance': loss_tolerance,
                'use_rt_score': use_rt_score,
                'cache_index': cache_index,
            }
            st.success("参数已保存")
    with col2:
        if 'saved_params' in st.session_state and st.button("🔄 加载已保存参数", use_container_width=True):
            params = st.session_state['saved_params']
            tolerance_ppm = params['tolerance_ppm']
            fragment_tolerance_no_ppm = params['fragment_tolerance_no_ppm']
            fragment_tolerance = params['fragment_tolerance']
            fragment_tolerance_ppm = params['fragment_tolerance_ppm']
            tolerance_type = params['tolerance_type']
            min_intensity_abs = params['min_intensity_abs']
            intensity_rel_threshold = params['intensity_rel_threshold']
            rt_tolerance = params['rt_tolerance']
            loss_tolerance = params['loss_tolerance']
            use_rt_score = params['use_rt_score']
            cache_index = params['cache_index']
            st.rerun()

    st.markdown("---")

    if ms_positive_file or ms_negative_file:
        if st.button("🚀 开始化合物鉴定", type="primary", use_container_width=True):
            with st.spinner("正在初始化鉴定程序..."):
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
                    if ms_negative_file:
                        neg_path = os.path.join(temp_dir, ms_negative_file.name)
                        with open(neg_path, 'wb') as f:
                            f.write(ms_negative_file.getbuffer())
                    if diagnostic_file:
                        diag_path = os.path.join(temp_dir, diagnostic_file.name)
                        with open(diag_path, 'wb') as f:
                            f.write(diagnostic_file.getbuffer())
                    if custom_db_file:
                        custom_db_path = os.path.join(temp_dir, custom_db_file.name)
                        with open(custom_db_path, 'wb') as f:
                            f.write(custom_db_file.getbuffer())

                    db_path = find_database_path()
                    if not db_path:
                        st.error("未找到主数据库文件！请将 TCM-SM-MS DB.xlsx 放在项目目录下。")
                        return

                    english_db_path = find_english_database_path()
                    if english_db_path:
                        st.info(f"✅ 已找到主数据库文件: {db_path}")
                        st.info(f"✅ 已找到英文数据库文件: {english_db_path}")
                    else:
                        st.info(f"✅ 已找到主数据库文件: {db_path}")
                        st.warning("⚠️ 未找到英文数据库文件，将仅使用主数据库")

                    progress_bar = st.progress(0)
                    status_text = st.empty()

                    status_text.text("【1/9】正在加载数据库...")
                    progress_bar.progress(10)

                    config = {
                        'min_intensity': min_intensity_abs,
                        'fragment_tolerance': fragment_tolerance,
                        'fragment_tolerance_ppm': fragment_tolerance_ppm,
                        'fragment_tolerance_no_ppm': fragment_tolerance_no_ppm,
                        'tolerance_ppm': tolerance_ppm,
                        'max_candidates': max_candidates,
                    }

                    identifier = UltimateGardeniaIdentifier(
                        database_path=db_path,
                        ms_positive_path=pos_path,
                        ms_negative_path=neg_path,
                        herb_name=herb_name,
                        config=config,
                        use_parallel=use_parallel,
                        rt_tolerance=rt_tolerance,
                        loss_tolerance=loss_tolerance,
                        external_diagnostic_file=diag_path,
                        rt_fusion_tolerance=rt_tolerance,
                        intensity_relative_threshold=intensity_rel_threshold,
                        tolerance_type=tolerance_type,
                        use_rt_score=use_rt_score,
                        custom_db_path=custom_db_path,
                        cache_index=cache_index,
                        english_db_path=english_db_path,
                        fragment_tolerance_no_ppm=fragment_tolerance_no_ppm
                    )

                    status_text.text("【7/9】正在处理质谱数据...")
                    progress_bar.progress(40)

                    report = identifier.generate_report('样品')

                    progress_bar.progress(80)
                    status_text.text("【8/9】生成报告...")

                    st.session_state['analysis_results'] = report
                    st.session_state['identifier'] = identifier

                    progress_bar.progress(100)
                    status_text.text("鉴定完成！")

                    st.success("✅ 化合物鉴定完成！")

                    if st.button("📊 查看鉴定结果 →"):
                        st.session_state['page'] = '结果分析'
                        st.rerun()

                except Exception as e:
                    st.error(f"鉴定过程中出错：{str(e)}")
                    st.exception(e)
    else:
        st.markdown("""
        <div class="empty-state">
            <h3>暂无上传文件</h3>
            <p>请上传质谱数据文件以开始鉴定</p>
        </div>
        """, unsafe_allow_html=True)


def show_diagnostic_ion_page():
    """诊断离子筛查页面"""
    create_header()
    st.markdown("## 🔬 诊断离子筛查")
    st.markdown("根据输入的m/z值，在诊断离子数据库中查找匹配的化合物特征离子")

    diagnostic_df = load_diagnostic_ions_cached()

    if diagnostic_df.empty:
        st.warning("⚠️ 未找到诊断离子数据库文件（诊断离子.xlsx），请将文件放在项目目录下。")
        return

    st.success(f"✅ 已加载诊断离子数据库，包含 {len(diagnostic_df)} 条记录")

    st.markdown("---")
    st.markdown("### 📥 输入m/z值")

    mz_input = st.text_area("输入m/z值（每行一个值，或用逗号分隔）", placeholder="例如：\n151.003\n137.024, 121.029", height=150)
    tolerance_ppm = st.number_input("ppm容差", min_value=1, max_value=100, value=10)

    def parse_mz_values(input_text):
        if not input_text or not input_text.strip():
            return []
        for sep in [',', ';', '\t', ' ']:
            if sep in input_text:
                input_text = input_text.replace(sep, '\n')
        mz_values = []
        for line in input_text.strip().split('\n'):
            line = line.strip()
            if line:
                try:
                    mz_values.append(float(line))
                except ValueError:
                    continue
        return mz_values

    if mz_input:
        user_mz_values = parse_mz_values(mz_input)
        if not user_mz_values:
            st.warning("⚠️ 无法解析输入的m/z值，请检查输入格式。")
        else:
            with st.spinner("正在匹配诊断离子..."):
                results_df = match_diagnostic_ions(user_mz_values, diagnostic_df, tolerance_ppm=tolerance_ppm)

            st.markdown("### 📊 匹配结果统计")
            cols = st.columns(3)
            with cols[0]:
                st.metric("输入离子数", len(user_mz_values))
            with cols[1]:
                matched_ions = results_df['输入m/z'].nunique() if not results_df.empty else 0
                st.metric("匹配离子数", matched_ions)
            with cols[2]:
                matched_compounds = results_df['化合物类型'].nunique() if not results_df.empty else 0
                st.metric("化合物类型数", matched_compounds)

            if not results_df.empty:
                st.markdown("### 📋 匹配结果详情")
                st.dataframe(results_df, use_container_width=True, hide_index=True)
            else:
                st.info("未找到匹配的诊断离子，请尝试增大ppm容差。")


def show_guide_page():
    """使用指南页面"""
    create_header()
    st.markdown("## 📖 使用指南（v5.13 文献-碎片详细映射）")
    st.markdown("""
    <div style="background: white; border-radius: 16px; padding: 2rem; box-shadow: 0 4px 20px rgba(0,0,0,0.08);">
        <h3>🎯 匹配策略</h3>
        <p><strong>一级匹配（准分子离子）</strong>：使用ppm容差控制，筛选候选化合物</p>
        <p><strong>二级匹配（碎片离子）</strong>：使用绝对Da容差控制，在一级匹配结果中进行碎片比对</p>
        <p><strong>文献-碎片映射</strong>：每个碎片自动绑定其来源文献，输出格式 <code>139.05(文献1,文献2); 134.03(文献1,文献3)</code></p>

        <h3>📊 综合评分规则</h3>
        <ul>
            <li>基础分：1级85分、2级65分、3级45分、4级25分、5级0分</li>
            <li>ppm调整：≤5ppm+5分、5-10ppm0分、10-20ppm-5分、20-30ppm-10分、30-50ppm-15分</li>
            <li>匹配碎片：每个+2分（上限20）</li>
            <li>诊断离子：根据权重累计（上限15）</li>
            <li>中性丢失：每个+2分（上限10）</li>
            <li>RT偏差<0.2：+5分，RT偏差<0.5：+2分</li>
        </ul>
    </div>
    """, unsafe_allow_html=True)


def show_database_page():
    """数据库预览页面"""
    create_header()
    st.markdown("## 🗃️ 数据库预览")

    db_path = find_database_path()
    if not db_path:
        st.warning("⚠️ 未找到主数据库文件！请将 TCM-SM-MS DB.xlsx 放在项目目录下。")
        return

    try:
        df = load_database_cached()
        if df.empty:
            st.error("主数据库文件为空或无法读取！")
            return

        st.success(f"✅ 成功加载主数据库，包含 {len(df)} 条化合物记录")
        st.dataframe(df.head(10), use_container_width=True)

        st.markdown("---")
        st.markdown("### 🌍 英文数据库")
        english_db_path = find_english_database_path()
        if english_db_path:
            english_df = pd.read_excel(english_db_path)
            st.success(f"✅ 英文数据库已加载: {english_db_path}，共 {len(english_df)} 条记录")
            st.dataframe(english_df.head(10), use_container_width=True)
        else:
            st.info("ℹ️ 未找到英文数据库文件")

    except Exception as e:
        st.error(f"加载数据库时出错：{str(e)}")


def show_results_page():
    """结果分析页面"""
    create_header()

    if 'analysis_results' not in st.session_state:
        st.markdown("""
        <div class="empty-state">
            <h3>暂无鉴定结果</h3>
            <p>请先进行化合物鉴定</p>
        </div>
        """, unsafe_allow_html=True)
        if st.button("前往鉴定页面"):
            st.session_state['page'] = '开始鉴定'
            st.rerun()
        return

    report = st.session_state['analysis_results']

    st.markdown("## 📊 鉴定结果分析")

    if report.empty:
        st.warning("鉴定结果为空，可能是因为没有匹配的化合物。")
        return

    cols = st.columns(4)
    with cols[0]:
        st.metric("鉴定化合物总数", len(report))
    with cols[1]:
        confirmed = (report['评级名称'] == '确证级').sum()
        st.metric("确证级化合物", confirmed)
    with cols[2]:
        high_score = (report['综合得分'] >= 90).sum()
        st.metric("90分以上", high_score)
    with cols[3]:
        avg_score = report['综合得分'].mean()
        st.metric("平均综合得分", f"{avg_score:.1f}")

    st.markdown("### 📋 完整鉴定结果")

    all_columns = report.columns.tolist()
    default_cols = ['序号', '化合物中文名', '分子式', 'ppm', '评级名称', '药材名称', '主要碎片离子(文献映射)', '综合得分']
    selected_cols = st.multiselect("选择显示的列", all_columns, default=[c for c in default_cols if c in all_columns])

    display_df = report[selected_cols] if selected_cols else report

    st.dataframe(display_df, use_container_width=True, hide_index=True)

    st.markdown("---")
    st.markdown("### 📥 导出报告")

    col1, col2 = st.columns(2)
    with col1:
        csv = report.to_csv(index=False, encoding='utf-8-sig')
        st.download_button(label="📥 导出CSV", data=csv, file_name=f"鉴定报告_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv", mime="text/csv")
    with col2:
        buffer = BytesIO()
        with pd.ExcelWriter(buffer, engine='openpyxl') as writer:
            report.to_excel(writer, index=False, sheet_name='鉴定结果')
        st.download_button(label="📥 导出Excel", data=buffer.getvalue(), file_name=f"鉴定报告_{datetime.now().strftime('%Y%m%d_%H%M%S')}.xlsx", mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")


def main():
    load_optimized_css()

    if 'logged_in' not in st.session_state:
        st.session_state.logged_in = False

    if not st.session_state.logged_in:
        login_page()
        return

    if 'page' not in st.session_state:
        st.session_state['page'] = '首页'

    logout_button()
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
