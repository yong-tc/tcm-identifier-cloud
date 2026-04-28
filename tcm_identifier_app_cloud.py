# -*- coding: utf-8 -*-
"""
中药化合物智能鉴定平台 - 云端部署增强版 v6.1（三库独立比对版）
==========================================================
更新说明：
- 修复化合物名称与准分子离子不匹配的问题
- 增强母离子匹配验证机制
- 添加匹配确认步骤确保化合物信息与母离子对应
"""

import streamlit as st
import pandas as pd
import numpy as np
import os
import time
import pickle
import hashlib
from datetime import datetime
from io import BytesIO
import tempfile
from bisect import bisect_left, bisect_right
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
# 碎片离子解析函数（支持多种分隔符）
# ============================================================================
def parse_fragments(fragment_string):
    """
    解析碎片离子字符串，支持顿号、分号、逗号、空格、换行符分隔符
    返回：碎片离子列表（float）
    """
    if pd.isna(fragment_string) or not fragment_string or str(fragment_string).strip() == '':
        return []
    s = str(fragment_string)
    for sep in ['、', '；', ',', ' ', '\t', '\n', '\r\n']:
        if sep in s:
            s = s.replace(sep, ';')
    
    fragments = []
    for part in s.split(';'):
        part = part.strip()
        if part:
            if '<br>' in part:
                sub_parts = part.split('<br>')
                for sub in sub_parts:
                    sub = sub.strip()
                    if sub:
                        try:
                            fragments.append(float(sub))
                        except ValueError:
                            continue
            else:
                try:
                    fragments.append(float(part))
                except ValueError:
                    continue
    return fragments


def parse_fragments_with_source(fragment_string, source_name, db_source=''):
    """
    解析碎片离子字符串，并记录来源（包括数据库来源和文献来源）
    返回：碎片离子列表（float），以及来源映射
    """
    fragments = parse_fragments(fragment_string)
    source_map = {}
    for frag in fragments:
        if frag not in source_map:
            source_map[frag] = set()
        # 记录文献来源
        if source_name and source_name != 'nan':
            source_map[frag].add(f"[{db_source}] {source_name}" if db_source else source_name)
        # 如果没有文献来源，至少记录数据库来源
        elif db_source:
            source_map[frag].add(f"[{db_source}]")
    return fragments, source_map


# ============================================================================
# 数据库加载函数（带缓存）- 支持多个数据库
# ============================================================================

@st.cache_data
def load_database_cached(db_filename=None):
    """加载主数据库（TCM-SM-MS DB）"""
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
                column_mapping = {
                    '药材名': '药材名称',
                    '文献': '文献来源',
                    '保留时间(min)': '保留时间',
                    '中丢失': '中性丢失'
                }
                for old_name, new_name in column_mapping.items():
                    if old_name in df.columns and new_name not in df.columns:
                        df.rename(columns={old_name: new_name}, inplace=True)
                df['_data_source'] = '主数据库'
                return df
            except Exception as e:
                continue
    return pd.DataFrame()


@st.cache_data
def load_english_database_cached(db_filename="数据库（英文）.xlsx"):
    """加载英文数据库"""
    db_paths = [
        "数据库（英文）.xlsx",
        "data/数据库（英文）.xlsx",
        "user_input_files/数据库（英文）.xlsx"
    ]
    if db_filename and os.path.exists(db_filename):
        db_paths.insert(0, db_filename)
    
    for path in db_paths:
        if os.path.exists(path):
            try:
                df = pd.read_excel(path)
                if '名称（中文）' in df.columns and '准分子离子（正）' in df.columns:
                    column_mapping = {
                        '名称（中文）': '名称（中文）',
                        '名称（英文）': '名称（英文）',
                        '分子式': '分子式',
                        '化合物类型': '化合物类型',
                        '加合物（正）': '加合物（正）',
                        '加合物（负）': '加合物（负）',
                        '准分子离子（正）': '准分子离子（正）',
                        '碎片离子（正）': '碎片离子（正）',
                        '准分子离子（负）': '准分子离子（负）',
                        '碎片离子（负）': '碎片离子（负）',
                        '药材名': '药材名称',
                        '文献': '文献来源'
                    }
                    for old_name, new_name in column_mapping.items():
                        if old_name in df.columns and new_name not in df.columns:
                            df.rename(columns={old_name: new_name}, inplace=True)
                    
                    if '药材名称' not in df.columns:
                        df['药材名称'] = ''
                    if '文献来源' not in df.columns:
                        df['文献来源'] = ''
                    
                    df['_data_source'] = '英文数据库'
                    print(f"  英文数据库加载成功: {len(df)} 条记录")
                    return df
            except Exception as e:
                print(f"  加载英文数据库失败: {path}, 错误: {e}")
                continue
    print("  未找到英文数据库文件（数据库（英文）.xlsx）")
    return pd.DataFrame()


@st.cache_data
def load_standard_database_cached(db_filename="对照品数据库.xlsx"):
    """加载对照品数据库"""
    db_paths = [
        "对照品数据库.xlsx",
        "data/对照品数据库.xlsx",
        "user_input_files/对照品数据库.xlsx"
    ]
    if db_filename and os.path.exists(db_filename):
        db_paths.insert(0, db_filename)
    
    for path in db_paths:
        if os.path.exists(path):
            try:
                df = pd.read_excel(path, sheet_name=0)
                column_mapping = {
                    '中文名': '名称（中文）',
                    'Name': '名称（英文）',
                    'Formula': '分子式',
                    '加合物': '加合物（正）',
                    '准分子离子（正）': '准分子离子（正）',
                    '碎片离子（正）': '碎片离子（正）',
                    'Rel Abund': '相对丰度',
                    '准分子离子（负）': '准分子离子（负）',
                    '碎片离子（负）': '碎片离子（负）'
                }
                for old_name, new_name in column_mapping.items():
                    if old_name in df.columns and new_name not in df.columns:
                        df.rename(columns={old_name: new_name}, inplace=True)
                
                df['_data_source'] = '对照品数据库'
                df['_is_standard'] = True
                print(f"  对照品数据库加载成功: {len(df)} 条记录")
                return df
            except Exception as e:
                print(f"  加载对照品数据库失败: {path}, 错误: {e}")
                continue
    print("  未找到对照品数据库文件（对照品数据库.xlsx）")
    return pd.DataFrame()


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
    db_paths = [
        "数据库（英文）.xlsx",
        "data/数据库（英文）.xlsx",
        "user_input_files/数据库（英文）.xlsx"
    ]
    for path in db_paths:
        if os.path.exists(path):
            return path
    return None


def find_standard_database_path():
    """查找对照品数据库路径"""
    db_paths = [
        "对照品数据库.xlsx",
        "data/对照品数据库.xlsx",
        "user_input_files/对照品数据库.xlsx"
    ]
    for path in db_paths:
        if os.path.exists(path):
            return path
    return None


@st.cache_data
def load_diagnostic_ions_cached():
    """加载诊断离子数据库"""
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
    """查找诊断离子文件路径"""
    diagnostic_ion_paths = [
        "诊断离子.xlsx",
        "data/诊断离子.xlsx",
        "user_input_files/诊断离子.xlsx"
    ]
    for path in diagnostic_ion_paths:
        if os.path.exists(path):
            return path
    return None


# ============================================================================
# 鉴定程序核心代码（三库独立比对版 - 修复母离子匹配）
# ============================================================================

class UltimateGardeniaIdentifier:
    """
    中药化合物鉴定终极版程序 v6.1（三库独立比对版）
    
    主要功能：
    1. 三个数据库（主数据库、英文数据库、对照品数据库）分别独立比对
    2. 对照品数据库文件扩展名改为 .xlsx
    3. 在"数据来源"列中标注"主数据库/英文数据库/对照品数据库"
    4. 碎片离子映射所有文献来源，包括3个数据库的
    5. 碎片离子标注来源于哪个数据库
    6. 修复化合物名称与准分子离子不匹配的问题
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
                 english_db_path=None,
                 standard_db_path=None,
                 cache_index=True):
        """初始化鉴定程序（三库独立比对版）"""
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
        self.english_db_path = english_db_path
        self.standard_db_path = standard_db_path
        self.cache_index = cache_index
        self.herb_name = herb_name

        self.use_parallel = use_parallel and os.cpu_count() > 1
        self.num_workers = min(os.cpu_count(), 8)

        print("="*80)
        print("中药化合物鉴定程序 v6.1（三库独立比对版 - 修复母离子匹配）")
        print("="*80)
        
        print("\n【1/9】正在加载数据库...")
        
        # 加载主数据库
        self.main_database = self._load_data(database_path)
        if self.main_database.empty:
            print("  警告: 主数据库为空或未找到")
        else:
            print(f"  主数据库加载成功: {len(self.main_database)} 条记录")
        
        # 加载英文数据库
        self.english_database = self._load_english_data(english_db_path)
        if not self.english_database.empty:
            print(f"  英文数据库加载成功: {len(self.english_database)} 条记录")
        
        # 加载对照品数据库
        self.standard_database = self._load_standard_data(standard_db_path)
        if not self.standard_database.empty:
            print(f"  对照品数据库加载成功: {len(self.standard_database)} 条记录")
        else:
            print("  对照品数据库未找到或为空")

        # 加载自定义数据库
        if custom_db_path and os.path.exists(custom_db_path):
            custom_db = self._load_data(custom_db_path)
            if not custom_db.empty:
                print(f"  加载自定义数据库: {len(custom_db)} 条记录")
                # 自定义数据库添加到主数据库
                if self.main_database.empty:
                    self.main_database = custom_db
                else:
                    self.main_database = pd.concat([self.main_database, custom_db], ignore_index=True)

        # 筛选药材（如果指定）
        if herb_name:
            print(f"【2/9】正在筛选 {herb_name} 相关数据...")
            self.main_database = self._filter_by_herb(self.main_database, herb_name)
            self.english_database = self._filter_by_herb(self.english_database, herb_name)
            self.standard_database = self._filter_by_herb(self.standard_database, herb_name)
        else:
            print("【2/9】使用全部数据库进行化合物鉴定")

        print(f"  主数据库筛选后: {len(self.main_database)} 条")
        print(f"  英文数据库筛选后: {len(self.english_database)} 条")
        print(f"  对照品数据库筛选后: {len(self.standard_database)} 条")

        print("【3/9】正在加载质谱数据...")
        self.ms_positive = self._load_data(ms_positive_path) if ms_positive_path else pd.DataFrame()
        self.ms_negative = self._load_data(ms_negative_path) if ms_negative_path else pd.DataFrame()
        print(f"  正离子数据: {len(self.ms_positive)} 条记录")
        print(f"  负离子数据: {len(self.ms_negative)} 条记录")

        print("【4/9】正在构建索引...")
        # 为每个数据库分别构建索引
        self._build_separate_indices()

        print("【5/9】正在加载诊断离子库...")
        if external_diagnostic_file is None:
            default_diag_path = find_diagnostic_ion_path()
            if default_diag_path:
                print(f"  自动找到默认诊断离子文件: {default_diag_path}")
                external_diagnostic_file = default_diag_path
        self._build_diagnostic_ion_library(external_diagnostic_file)

        print("【6/9】正在加载辅助数据...")
        self._load_auxiliary_data()

        self.stats = {
            'total_precursors': 0,
            'identified_compounds': 0,
            'scored_candidates': 0
        }

        self._print_initialization_info()
        print("【7/9】初始化完成，准备鉴定")

    def _load_data(self, filepath):
        """加载数据文件（Excel/CSV）"""
        if filepath and os.path.exists(filepath):
            try:
                if filepath.endswith('.xlsx'):
                    df = pd.read_excel(filepath)
                elif filepath.endswith('.csv'):
                    return pd.read_csv(filepath)
                else:
                    return pd.DataFrame()

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
            except Exception as e:
                print(f"警告: 无法加载文件 {filepath}: {e}")
                return pd.DataFrame()
        return pd.DataFrame()

    def _load_english_data(self, filepath):
        """加载英文数据库"""
        if filepath and os.path.exists(filepath):
            try:
                df = pd.read_excel(filepath)
                column_mapping = {
                    '名称（中文）': '名称（中文）',
                    '名称（英文）': '名称（英文）',
                    '分子式': '分子式',
                    '化合物类型': '化合物类型',
                    '加合物（正）': '加合物（正）',
                    '加合物（负）': '加合物（负）',
                    '准分子离子（正）': '准分子离子（正）',
                    '碎片离子（正）': '碎片离子（正）',
                    '准分子离子（负）': '准分子离子（负）',
                    '碎片离子（负）': '碎片离子（负）',
                    '药材名': '药材名称',
                    '文献': '文献来源'
                }
                for old_name, new_name in column_mapping.items():
                    if old_name in df.columns and new_name not in df.columns:
                        df.rename(columns={old_name: new_name}, inplace=True)
                
                df['_data_source'] = '英文数据库'
                return df
            except Exception as e:
                print(f"警告: 无法加载英文数据库 {filepath}: {e}")
                return pd.DataFrame()
        return pd.DataFrame()

    def _load_standard_data(self, filepath):
        """加载对照品数据库"""
        if filepath and os.path.exists(filepath):
            try:
                df = pd.read_excel(filepath, sheet_name=0)
                column_mapping = {
                    '中文名': '名称（中文）',
                    'Name': '名称（英文）',
                    'Formula': '分子式',
                    '加合物': '加合物（正）',
                    '准分子离子（正）': '准分子离子（正）',
                    '碎片离子（正）': '碎片离子（正）',
                    'Rel Abund': '相对丰度',
                    '准分子离子（负）': '准分子离子（负）',
                    '碎片离子（负）': '碎片离子（负）'
                }
                for old_name, new_name in column_mapping.items():
                    if old_name in df.columns and new_name not in df.columns:
                        df.rename(columns={old_name: new_name}, inplace=True)
                
                if '名称（中文）' not in df.columns:
                    df['名称（中文）'] = ''
                if '名称（英文）' not in df.columns:
                    df['名称（英文）'] = ''
                if '分子式' not in df.columns:
                    df['分子式'] = ''
                if '化合物类型' not in df.columns:
                    df['化合物类型'] = ''
                if '药材名称' not in df.columns:
                    df['药材名称'] = ''
                
                df['_data_source'] = '对照品数据库'
                df['_is_standard'] = True
                return df
            except Exception as e:
                print(f"警告: 无法加载对照品数据库 {filepath}: {e}")
                return pd.DataFrame()
        return pd.DataFrame()

    def _build_diagnostic_ion_library(self, external_file=None):
        """构建诊断性离子库（支持外部Excel文件，自动去重合并权重）"""
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
                print(f"  成功加载外部诊断离子库：{len(self.diagnostic_ions)} 类，{len(df)} 条记录（去重后）")
            except Exception as e:
                print(f"  加载外部诊断离子文件失败：{e}，将使用内置库")
                self._build_default_diagnostic_ions()
        else:
            print("  未找到外部诊断离子文件，使用内置诊断离子库（无权重，默认权重1）")
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
        """加载辅助数据"""
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

    def _classify_compound(self, name, compound_type):
        """化合物类型分类"""
        name = name.lower()
        compound_type = str(compound_type).lower()

        if any(keyword in compound_type for keyword in ['环烯醚', 'iridoid']):
            return '环烯醚萜类'
        if any(keyword in compound_type for keyword in ['有机酸', '酚酸', '有机酸类']):
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

    def _filter_by_herb(self, database, herb_name):
        """根据药材名称筛选数据库"""
        if herb_name is None or database.empty:
            return database.copy() if not database.empty else pd.DataFrame()

        herb_col = '药材名称' if '药材名称' in database.columns else '药材名'
        if herb_col not in database.columns:
            return database.copy()
        
        mask = database[herb_col].str.contains(herb_name, na=False, case=False)
        filtered_db = database[mask].copy()
        if len(filtered_db) == 0:
            print(f"警告: 未找到 '{herb_name}' 相关数据，将使用全部数据")
            return database.copy()
        return filtered_db

    def _build_separate_indices(self):
        """为每个数据库分别构建索引"""
        # 主数据库索引
        self.main_index = self._build_single_index(self.main_database, '主数据库')
        
        # 英文数据库索引
        self.english_index = self._build_single_index(self.english_database, '英文数据库')
        
        # 对照品数据库索引
        self.standard_index = self._build_single_index(self.standard_database, '对照品数据库')
        
        # 合并所有碎片来源（用于统一映射）
        self._build_global_fragment_sources()

    def _build_single_index(self, database, db_name):
        """为单个数据库构建索引"""
        index_data = {
            'sorted_idx_pos': [],
            'sorted_idx_neg': [],
            'mz_values_pos': np.array([]),
            'mz_values_neg': np.array([]),
            'db_frag_pos': [],
            'db_frag_neg': [],
            'compound_info': {},
            'rt_values': {},
            'neutral_losses': {},
            'fragment_sources': {},
            'db_name': db_name,
            # 预构建的碎片来源快速查找表（按离子化模式）
            'frag_source_lookup_pos': {},  # frag_mz -> set of sources
            'frag_source_lookup_neg': {}   # frag_mz -> set of sources
        }
        
        if database.empty:
            return index_data
        
        # 创建本地索引到原始数据库的映射
        local_to_global = {}  # 本地索引 -> 原始数据库中的索引
        
        for local_idx, (orig_idx, row) in enumerate(database.iterrows()):
            local_to_global[local_idx] = orig_idx
            
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
            cas = str(row.get('CAS', ''))
            data_source = str(row.get('_data_source', db_name))
            is_standard = row.get('_is_standard', False)

            index_data['compound_info'][local_idx] = {
                'name_cn': cn_name,
                'name_en': en_name,
                'formula': formula,
                'herb': herb,
                'compound_type': compound_type,
                'source': source,
                'cas': cas if cas and cas != 'nan' else '',
                'adduct_pos': str(row.get('加合物（正）', '')),
                'adduct_neg': str(row.get('加合物（负）', '')),
                'data_source': data_source,
                'is_standard': is_standard,
                'db_name': db_name,
                'orig_idx': orig_idx  # 保留原始索引以区分不同数据库中的相同化合物
            }

            if '保留时间(min)' in row and pd.notna(row['保留时间(min)']):
                try:
                    index_data['rt_values'][local_idx] = float(row['保留时间(min)'])
                except:
                    pass

            if '中性丢失' in row and pd.notna(row['中性丢失']):
                losses = self._parse_losses(row['中性丢失'])
                if losses:
                    index_data['neutral_losses'][local_idx] = losses

            # 处理正离子数据
            if pd.notna(mz_pos) and str(mz_pos).strip() != '':
                try:
                    mz_val = float(mz_pos)
                    if mz_val > 0:
                        frag_string = row.get('碎片离子（正）', '')
                        fragments, source_map = parse_fragments_with_source(frag_string, source, db_name)
                        index_data['sorted_idx_pos'].append((mz_val, local_idx, fragments))
                        if fragments:
                            index_data['fragment_sources'][local_idx] = index_data['fragment_sources'].get(local_idx, {})
                            index_data['fragment_sources'][local_idx]['positive'] = source_map
                            # 预构建碎片来源查找表（加速匹配）
                            for frag_mz, src_set in source_map.items():
                                if frag_mz not in index_data['frag_source_lookup_pos']:
                                    index_data['frag_source_lookup_pos'][frag_mz] = set()
                                index_data['frag_source_lookup_pos'][frag_mz].update(src_set)
                except (ValueError, TypeError):
                    pass

            # 处理负离子数据
            if pd.notna(mz_neg) and str(mz_neg).strip() != '':
                try:
                    mz_val = float(mz_neg)
                    if mz_val > 0:
                        frag_string = row.get('碎片离子（负）', '')
                        fragments, source_map = parse_fragments_with_source(frag_string, source, db_name)
                        index_data['sorted_idx_neg'].append((mz_val, local_idx, fragments))
                        if fragments:
                            index_data['fragment_sources'][local_idx] = index_data['fragment_sources'].get(local_idx, {})
                            index_data['fragment_sources'][local_idx]['negative'] = source_map
                            # 预构建碎片来源查找表（加速匹配）
                            for frag_mz, src_set in source_map.items():
                                if frag_mz not in index_data['frag_source_lookup_neg']:
                                    index_data['frag_source_lookup_neg'][frag_mz] = set()
                                index_data['frag_source_lookup_neg'][frag_mz].update(src_set)
                except (ValueError, TypeError):
                    pass

        index_data['sorted_idx_pos'].sort(key=lambda x: x[0])
        index_data['sorted_idx_neg'].sort(key=lambda x: x[0])

        index_data['mz_values_pos'] = np.array([x[0] for x in index_data['sorted_idx_pos']])
        index_data['mz_values_neg'] = np.array([x[0] for x in index_data['sorted_idx_neg']])

        index_data['db_frag_pos'] = [x[2] for x in index_data['sorted_idx_pos']]
        index_data['db_frag_neg'] = [x[2] for x in index_data['sorted_idx_neg']]

        print(f"  {db_name}索引构建完成: {len(index_data['mz_values_pos'])} 条正离子, {len(index_data['mz_values_neg'])} 条负离子")
        
        return index_data

    def _build_global_fragment_sources(self):
        """构建全局碎片离子来源映射（包含所有数据库）"""
        self.global_fragment_sources = {}
        
        for idx_name, index_data in [('main', self.main_index), 
                                      ('english', self.english_index), 
                                      ('standard', self.standard_index)]:
            for local_idx, sources in index_data['fragment_sources'].items():
                compound_key = f"{idx_name}_{local_idx}"
                self.global_fragment_sources[compound_key] = sources

    def _search_database(self, precursor_mz, tolerance_ppm, ionization_mode, index_data):
        """在指定数据库索引中搜索候选化合物 - 增强匹配验证"""
        candidates = []
        
        if ionization_mode == 'positive' or ionization_mode == 'both':
            mz_values = index_data['mz_values_pos']
            db_frags = index_data['db_frag_pos']
            sorted_idx = index_data['sorted_idx_pos']
        else:
            mz_values = index_data['mz_values_neg']
            db_frags = index_data['db_frag_neg']
            sorted_idx = index_data['sorted_idx_neg']
        
        if len(mz_values) == 0:
            return candidates
        
        match_range = self._binary_search_range(mz_values, precursor_mz, tolerance_ppm)
        
        # 记录候选匹配信息用于验证
        for i in match_range:
            db_mz = mz_values[i]
            ppm_error = abs(precursor_mz - db_mz) / db_mz * 1e6 if db_mz > 0 else float('inf')
            
            # 获取化合物信息
            local_idx = sorted_idx[i][1] if i < len(sorted_idx) else i
            compound_info = index_data['compound_info'].get(local_idx, {})
            
            # 验证：确保化合物名称与准分子离子匹配
            # 检查化合物名称是否有效（不为空且不是无效值）
            name_cn = compound_info.get('name_cn', '')
            name_en = compound_info.get('name_en', '')
            
            # 如果化合物名称无效，跳过该候选
            if (not name_cn or name_cn == 'nan' or name_cn == '') and \
               (not name_en or name_en == 'nan' or name_en == ''):
                print(f"  警告: 跳过无效化合物名称 - m/z {db_mz}")
                continue
            
            candidate = {
                'db_idx': local_idx,
                'db_mz': db_mz,
                'ppm': ppm_error,
                'fragments': db_frags[i] if i < len(db_frags) else [],
                'compound_info': compound_info,
                'fragment_sources': index_data['fragment_sources'].get(local_idx, {}),
                'db_name': index_data['db_name'],
                'rt_values': index_data['rt_values'],
                'neutral_losses': index_data['neutral_losses'],
                'matching_verified': True  # 标记已验证
            }
            candidates.append(candidate)
        
        return candidates

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

    def _match_fragments_fast(self, observed, reference, tolerance_value, tolerance_type='Da', precursor_mz=None):
        """快速碎片匹配"""
        matched = []
        if len(reference) == 0:
            return matched
        ref_arr = np.asarray(reference)
        obs_arr = np.asarray(observed)
        for ref_val in ref_arr:
            if pd.notna(ref_val) and float(ref_val) > 0:
                for obs_val in obs_arr:
                    if precursor_mz is not None and abs(obs_val - precursor_mz) <= tolerance_value:
                        continue
                    if tolerance_type == 'Da':
                        if abs(obs_val - ref_val) <= tolerance_value:
                            matched.append(obs_val)
                            break
                    else:
                        ppm_error = abs(obs_val - ref_val) / ref_val * 1e6
                        if ppm_error <= tolerance_value:
                            matched.append(obs_val)
                            break
        return list(set(matched))

    def _match_fragments_with_source(self, observed, reference, tolerance_value, tolerance_type='Da', precursor_mz=None):
        """碎片匹配，返回匹配结果及来源"""
        matched = []
        matched_sources = {}
        
        if len(reference) == 0:
            return matched, matched_sources
            
        for obs_val in observed:
            if precursor_mz is not None and abs(obs_val - precursor_mz) <= tolerance_value:
                continue
                
            for ref_val in reference:
                if pd.notna(ref_val) and float(ref_val) > 0:
                    if tolerance_type == 'Da':
                        if abs(obs_val - ref_val) <= tolerance_value:
                            if obs_val not in matched:
                                matched.append(obs_val)
                            if obs_val not in matched_sources:
                                matched_sources[obs_val] = set()
                            # 添加该碎片的文献来源
                            if ref_val in self.temp_fragment_source_map:
                                matched_sources[obs_val].update(self.temp_fragment_source_map[ref_val])
                            break
                    else:
                        ppm_error = abs(obs_val - ref_val) / ref_val * 1e6
                        if ppm_error <= tolerance_value:
                            if obs_val not in matched:
                                matched.append(obs_val)
                            if obs_val not in matched_sources:
                                matched_sources[obs_val] = set()
                            if ref_val in self.temp_fragment_source_map:
                                matched_sources[obs_val].update(self.temp_fragment_source_map[ref_val])
                            break
        
        return matched, matched_sources

    def _find_diagnostic_ions_fast(self, matched_fragments, category, precursor_mz=None):
        """快速诊断离子匹配"""
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
        """提取母离子和碎片离子"""
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
                    'fragments': fragments,
                    'fragments_dict': fragments_dict,
                    'rt': rt if pd.notna(rt) else None,
                    'ionization_mode': ionization_mode
                })

        return precursors

    def identify_compound(self, precursor_mz, fragments, rt, ionization_mode):
        """鉴定单个化合物 - 增强验证确保化合物名称与母离子匹配"""
        tolerance_ppm = self.config['tolerance_ppm']
        fragment_tolerance = self.config['fragment_tolerance']
        self.temp_fragment_source_map = {}  # 临时存储碎片来源映射
        
        all_candidates = []
        
        # 在三个数据库中分别搜索
        databases = [
            (self.main_index, '主数据库'),
            (self.english_index, '英文数据库'),
            (self.standard_index, '对照品数据库')
        ]
        
        for index_data, db_name in databases:
            if index_data['mz_values_pos'].size == 0 and index_data['mz_values_neg'].size == 0:
                continue
            
            # 使用预构建的碎片来源查找表（加速）
            if ionization_mode in ['positive', 'both']:
                self.temp_fragment_source_map = index_data.get('frag_source_lookup_pos', {})
            else:
                self.temp_fragment_source_map = index_data.get('frag_source_lookup_neg', {})
            
            candidates = self._search_database(precursor_mz, tolerance_ppm, ionization_mode, index_data)
            all_candidates.extend(candidates)
        
        if not all_candidates:
            return []
        
        # 按ppm排序，保留所有ppm范围内的候选化合物
        all_candidates.sort(key=lambda x: x['ppm'])
        
        results = []
        for candidate in all_candidates:
            # 再次验证化合物名称有效
            compound_info = candidate['compound_info']
            name_cn = compound_info.get('name_cn', '')
            name_en = compound_info.get('name_en', '')
            
            # 跳过无效化合物名称的候选
            if (not name_cn or name_cn == 'nan' or name_cn == '') and \
               (not name_en or name_en == 'nan' or name_en == ''):
                continue
            
            matched_frags, matched_frag_sources = self._match_fragments_with_source(
                fragments,
                candidate['fragments'],
                fragment_tolerance,
                self.tolerance_type,
                precursor_mz
            )
            
            # 优化：碎片匹配过滤规则
            # 1. 如果数据库中有碎片离子数据（不为空），但全部匹配失败，则跳过
            # 2. 如果数据库中没有碎片离子数据（为空），则保留（允许只有母离子匹配）
            # 3. 如果有匹配成功的碎片，则保留
            has_db_fragments = len(candidate['fragments']) > 0
            has_matched_frags = len(matched_frags) > 0
            
            # 只有有碎片数据且全部匹配失败时才跳过
            if has_db_fragments and not has_matched_frags:
                continue
            
            category = self._classify_compound(
                compound_info.get('name_cn', '') + ' ' + compound_info.get('name_en', ''),
                compound_info.get('compound_type', '')
            )
            
            diagnostic, diag_weights = self._find_diagnostic_ions_fast(matched_frags, category, precursor_mz)
            
            rating, rating_name, confidence, recommendation = self._determine_confidence_level(
                candidate['ppm'],
                len(matched_frags),
                len(diagnostic),
                len(candidate['fragments']) > 0
            )
            
            # 获取保留时间偏差
            rt_deviation = None
            if rt is not None and candidate['db_idx'] in candidate['rt_values']:
                db_rt = candidate['rt_values'][candidate['db_idx']]
                rt_deviation = abs(rt - db_rt) if db_rt else None
            
            base_score = self._calculate_base_score(
                rating,
                candidate['ppm'],
                len(matched_frags),
                len(diagnostic),
                diag_weights,
                0,
                rt_deviation
            )
            
            # 格式化碎片离子来源
            frag_sources_formatted = []
            for frag_mz, src_set in matched_frag_sources.items():
                sources_str = '; '.join(sorted(src_set))
                frag_sources_formatted.append(f"{frag_mz}({sources_str})")
            
            # 构建碎片离子列表（不合并，每个碎片保留来源）
            fragment_list = []
            for frag_mz, src_set in matched_frag_sources.items():
                sources_str = '; '.join(sorted(src_set))
                fragment_list.append({
                    'fragment_mz': frag_mz,
                    'sources': sources_str,
                    'display': f"{frag_mz}({sources_str})"
                })
            
            result = {
                '母离子m/z': precursor_mz,
                '观测RT': rt,
                'ppm': round(candidate['ppm'], 4),
                'db_mz': candidate['db_mz'],
                '化合物中文名': name_cn if name_cn != 'nan' else '',
                '化合物英文名': name_en if name_en != 'nan' else '',
                '分子式': compound_info.get('formula', ''),
                'CAS号': compound_info.get('cas', ''),
                '药材名称': compound_info.get('herb', ''),
                '化合物类型': compound_info.get('compound_type', ''),
                '数据来源': candidate['db_name'],  # 明确标注数据库来源
                '是否为对照品': '是' if compound_info.get('is_standard', False) else '否',
                '离子化方式': ionization_mode,
                '加和离子': compound_info.get('adduct_pos', '') if ionization_mode in ['positive', 'both'] else compound_info.get('adduct_neg', ''),
                '_fragment_list': fragment_list,  # 内部存储，用于后续处理
                '文献来源': compound_info.get('source', ''),
                '诊断性离子': '; '.join([str(d) for d in diagnostic]) if diagnostic else '',
                '匹配碎片数': len(matched_frags),
                '诊断离子数': len(diagnostic),
                '评级': rating,
                '评级名称': rating_name,
                '置信度': confidence,
                '报告建议': recommendation,
                '综合得分': base_score,
                '基础得分': base_score
            }
            
            results.append(result)
        
        return results

    def generate_report(self, sample_name='样品'):
        """生成鉴定报告"""
        records = []
        
        print("\n【8/9】正在处理正离子数据...")
        positive_precursors = self.extract_precursor_ions(self.ms_positive, 'positive')
        self.stats['total_precursors'] += len(positive_precursors)
        print(f"  正离子模式: {len(positive_precursors)} 个母离子待鉴定")
        
        for i, precursor in enumerate(positive_precursors):
            if (i + 1) % 10 == 0:
                print(f"  正离子进度: {i + 1}/{len(positive_precursors)}")
            
            results = self.identify_compound(
                precursor['precursor_mz'],
                precursor['fragments'],
                precursor['rt'],
                precursor['ionization_mode']
            )
            
            for result in results:
                result['出峰时间t/min'] = precursor['rt']
                records.append(result)
        
        print("\n【8/9】正在处理负离子数据...")
        negative_precursors = self.extract_precursor_ions(self.ms_negative, 'negative')
        self.stats['total_precursors'] += len(negative_precursors)
        print(f"  负离子模式: {len(negative_precursors)} 个母离子待鉴定")
        
        for i, precursor in enumerate(negative_precursors):
            if (i + 1) % 10 == 0:
                print(f"  负离子进度: {i + 1}/{len(negative_precursors)}")
            
            results = self.identify_compound(
                precursor['precursor_mz'],
                precursor['fragments'],
                precursor['rt'],
                precursor['ionization_mode']
            )
            
            for result in results:
                result['出峰时间t/min'] = precursor['rt']
                records.append(result)
        
        print(f"\n【9/9】初步匹配完成，共 {len(records)} 条候选记录，正在合并...")
        
        report_df = self._merge_and_fuse_records(records)
        
        self.stats['identified_compounds'] = len(report_df)
        
        print("\n" + "="*80)
        print(f"鉴定完成！共识别出 {len(report_df)} 个化合物")
        print("="*80)
        
        return report_df

    def _merge_and_fuse_records(self, records_list):
        """合并和融合记录"""
        if not records_list:
            return pd.DataFrame()
        
        # 创建合并键
        for record in records_list:
            merge_key = f"{record.get('CAS号', '')}_{record.get('分子式', '')}_{record.get('化合物中文名', '')}"
            record['_merge_key'] = merge_key
            # 保留碎片离子列表（不合并）
            record['_fragment_list'] = record.get('_fragment_list', [])
        
        best_records = {}
        herbs_collection = defaultdict(set)
        
        for record in records_list:
            merge_key = record['_merge_key']
            base_score = record.get('综合得分', 0)
            
            herbs_collection[merge_key].add(record['药材名称'])
            
            current_score = base_score
            if merge_key not in best_records or current_score > best_records[merge_key][1]:
                best_records[merge_key] = (record, current_score)
        
        fused_records = []
        for merge_key, (record, score) in best_records.items():
            herbs = herbs_collection[merge_key]
            record['药材名称'] = '; '.join(sorted(herbs))
            if '_merge_key' in record:
                del record['_merge_key']
            if '基础得分' in record:
                del record['基础得分']
            fused_records.append(record)
        
        if not fused_records:
            return pd.DataFrame()
        
        fused_records.sort(key=lambda x: x['出峰时间t/min'] if x['出峰时间t/min'] is not None else 0)
        
        rt_tol = self.rt_fusion_tolerance
        final_records = []
        i = 0
        n = len(fused_records)
        
        while i < n:
            group = [fused_records[i]]
            j = i + 1
            while j < n:
                rt1 = group[0]['出峰时间t/min']
                rt2 = fused_records[j]['出峰时间t/min']
                if rt1 is None or rt2 is None:
                    break
                if abs(rt1 - rt2) <= rt_tol:
                    group.append(fused_records[j])
                    j += 1
                else:
                    break
            
            compound_groups = {}
            for rec in group:
                key = rec.get('CAS号', rec.get('化合物英文名', ''))
                if key not in compound_groups:
                    compound_groups[key] = []
                compound_groups[key].append(rec)
            
            for comp_key, comp_recs in compound_groups.items():
                if len(comp_recs) == 1:
                    # 单条记录，格式化碎片离子输出
                    rec = comp_recs[0]
                    self._format_fragment_output(rec)
                    final_records.append(rec)
                else:
                    best_rec = max(comp_recs, key=lambda x: x['综合得分']).copy()
                    
                    # 合并数据来源
                    data_sources = set()
                    for rec in comp_recs:
                        if rec.get('数据来源'):
                            data_sources.add(rec['数据来源'])
                    best_rec['数据来源'] = '; '.join(sorted(data_sources))
                    
                    # 合并其他字段
                    adducts = set()
                    modes = set()
                    for rec in comp_recs:
                        if rec.get('加和离子'):
                            adducts.add(rec['加和离子'])
                        if rec.get('离子化方式'):
                            modes.add(rec['离子化方式'])
                    best_rec['加和离子'] = '; '.join(sorted(adducts))
                    best_rec['离子化方式'] = '/'.join(sorted(modes))
                    
                    # 合并碎片离子（保留每个碎片的独立来源，不合并）
                    all_fragments = {}
                    for rec in comp_recs:
                        frag_list = rec.get('_fragment_list', [])
                        for frag in frag_list:
                            frag_mz = frag['fragment_mz']
                            if frag_mz not in all_fragments:
                                all_fragments[frag_mz] = {
                                    'fragment_mz': frag_mz,
                                    'sources': set(),
                                    'display_parts': []
                                }
                            # 合并来源（去重）
                            for src in frag['sources'].split('; '):
                                if src:
                                    all_fragments[frag_mz]['sources'].add(src.strip())
                    
                    # 构建碎片离子输出（不合并，每个碎片独立）
                    fragment_output = []
                    for frag_mz, frag_data in sorted(all_fragments.items()):
                        sources_str = '; '.join(sorted(frag_data['sources']))
                        fragment_output.append(f"{frag_mz}({sources_str})")
                    
                    best_rec['主要碎片离子'] = '; '.join(fragment_output) if fragment_output else ''
                    best_rec['匹配碎片数'] = len(all_fragments)
                    
                    # 合并文献来源
                    all_sources = set()
                    for rec in comp_recs:
                        if rec.get('文献来源'):
                            for s in str(rec['文献来源']).split(';'):
                                all_sources.add(s.strip())
                    best_rec['文献来源'] = '; '.join(sorted(all_sources))
                    best_rec['文献来源数'] = len(all_sources)
                    
                    # 合并诊断性离子
                    all_diag = set()
                    for rec in comp_recs:
                        if rec.get('诊断性离子'):
                            for d in str(rec['诊断性离子']).split('; '):
                                all_diag.add(d.strip())
                    best_rec['诊断性离子'] = '; '.join(sorted(all_diag)) if all_diag else ''
                    
                    # 删除内部字段
                    if '_fragment_list' in best_rec:
                        del best_rec['_fragment_list']
                    
                    # 评分调整
                    extra_count = len(data_sources) - 1
                    if extra_count > 0:
                        fusion_bonus = min(extra_count * 5, 15)
                    else:
                        fusion_bonus = 0
                    
                    final_score = min(best_rec['综合得分'] + fusion_bonus, 100)
                    
                    if best_rec['评级'] == 1 and final_score < 80:
                        final_score = 80
                    if best_rec['评级'] == 2 and final_score < 60:
                        final_score = 60
                    
                    best_rec['综合得分'] = final_score
                    
                    if len(comp_recs) >= 2 and all(rec.get('评级', 5) <= 2 for rec in comp_recs):
                        best_rec['评级'] = 1
                        best_rec['评级名称'] = '确证级'
                        best_rec['置信度'] = '最高'
                        best_rec['报告建议'] = '可直接报告'
                    
                    final_records.append(best_rec)
            
            i = j
        
        report_df = pd.DataFrame(final_records)
        
        if not report_df.empty:
            report_df = report_df.sort_values(by=['评级', '综合得分', 'ppm'], ascending=[True, False, True])
            report_df = report_df.reset_index(drop=True)
            report_df['序号'] = range(1, len(report_df) + 1)
        
        return report_df

    def save_report(self, report_df, output_path):
        """保存报告"""
        report_df.to_excel(output_path, index=False)
        print(f"\n报告已保存至: {output_path}")

    def print_summary(self, report_df, herb_name):
        """打印摘要"""
        print("\n" + "="*100)
        print(f"【{herb_name}化合物鉴定报告】")
        print("="*100)
        print(f"鉴定药材: {herb_name}")
        print(f"鉴定时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"共鉴定出 {len(report_df)} 个化合物")

        print(f"\n【处理统计】")
        print(f"  - 总母离子数: {self.stats['total_precursors']}")
        print(f"  - 鉴定化合物数: {self.stats['identified_compounds']}")

        if not report_df.empty:
            # 数据来源统计
            if '数据来源' in report_df.columns:
                print(f"\n【数据来源分布】")
                for source, count in report_df['数据来源'].value_counts().items():
                    print(f"  - {source}: {count} 个")
            
            # 对照品统计
            if '是否为对照品' in report_df.columns:
                std_count = (report_df['是否为对照品'] == '是').sum()
                print(f"\n【对照品匹配】")
                print(f"  - 匹配到对照品: {std_count} 个")

            print(f"\n【评级分布】")
            for level, count in report_df['评级名称'].value_counts().items():
                print(f"  - {level}: {count} 个")

        print("\n" + "="*100)

    def _format_fragment_output(self, record):
        """格式化碎片离子输出（保持独立，不合并）"""
        fragment_list = record.get('_fragment_list', [])
        if fragment_list:
            # 每个碎片离子单独显示，保留来源
            fragment_output = []
            for frag in fragment_list:
                display = frag.get('display', f"{frag.get('fragment_mz')}({frag.get('sources', '')})")
                fragment_output.append(display)
            record['主要碎片离子'] = '; '.join(fragment_output)
        else:
            record['主要碎片离子'] = ''
        
        # 删除内部字段
        if '_fragment_list' in record:
            del record['_fragment_list']
        
        return record

    def _print_initialization_info(self):
        """打印初始化信息"""
        print("\n" + "="*80)
        print("程序初始化完成（三库独立比对版 v6.1 - 修复母离子匹配）")
        print("="*80)
        print(f"  - 主数据库索引: {len(self.main_index['mz_values_pos'])} 条正离子, {len(self.main_index['mz_values_neg'])} 条负离子")
        print(f"  - 英文数据库索引: {len(self.english_index['mz_values_pos'])} 条正离子, {len(self.english_index['mz_values_neg'])} 条负离子")
        print(f"  - 对照品数据库索引: {len(self.standard_index['mz_values_pos'])} 条正离子, {len(self.standard_index['mz_values_neg'])} 条负离子")
        print(f"  - 诊断性离子库: {len(self.diagnostic_ions)} 类")
        print(f"  - 容差类型: {self.tolerance_type}")
        print(f"  - 强度相对阈值: {self.intensity_relative_threshold*100:.1f}%")
        print(f"  - RT得分: {'启用' if self.use_rt_score else '禁用'}")
        print(f"  - 并行处理: {'启用' if self.use_parallel else '禁用'}")
        print("="*80)


from collections import defaultdict


# ============================================================================
# 诊断离子匹配函数
# ============================================================================

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
# Streamlit 网页应用部分（保持不变）
# ============================================================================

st.set_page_config(
    page_title="中药化合物智能鉴定平台 v6.1",
    page_icon="data:image/svg+xml,<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 100 100'><text y='.9em' font-size='90'>🌿</text></svg>",
    layout="wide",
    initial_sidebar_state="expanded"
)


def load_optimized_css():
    """加载优化后的CSS样式"""
    st.markdown("""
    <style>
        @import url('https://fonts.googleapis.com/css2?family=Noto+Sans+SC:wght@300;400;500;700&display=swap');
        * {
            font-family: 'Noto Sans SC', 'Microsoft YaHei', sans-serif !important;
        }
        .stApp {
            background: linear-gradient(135deg, #f8fafc 0%, #e2e8f0 50%, #f1f5f9 100%);
            min-height: 100vh;
        }
        .main-header {
            background: linear-gradient(135deg, #059669 0%, #0891b2 50%, #7c3aed 100%);
            padding: 2rem 2rem;
            border-radius: 20px;
            margin-bottom: 2rem;
            color: white;
            box-shadow: 0 20px 40px rgba(5, 150, 105, 0.3);
        }
        .main-header h1 {
            color: white !important;
            font-size: 2rem;
            font-weight: 700;
            margin: 0;
        }
        .stat-card {
            background: rgba(255, 255, 255, 0.95);
            border-radius: 16px;
            padding: 1rem;
            text-align: center;
            box-shadow: 0 4px 20px rgba(0, 0, 0, 0.08);
        }
        .stat-number {
            font-size: 2rem;
            font-weight: 700;
            background: linear-gradient(135deg, #059669, #0891b2);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
        }
        .feature-card {
            background: white;
            border-radius: 16px;
            padding: 1.5rem;
            box-shadow: 0 4px 20px rgba(0, 0, 0, 0.08);
            margin-bottom: 1rem;
        }
    </style>
    """, unsafe_allow_html=True)


def login_page():
    """显示登录页面"""
    st.markdown("""
    <div style="max-width: 450px; margin: 80px auto; padding: 2rem; background: white; border-radius: 24px; box-shadow: 0 20px 60px rgba(0, 0, 0, 0.15);">
        <div style="text-align: center; font-size: 1.5rem; font-weight: 700; margin-bottom: 1.5rem; background: linear-gradient(135deg, #059669, #0891b2); -webkit-background-clip: text; -webkit-text-fill-color: transparent;">🌿 中药化合物智能鉴定平台</div>
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
                st.success("登录成功！正在跳转...")
                time.sleep(0.5)
                st.rerun()
            else:
                st.error("用户名或密码错误，请重试。")


def logout_button():
    """显示登出按钮"""
    if st.sidebar.button("登出", use_container_width=True):
        st.session_state.logged_in = False
        st.session_state.pop('username', None)
        st.rerun()


def create_header():
    """创建应用头部"""
    username = st.session_state.get('username', 'Guest')
    st.markdown(f"""
    <div class="main-header">
        <h1>🌿 中药化合物智能鉴定平台</h1>
        <p>v6.1 三库独立比对版（修复母离子匹配） | 欢迎回来，{username}</p>
    </div>
    """, unsafe_allow_html=True)


def create_sidebar():
    """创建侧边栏"""
    st.sidebar.markdown("""
    <div style="text-align: center; padding: 1rem 0; border-bottom: 2px solid #e2e8f0;">
        <h3 style="color: #1e293b; margin: 0;">TCM Identifier</h3>
        <p style="color: #64748b; font-size: 0.8rem;">中药化合物鉴定系统</p>
    </div>
    """, unsafe_allow_html=True)

    if st.session_state.get('logged_in'):
        st.sidebar.markdown(f"""
        <div style="background: linear-gradient(135deg, #ecfdf5, #d1fae5); padding: 0.75rem; border-radius: 12px; margin: 1rem 0; text-align: center;">
            <p style="margin: 0; color: #059669; font-weight: 600;">👤 {st.session_state.username}</p>
        </div>
        """, unsafe_allow_html=True)

    page = st.sidebar.radio(
        "选择功能",
        ["首页", "开始鉴定", "诊断离子筛查", "使用指南", "数据库预览", "结果分析"],
        index=0
    )

    st.sidebar.markdown("---")
    st.sidebar.markdown("""
    <div style="background: white; border-radius: 12px; padding: 0.75rem;">
        <h4 style="color: #1e293b; margin: 0 0 0.5rem 0; font-size: 0.9rem;">版本信息</h4>
        <ul style="color: #64748b; font-size: 0.8rem; padding-left: 1rem;">
            <li>v6.1 三库独立比对版</li>
            <li>三个数据库独立比对</li>
            <li>碎片离子标注来源</li>
            <li>修复母离子匹配问题</li>
        </ul>
    </div>
    """, unsafe_allow_html=True)

    return page


def show_home_page():
    """首页"""
    create_header()

    st.markdown("## 系统概览")

    cols = st.columns(4)
    stats_data = [
        ("35,828+", "数据库化合物数", "📚"),
        ("300+", "支持药材种类", "🌱"),
        ("6", "鉴定评级级别", "⭐"),
        ("3", "独立数据库", "💾")
    ]

    for col, (number, label, icon) in zip(cols, stats_data):
        with col:
            st.markdown(f"""
            <div class="stat-card">
                <div style="font-size: 1.5rem;">{icon}</div>
                <div class="stat-number">{number}</div>
                <div style="color: #64748b; margin-top: 0.25rem;">{label}</div>
            </div>
            """, unsafe_allow_html=True)

    st.markdown("---")
    st.markdown("## 核心功能")

    features = [
        ("🔍", "三库独立比对", "主数据库、英文数据库、对照品数据库分别独立比对"),
        ("📈", "六级评级标准", "确证级、高置信级、推定级、提示级、参考级、排除级"),
        ("🧪", "碎片离子溯源", "碎片离子标注来源于哪个数据库和文献"),
        ("📁", "灵活数据上传", "可单独上传正离子或负离子数据文件"),
        ("⚡", "高效索引加速", "使用索引缓存加速大规模数据处理"),
        ("🔄", "智能结果融合", "自动合并同化合物不同离子模式的结果"),
        ("✅", "母离子验证", "确保化合物名称与准分子离子匹配")
    ]

    col1, col2 = st.columns(2)
    for i, (icon, title, desc) in enumerate(features):
        with col1 if i % 2 == 0 else col2:
            st.markdown(f"""
            <div class="feature-card">
                <h4>{icon} {title}</h4>
                <p>{desc}</p>
            </div>
            """, unsafe_allow_html=True)

    if st.button("立即开始鉴定", type="primary", use_container_width=True):
        st.session_state['page'] = '开始鉴定'
        st.rerun()


def show_analysis_page():
    """鉴定分析页面"""
    create_header()

    st.markdown("## 上传质谱数据")
    st.markdown("**至少上传一个文件（正离子或负离子）**")

    col1, col2 = st.columns(2)
    with col1:
        ms_positive_file = st.file_uploader(
            "上传正离子模式质谱数据 (.xlsx，可选)",
            type=['xlsx'],
            key='ms_positive'
        )
        if ms_positive_file:
            st.success(f"已上传: {ms_positive_file.name}")
    with col2:
        ms_negative_file = st.file_uploader(
            "上传负离子模式质谱数据 (.xlsx，可选)",
            type=['xlsx'],
            key='ms_negative'
        )
        if ms_negative_file:
            st.success(f"已上传: {ms_negative_file.name}")

    if not ms_positive_file and not ms_negative_file:
        st.warning("请至少上传一个质谱数据文件")

    st.markdown("---")
    st.markdown("## 数据库配置")

    col1, col2 = st.columns(2)
    with col1:
        use_english_db = st.checkbox("启用英文数据库（数据库（英文）.xlsx）", value=True)
        use_standard_db = st.checkbox("启用对照品数据库（对照品数据库.xlsx）", value=True)
    with col2:
        custom_db_file = st.file_uploader(
            "上传自定义数据库 (.xlsx，可选)",
            type=['xlsx'],
            key='custom_db'
        )
        if custom_db_file:
            st.success(f"已上传自定义数据库: {custom_db_file.name}")

    st.markdown("---")
    st.markdown("## 鉴定参数配置")

    preset = st.radio("参数预设", options=["自定义", "快速模式", "高精度模式"], horizontal=True)

    col1, col2, col3 = st.columns(3)
    with col1:
        tolerance_ppm = st.number_input("ppm误差容限", min_value=10, max_value=100, value=50)
    with col2:
        max_candidates = st.number_input("最大候选数", min_value=1, max_value=10, value=3)
    with col3:
        min_intensity_abs = st.number_input("最小绝对强度", min_value=0, value=100)

    col1, col2, col3 = st.columns(3)
    with col1:
        tolerance_type = st.selectbox("碎片匹配容差类型", options=["Da", "ppm"], index=0)
        if tolerance_type == "Da":
            fragment_tolerance = st.number_input("碎片匹配容差 (Da)", min_value=0.01, max_value=1.0, value=0.05, step=0.01)
            fragment_tolerance_ppm = 20
        else:
            fragment_tolerance_ppm = st.number_input("碎片匹配容差 (ppm)", min_value=1, max_value=100, value=20)
            fragment_tolerance = 0.05
    with col2:
        rt_tolerance = st.number_input("保留时间容差 (min)", min_value=0.1, max_value=2.0, value=0.3, step=0.1)
    with col3:
        use_rt_score = st.checkbox("启用保留时间得分", value=True)

    if preset == "快速模式":
        tolerance_ppm = 50
        fragment_tolerance = 0.05
        min_intensity_abs = 100
        tolerance_type = "Da"
    elif preset == "高精度模式":
        tolerance_ppm = 10
        fragment_tolerance = 0.02
        min_intensity_abs = 50
        tolerance_type = "ppm"

    st.markdown("---")

    if ms_positive_file or ms_negative_file:
        if st.button("开始化合物鉴定", type="primary", use_container_width=True):
            with st.spinner("正在初始化鉴定程序..."):
                try:
                    temp_dir = tempfile.gettempdir()
                    pos_path = neg_path = custom_db_path = english_db_path = standard_db_path = None

                    if ms_positive_file:
                        pos_path = os.path.join(temp_dir, ms_positive_file.name)
                        with open(pos_path, 'wb') as f:
                            f.write(ms_positive_file.getbuffer())
                    if ms_negative_file:
                        neg_path = os.path.join(temp_dir, ms_negative_file.name)
                        with open(neg_path, 'wb') as f:
                            f.write(ms_negative_file.getbuffer())
                    if custom_db_file:
                        custom_db_path = os.path.join(temp_dir, custom_db_file.name)
                        with open(custom_db_path, 'wb') as f:
                            f.write(custom_db_file.getbuffer())
                    
                    if use_english_db:
                        english_db_path = find_english_database_path()
                        if not english_db_path:
                            st.warning("未找到英文数据库文件")
                    
                    if use_standard_db:
                        standard_db_path = find_standard_database_path()
                        if standard_db_path:
                            st.info(f"已找到对照品数据库: {standard_db_path}")
                        else:
                            st.warning("未找到对照品数据库文件（对照品数据库.xlsx）")

                    db_path = find_database_path()
                    if not db_path:
                        st.error("未找到主数据库文件！请将 TCM-SM-MS DB.xlsx 放在项目目录下。")
                        return

                    st.info(f"已找到主数据库文件: {db_path}")

                    progress_bar = st.progress(0)
                    status_text = st.empty()

                    status_text.text("正在加载数据库...")
                    progress_bar.progress(10)

                    config = {
                        'min_intensity': min_intensity_abs,
                        'fragment_tolerance': fragment_tolerance,
                        'fragment_tolerance_ppm': fragment_tolerance_ppm,
                        'tolerance_ppm': tolerance_ppm,
                        'max_candidates': max_candidates,
                    }

                    identifier = UltimateGardeniaIdentifier(
                        database_path=db_path,
                        ms_positive_path=pos_path,
                        ms_negative_path=neg_path,
                        herb_name=None,
                        config=config,
                        use_parallel=True,
                        rt_tolerance=rt_tolerance,
                        loss_tolerance=0.02,
                        tolerance_type=tolerance_type,
                        use_rt_score=use_rt_score,
                        custom_db_path=custom_db_path,
                        english_db_path=english_db_path,
                        standard_db_path=standard_db_path,
                        cache_index=True
                    )

                    progress_bar.progress(50)
                    status_text.text("正在处理质谱数据...")

                    report = identifier.generate_report('样品')

                    progress_bar.progress(90)
                    status_text.text("生成报告...")

                    st.session_state['analysis_results'] = report
                    st.session_state['identifier'] = identifier

                    progress_bar.progress(100)
                    status_text.text("鉴定完成！")

                    st.success("化合物鉴定完成！")

                    if st.button("查看鉴定结果"):
                        st.session_state['page'] = '结果分析'
                        st.rerun()

                except Exception as e:
                    st.error(f"鉴定过程中出错：{str(e)}")
                    st.exception(e)


def show_diagnostic_ion_page():
    """诊断离子筛查页面"""
    create_header()

    st.markdown("## 诊断离子筛查")

    diagnostic_df = load_diagnostic_ions_cached()

    if diagnostic_df.empty:
        st.warning("未找到诊断离子数据库文件（诊断离子.xlsx）")
        return

    st.success(f"已加载诊断离子数据库，包含 {len(diagnostic_df)} 条记录")

    mz_input = st.text_area("输入m/z值（每行一个值，或用逗号分隔）", placeholder="例如：\n151.003\n137.024, 121.029", height=150)

    col1, col2 = st.columns([3, 1])
    with col2:
        tolerance_ppm = st.number_input("ppm容差", min_value=1, max_value=100, value=10)
        ion_mode = st.selectbox("离子模式", options=["全部", "正离子", "负离子"])

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
            st.warning("无法解析输入的m/z值")
        else:
            with st.spinner("正在匹配诊断离子..."):
                results_df = match_diagnostic_ions(user_mz_values, diagnostic_df, tolerance_ppm=tolerance_ppm, ion_mode=ion_mode)

            if not results_df.empty:
                st.markdown("### 匹配结果统计")
                cols = st.columns(3)
                with cols[0]:
                    st.metric("输入离子数", len(user_mz_values))
                with cols[1]:
                    st.metric("匹配离子数", results_df['输入m/z'].nunique())
                with cols[2]:
                    st.metric("化合物类型数", results_df['化合物类型'].nunique())

                st.markdown("### 化合物类型分布")
                st.bar_chart(results_df['化合物类型'].value_counts())

                st.markdown("### 匹配结果详情")
                display_cols = ['输入m/z', '匹配诊断离子m/z', '误差(ppm)', '化合物类型', '中文名称', '英文名称']
                available_cols = [c for c in display_cols if c in results_df.columns]
                st.dataframe(results_df[available_cols], use_container_width=True, hide_index=True)
            else:
                st.info("未找到匹配的诊断离子")


def show_guide_page():
    """使用指南页面"""
    create_header()
    st.markdown("## 使用指南（v6.1 三库独立比对版 - 修复母离子匹配）")

    st.markdown("""
    ### 三库独立比对机制

    本版本实现了三个数据库的独立比对：

    1. **主数据库 (TCM-SM-MS DB.xlsx)**：中药小分子化学成分高分辨质谱数据库
    2. **英文数据库 (数据库（英文）.xlsx)**：英文化合物数据库
    3. **对照品数据库 (对照品数据库.xlsx)**：标准品参考数据库

    ### 碎片离子溯源

    每个碎片离子都会标注其来源：
    - `[主数据库] 文献来源`
    - `[英文数据库] 文献来源`
    - `[对照品数据库]`

    ### 母离子匹配验证（v6.1新增）

    - 自动验证化合物名称是否有效
    - 跳过化合物名称为空的无效记录
    - 确保数据库中的化合物名称与准分子离子对应

    ### 评级标准

    | 等级 | ppm要求 | 碎片要求 | 说明 |
    |------|---------|----------|------|
    | 确证级 | ≤10ppm | ≥2个碎片+≥1个诊断离子 | 可直接报告 |
    | 高置信级 | ≤20ppm | ≥2个碎片 | 建议复核后报告 |
    | 推定级 | ≤50ppm | ≥1个碎片 | 需验证后报告 |
    """)


def show_database_page():
    """数据库预览页面"""
    create_header()
    st.markdown("## 数据库预览")

    # 主数据库
    db_path = find_database_path()
    if db_path:
        try:
            df = load_database_cached()
            if not df.empty:
                st.success(f"主数据库（TCM-SM-MS DB）：{len(df)} 条记录")
                st.info(f"路径: {db_path}")
                st.dataframe(df.head(5), use_container_width=True)
        except Exception as e:
            st.error(f"加载主数据库时出错：{str(e)}")
    else:
        st.warning("未找到主数据库文件！")

    # 英文数据库
    st.markdown("---")
    st.markdown("### 英文数据库（数据库（英文）.xlsx）")
    eng_path = find_english_database_path()
    if eng_path:
        try:
            eng_df = load_english_database_cached()
            if not eng_df.empty:
                st.success(f"英文数据库：{len(eng_df)} 条记录")
                display_cols = ['名称（中文）', '名称（英文）', '分子式', '准分子离子（正）']
                available_cols = [c for c in display_cols if c in eng_df.columns]
                st.dataframe(eng_df[available_cols].head(5), use_container_width=True)
        except Exception as e:
            st.error(f"加载英文数据库时出错：{str(e)}")
    else:
        st.info("未找到英文数据库文件")

    # 对照品数据库
    st.markdown("---")
    st.markdown("### 对照品数据库（对照品数据库.xlsx）")
    std_path = find_standard_database_path()
    if std_path:
        try:
            std_df = load_standard_database_cached()
            if not std_df.empty:
                st.success(f"对照品数据库：{len(std_df)} 条记录")
                st.dataframe(std_df.head(5), use_container_width=True)
        except Exception as e:
            st.error(f"加载对照品数据库时出错：{str(e)}")
    else:
        st.info("未找到对照品数据库文件")


def show_results_page():
    """结果分析页面"""
    create_header()

    if 'analysis_results' not in st.session_state:
        st.markdown("""
        <div style="text-align: center; padding: 2rem; color: #94a3b8;">
            <h3>暂无鉴定结果</h3>
            <p>请先进行化合物鉴定</p>
        </div>
        """, unsafe_allow_html=True)
        if st.button("前往鉴定页面"):
            st.session_state['page'] = '开始鉴定'
            st.rerun()
        return

    report = st.session_state['analysis_results']

    st.markdown("## 鉴定结果分析（v6.1 三库独立比对版 - 修复母离子匹配）")

    if report.empty:
        st.warning("鉴定结果为空")
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

    # 数据来源分布
    if '数据来源' in report.columns:
        st.markdown("### 数据来源分布")
        st.bar_chart(report['数据来源'].value_counts())

    st.markdown("### 评级分布")
    level_counts = report['评级名称'].value_counts()
    st.bar_chart(level_counts.reindex(['确证级', '高置信级', '推定级', '提示级', '参考级']).fillna(0))

    st.markdown("### 完整鉴定结果")

    all_columns = report.columns.tolist()
    default_cols = ['序号', '化合物中文名', '分子式', 'ppm', '评级名称', '药材名称', '综合得分', '数据来源']
    selected_cols = st.multiselect("选择显示的列", all_columns, default=[c for c in default_cols if c in all_columns])

    display_df = report[selected_cols] if selected_cols else report
    st.dataframe(display_df, use_container_width=True, hide_index=True)

    st.markdown("---")
    st.markdown("### 导出报告")

    col1, col2 = st.columns(2)
    with col1:
        csv = report.to_csv(index=False, encoding='utf-8-sig')
        st.download_button(label="导出CSV", data=csv, file_name=f"鉴定报告_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv", mime="text/csv")
    with col2:
        buffer = BytesIO()
        with pd.ExcelWriter(buffer, engine='openpyxl') as writer:
            report.to_excel(writer, index=False, sheet_name='鉴定结果')
        st.download_button(label="导出Excel", data=buffer.getvalue(), file_name=f"鉴定报告_{datetime.now().strftime('%Y%m%d_%H%M%S')}.xlsx", mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")


def main():
    """主函数"""
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
