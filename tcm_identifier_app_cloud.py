# -*- coding: utf-8 -*-
"""
中药化合物智能鉴定平台 - 云端部署增强版 v5.5（评分优化版）
==========================================================

功能特点：
- 自动加载项目目录下的诊断离子.xlsx（如果存在）
- 支持上传自定义诊断离子文件覆盖默认
- 重新设计的综合评分系统：确保1、2级化合物得分>60，级别越小得分越高
- 择优输出最佳候选，自动区分同分异构体
- 六级评级标准 + 综合得分辅助判断
- 智能去重基于得分，并合并所有药材来源
- 分子式标准化（去除Unicode下标）
- 加和离子字段清洗（多值时只取第一个）

作者：张永
日期：2026-03-12
版本：v5.5（评分优化版）
"""

import streamlit as st
import pandas as pd
import numpy as np
import os
import time
from datetime import datetime
from io import BytesIO
import tempfile
from bisect import bisect_left, bisect_right
import warnings
warnings.filterwarnings('ignore')


# ============================================================================
# 鉴定程序核心代码（评分优化版）
# ============================================================================

class UltimateGardeniaIdentifier:
    """
    中药化合物鉴定终极版程序 v5.5（评分优化版）
    
    新增特性：
    1. 支持外部诊断离子文件（Excel格式），自动查找项目目录下的诊断离子.xlsx
    2. 重新设计的综合评分系统：基于评级、ppm、匹配碎片数、诊断离子个数计算，确保1、2级得分>60
    3. 保留时间匹配已移除（根据用户要求）
    4. 中性丢失暂不计分（可扩展）
    5. 择优输出：每个母离子仅返回得分最高的候选
    6. 去重基于得分，并合并所有药材来源
    7. 分子式标准化（去除Unicode下标）
    8. 加和离子字段清洗（多值时只取第一个）
    """
    
    def __init__(self, database_path, ms_positive_path, ms_negative_path, 
                 herb_name=None, config=None, use_parallel=True,
                 rt_tolerance=0.3, loss_tolerance=0.02,
                 external_diagnostic_file=None):
        """初始化鉴定程序（增强版）"""
        # 默认配置参数
        self.config = {
            'gradient_time': 30.0,
            'cid_min': 0.5,
            'cid_max': 1.0,
            'fragment_tolerance': 0.05,
            'tolerance_ppm': 50,
            'max_candidates': 3,          # 保留，用于调试
            'min_fragment_count': 1,
            'min_intensity': 100,
            'ppm_tier1': 10,
            'ppm_tier2': 20,
        }
        
        # 更新用户配置
        if config:
            self.config.update(config)
        
        # 新增参数
        self.rt_tolerance = rt_tolerance          # 保留时间容差（分钟）（虽已不使用，但保留参数兼容）
        self.loss_tolerance = loss_tolerance      # 中性丢失质量容差（Da）
        
        # 目标药材名称
        self.herb_name = herb_name
        
        # 多进程设置
        self.use_parallel = use_parallel and os.cpu_count() > 1
        self.num_workers = min(os.cpu_count(), 8)
        
        # 加载数据文件
        print("="*80)
        print("中药化合物鉴定程序 v5.5（评分优化版）")
        print("="*80)
        print("\n【1/7】正在加载数据库...")
        self.full_database = self._load_data(database_path)
        
        # 根据药材名称筛选数据库
        if herb_name:
            print(f"【2/7】正在筛选 {herb_name} 相关数据...")
            self.database = self._filter_by_herb(herb_name)
        else:
            self.database = self.full_database.copy()
            print("【2/7】使用全部数据库进行化合物鉴定")
        
        print(f"  数据库记录数: {len(self.database)}")
        
        # 加载质谱数据
        print("【3/7】正在加载质谱数据...")
        self.ms_positive = self._load_data(ms_positive_path)
        self.ms_negative = self._load_data(ms_negative_path)
        
        # 构建优化索引
        print("【4/7】正在构建索引...")
        self._build_optimized_index()
        
        # 构建诊断性离子库（支持外部文件）
        print("【5/7】正在加载诊断离子库...")
        # 如果未提供外部文件，尝试自动查找项目目录下的诊断离子.xlsx
        if external_diagnostic_file is None:
            default_diag_path = find_diagnostic_ion_path()
            if default_diag_path:
                print(f"  自动找到默认诊断离子文件: {default_diag_path}")
                external_diagnostic_file = default_diag_path
        self._build_diagnostic_ion_library(external_diagnostic_file)
        
        # 加载辅助数据（保留时间、裂解规则等）
        print("【6/7】正在加载辅助数据...")
        self._load_auxiliary_data()
        
        # 统计信息
        self.stats = {
            'total_precursors': 0,
            'identified_compounds': 0,
            'scored_candidates': 0
        }
        
        # 打印初始化信息
        self._print_initialization_info()
        print("【7/7】初始化完成，准备鉴定")
    
    def _load_data(self, filepath):
        """加载数据文件（支持Excel格式）"""
        if filepath and os.path.exists(filepath):
            try:
                if filepath.endswith('.xlsx'):
                    return pd.read_excel(filepath)
                elif filepath.endswith('.csv'):
                    return pd.read_csv(filepath)
            except Exception as e:
                print(f"警告: 无法加载文件 {filepath}: {e}")
                return pd.DataFrame()
        print(f"警告: 文件不存在或路径无效: {filepath}")
        return pd.DataFrame()
    
    def _filter_by_herb(self, herb_name):
        """按药材名称筛选数据库"""
        if herb_name is None:
            return self.full_database.copy()
        
        mask = self.full_database['药材名称'].str.contains(herb_name, na=False, case=False)
        filtered_db = self.full_database[mask].copy()
        
        if len(filtered_db) == 0:
            print(f"警告: 未找到 '{herb_name}' 相关数据，将使用全部数据")
            return self.full_database.copy()
        
        return filtered_db
    
    def _build_optimized_index(self):
        """构建优化索引（增强版：保留额外信息，并对分子式标准化）"""
        self.sorted_idx_pos = []
        self.sorted_idx_neg = []
        self.mz_values_pos = []
        self.mz_values_neg = []
        self.db_frag_pos = []
        self.db_frag_neg = []
        self.compound_info = {}
        
        # 额外存储用于评分的字段（此处保留中性丢失，但暂不使用）
        self.neutral_losses = {}      # 中性丢失列表，键为db_idx
        
        for idx, row in self.database.iterrows():
            mz_pos = row.get('准分子离子（正）', 0)
            mz_neg = row.get('准分子离子（负）', 0)
            
            cn_name = str(row.get('名称（中文）', ''))
            en_name = str(row.get('名称（英文）', ''))
            formula = str(row.get('分子式', ''))
            # 标准化分子式
            if formula and formula != 'nan':
                formula = normalize_formula(formula)
            herb = str(row.get('药材名称', ''))
            compound_type = str(row.get('化合物类型', ''))
            source = str(row.get('文献来源', ''))
            cas = str(row.get('CAS', ''))
            
            self.compound_info[idx] = {
                'name_cn': cn_name,
                'name_en': en_name,
                'formula': formula,
                'herb': herb,
                'compound_type': compound_type,
                'source': source,
                'cas': cas if cas and cas != 'nan' else '',
                'adduct_pos': str(row.get('加合物（正）', '')),
                'adduct_neg': str(row.get('加合物（负）', ''))
            }
            
            # 加载中性丢失（如果存在）
            if '中性丢失' in row and pd.notna(row['中性丢失']):
                losses = self._parse_losses(row['中性丢失'])
                if losses:
                    self.neutral_losses[idx] = losses
            
            # 正离子索引
            if pd.notna(mz_pos) and str(mz_pos).strip() != '':
                try:
                    mz_val = float(mz_pos)
                    if mz_val > 0:
                        fragments = self._parse_fragments_fast(row.get('碎片离子（正）', ''))
                        self.sorted_idx_pos.append((mz_val, idx, fragments))
                except (ValueError, TypeError):
                    pass
            
            # 负离子索引
            if pd.notna(mz_neg) and str(mz_neg).strip() != '':
                try:
                    mz_val = float(mz_neg)
                    if mz_val > 0:
                        fragments = self._parse_fragments_fast(row.get('碎片离子（负）', ''))
                        self.sorted_idx_neg.append((mz_val, idx, fragments))
                except (ValueError, TypeError):
                    pass
        
        # 排序并提取分子量数组和碎片数组
        self.sorted_idx_pos.sort(key=lambda x: x[0])
        self.sorted_idx_neg.sort(key=lambda x: x[0])
        
        self.mz_values_pos = np.array([x[0] for x in self.sorted_idx_pos])
        self.mz_values_neg = np.array([x[0] for x in self.sorted_idx_neg])
        
        self.db_frag_pos = [x[2] for x in self.sorted_idx_pos]
        self.db_frag_neg = [x[2] for x in self.sorted_idx_neg]
    
    def _parse_fragments_fast(self, fragment_string):
        """快速解析碎片离子"""
        if pd.isna(fragment_string) or not fragment_string or str(fragment_string).strip() == '':
            return np.array([])
        
        fragments = []
        try:
            for part in str(fragment_string).split(','):
                part = part.strip()
                if part:
                    try:
                        fragments.append(float(part))
                    except ValueError:
                        continue
        except Exception:
            pass
        return np.array(fragments)
    
    def _build_diagnostic_ion_library(self, external_file=None):
        """构建诊断性离子库（支持外部Excel文件）"""
        if external_file and os.path.exists(external_file):
            try:
                df = pd.read_excel(external_file)
                # 要求Excel文件至少包含 '化合物类型' 和 '诊断碎片离子m/z' 两列
                required_cols = ['化合物类型', '诊断碎片离子m/z']
                if not all(col in df.columns for col in required_cols):
                    raise ValueError(f"外部诊断离子文件缺少必要列：{required_cols}")
                
                self.diagnostic_ions = {}
                for category, group in df.groupby('化合物类型'):
                    ions = group['诊断碎片离子m/z'].dropna().tolist()
                    # 可选添加描述列，若无则使用默认描述
                    description = group['描述'].iloc[0] if '描述' in group.columns else f'来自外部文件，{len(ions)}个离子'
                    self.diagnostic_ions[category] = {
                        'ions': ions,
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
        """默认内置诊断离子库"""
        self.diagnostic_ions = {
            '环烯醚萜类': {
                'ions': [138.055, 124.039, 110.023, 96.008, 82.029, 67.029, 127.039],
                'description': '环烯醚萜类特征脱水碎片'
            },
            '有机酸类': {
                'ions': [191.056, 179.034, 173.045, 135.045, 93.034, 85.029],
                'description': '咖啡酰奎尼酸系列特征离子'
            },
            '黄酮类': {
                'ions': [151.003, 137.024, 121.029, 107.049, 81.034, 65.039],
                'description': '黄酮苷元特征碎片'
            },
            '萜类': {
                'ions': [127.076, 113.060, 99.044, 85.029, 71.013],
                'description': '萜类特征碎片'
            },
            '栀子特异': {
                'ions': [127.039, 113.024, 101.024, 69.034, 97.028],
                'description': '栀子特有成分特征离子'
            },
            '生物碱类': {
                'ions': [105.070, 91.054, 79.054, 65.039],
                'description': '生物碱特征碎片'
            },
            '酚酸类': {
                'ions': [137.024, 123.044, 109.028, 95.049],
                'description': '酚酸类特征碎片'
            }
        }
    
    def _load_auxiliary_data(self):
        """加载辅助数据（如外部保留时间表、裂解规则）"""
        # 预留扩展点
        pass
    
    def _parse_losses(self, loss_string):
        """解析中性丢失字符串（支持逗号分隔）"""
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
        """从碎片离子中找出可能的中性丢失质量数"""
        if len(fragments) == 0:
            return []
        observed_losses = []
        precursor_mz = float(precursor_mz)
        for frag in fragments:
            loss = precursor_mz - float(frag)
            if 10 < loss < 500:  # 合理的中性丢失范围
                observed_losses.append(round(loss, 4))
        return observed_losses
    
    def _classify_compound(self, name, compound_type):
        """分类化合物（用于匹配诊断离子库中的类别）"""
        name = name.lower()
        compound_type = str(compound_type).lower()
        
        # 优先根据compound_type判断
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
        
        # 根据名称判断栀子特异
        if any(keyword in name for keyword in ['京尼平', '栀子', 'genipos', 'gardenia']):
            return '栀子特异'
        
        # 如果都不匹配，返回'其他'，但诊断离子库中可能没有这个类别，匹配时会跳过
        return '其他'
    
    def _binary_search_range(self, mz_array, mz, tolerance_ppm):
        """二分查找快速定位候选化合物"""
        if len(mz_array) == 0:
            return range(0, 0)
        
        mz = float(mz)
        tolerance = mz * tolerance_ppm / 1e6
        mz_min = mz - tolerance
        mz_max = mz + tolerance
        
        left = bisect_left(mz_array, mz_min)
        right = bisect_right(mz_array, mz_max)
        
        return range(left, right)
    
    def _match_fragments_fast(self, observed, reference, tolerance=0.05):
        """快速匹配碎片离子"""
        matched = []
        if len(reference) == 0:
            return matched
        
        ref_arr = np.asarray(reference)
        obs_arr = np.asarray(observed)
        
        for ref_val in ref_arr:
            if pd.notna(ref_val) and float(ref_val) > 0:
                for obs_val in obs_arr:
                    if abs(float(obs_val) - float(ref_val)) <= tolerance:
                        matched.append(float(obs_val))
                        break
        
        return list(set(matched))
    
    def _find_diagnostic_ions_fast(self, matched_fragments, category):
        """快速查找诊断性离子（使用当前加载的诊断离子库）"""
        if len(matched_fragments) == 0:
            return []
        
        if category not in self.diagnostic_ions:
            return []
        
        diag_ions = np.array(self.diagnostic_ions[category]['ions'])
        matched_arr = np.asarray(matched_fragments)
        
        diagnostic = []
        for diag_val in diag_ions:
            for matched_val in matched_arr:
                if abs(float(matched_val) - float(diag_val)) <= self.config['fragment_tolerance']:
                    diagnostic.append(float(matched_val))
                    break
        
        return list(set(diagnostic))
    
    def _determine_confidence_level(self, ppm, matched_count, diagnostic_count, has_fragment_data):
        """确定置信度等级（保留原评级）"""
        ppm_tier1 = self.config.get('ppm_tier1', 10)
        ppm_tier2 = self.config.get('ppm_tier2', 20)
        
        if ppm > 50:
            return 6, '排除级', '未识别', 'ppm > 50，不符合报告要求'
        
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
        
        return 6, '排除级', '未识别', '不符合评级标准'
    
    def _calculate_score_by_rating(self, rating, ppm, matched_frag_count, diag_count):
        """
        根据评级和原始指标计算综合得分
        确保1级≥80，2级≥60，且级别越小得分越高
        """
        # 基础分
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
        
        # ppm加减
        if ppm <= 5:
            ppm_adj = 5
        elif ppm <= 10:
            ppm_adj = 0
        elif ppm <= 20:
            ppm_adj = -5
        elif ppm <= 30:
            ppm_adj = -10
        else:  # ppm <= 50
            ppm_adj = -15
        
        # 碎片加分（每个匹配碎片+2分，上限20）
        frag_adj = min(matched_frag_count * 2, 20)
        
        # 诊断离子加分（每个+5分，上限15）
        diag_adj = min(diag_count * 5, 15)
        
        total = base + ppm_adj + frag_adj + diag_adj
        
        # 保底机制
        if rating == 1 and total < 80:
            total = 80
        if rating == 2 and total < 60:
            total = 60
        
        return round(total, 2)
    
    def extract_precursor_ions(self, ms_data, ionization_mode):
        """从质谱数据中提取母离子信息"""
        precursors = []
        mz_columns = [col for col in ms_data.columns if 'Peak_' in col and '_m/z' in col]
        min_intensity = self.config['min_intensity']
        
        for idx, row in ms_data.iterrows():
            precursor_mz = row.get('Precursor M/z')
            if pd.notna(precursor_mz) and float(precursor_mz) > 0:
                fragments = []
                fragments_dict = {}
                
                for col in mz_columns:
                    fragment_mz = row[col]
                    if pd.notna(fragment_mz) and float(fragment_mz) > 0:
                        intensity_col = col.replace('_m/z', '_Intensity')
                        intensity = 0
                        if intensity_col in row.index and pd.notna(row[intensity_col]):
                            intensity = float(row[intensity_col])
                        
                        if intensity >= min_intensity:
                            mz_value = float(fragment_mz)
                            fragments.append(mz_value)
                            fragments_dict[mz_value] = intensity
                
                # 尝试获取保留时间：先找“出峰时间t/min”列，再找“出峰时间”
                rt = row.get('出峰时间t/min', np.nan)
                if pd.isna(rt):
                    rt = row.get('出峰时间', np.nan)
                
                # 如果仍为空且存在CID列，则模拟计算
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
                    'base_peak': max(fragments_dict.values()) if fragments_dict else 0
                })
        
        return precursors
    
    def identify_compound(self, precursor, return_top=1, score_threshold=0):
        """
        鉴定单个化合物，返回得分最高的候选（增强版）
        
        参数:
            precursor: 母离子信息
            return_top: 返回前几个候选（默认1，设为>1返回多个，用于调试）
            score_threshold: 得分阈值，低于此值不返回（原得分已弃用，设默认0）
        
        返回:
            若return_top==1: 返回最佳候选字典（或None）
            若return_top>1: 返回候选列表，按得分排序
        """
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
            
            matched_fragments = self._match_fragments_fast(precursor['fragments'], ref_frags, self.config['fragment_tolerance'])
            diagnostic_ions = self._find_diagnostic_ions_fast(matched_fragments, category)
            
            theoretical_mz = mz_val
            ppm = abs(float(precursor['precursor_mz']) - theoretical_mz) / theoretical_mz * 1e6
            
            candidate = {
                'name_cn': info['name_cn'],
                'name_en': info['name_en'],
                'formula': info['formula'],  # 已标准化
                'cas': info['cas'],
                'herb': info['herb'],
                'compound_type': info['compound_type'],
                'observed_mz': precursor['precursor_mz'],
                'theoretical_mz': theoretical_mz,
                'ppm': ppm,
                'adduct': info['adduct_pos'] if mode == '正离子' else info['adduct_neg'],
                'matched_fragments': matched_fragments,
                'diagnostic_ions': diagnostic_ions,
                'mode': mode,
                'source': info['source'],
                'category': category,
                'db_index': db_idx  # 保留索引以便后续获取额外信息
            }
            
            # 旧版评分已弃用，这里只保留候选，得分在生成报告时重新计算
            # 为保持排序，暂时用ppm和匹配碎片数作为临时排序依据
            candidate['temp_score'] = -ppm  # 临时排序用
            
            if True:  # 移除得分阈值过滤
                scored_candidates.append(candidate)
        
        # 按临时分数排序（仅用于调试时的返回多个）
        scored_candidates.sort(key=lambda x: (-len(x['matched_fragments']), -x['temp_score']))
        
        if return_top == 1:
            return scored_candidates[0] if scored_candidates else None
        else:
            return scored_candidates[:return_top]
    
    def generate_report(self, herb_name=None):
        """生成化合物鉴定报告（增强版：择优输出，合并所有药材来源，使用新评分）"""
        if herb_name is None:
            herb_name = self.herb_name if self.herb_name else '中药'
        
        # 用于存储最佳记录和所有药材集合
        best_records = {}      # 键: (name_cn, formula) -> (record, temp_score)
        herbs_collection = {}  # 键: (name_cn, formula) -> set of herbs
        
        # 处理正离子模式
        print(f"\n【6/7】正在处理 {herb_name} 正离子模式质谱数据...")
        pos_precursors = self.extract_precursor_ions(self.ms_positive, '正离子')
        print(f"  - 共提取 {len(pos_precursors)} 个母离子")
        
        # 处理负离子模式
        print(f"\n【7/7】正在处理 {herb_name} 负离子模式质谱数据...")
        neg_precursors = self.extract_precursor_ions(self.ms_negative, '负离子')
        print(f"  - 共提取 {len(neg_precursors)} 个母离子")
        
        all_precursors = pos_precursors + neg_precursors
        total = len(all_precursors)
        
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
            
            # 如果评级为排除级，跳过
            if lvl_id == 6:
                continue
            
            formula = best_candidate['formula'] if best_candidate['formula'] and best_candidate['formula'] != 'nan' else '待确定'
            en_name = best_candidate['name_en'] if best_candidate['name_en'] and best_candidate['name_en'] != 'nan' else ''
            
            # 处理文献来源：拆分并计数
            source_str = best_candidate.get('source', '')
            if pd.notna(source_str) and source_str:
                source_list = [s.strip() for s in str(source_str).split(';') if s.strip()]
                source_count = len(source_list)
                source_display = '; '.join(source_list)
            else:
                source_count = 0
                source_display = ''
            
            # 匹配上的碎片和诊断离子
            matched_frags = best_candidate['matched_fragments']
            diag_ions = best_candidate['diagnostic_ions']
            
            # 清洗加和离子：如果包含逗号或分号，只取第一个
            adduct = best_candidate['adduct']
            if adduct:
                for sep in [',', ';']:
                    if sep in adduct:
                        adduct = adduct.split(sep)[0].strip()
                        break
            
            # 计算新综合得分
            new_score = self._calculate_score_by_rating(
                lvl_id, 
                best_candidate['ppm'],
                len(matched_frags),
                len(diag_ions)
            )
            
            record = {
                '序号': 0,
                '出峰时间t/min': precursor['retention_time'],
                '化合物中文名': best_candidate['name_cn'],
                '化合物英文名': en_name,
                '分子式': formula,
                'CAS号': best_candidate['cas'],
                '药材名称': best_candidate['herb'],  # 临时值，后面会合并
                '化合物类型': best_candidate['compound_type'],
                '离子化方式': best_candidate['mode'],
                '加和离子': adduct,
                'm/z实际值': round(best_candidate['observed_mz'], 4),
                'm/z理论值': round(best_candidate['theoretical_mz'], 4),
                'ppm': round(best_candidate['ppm'], 4),
                '是否有碎片数据': '是' if precursor['fragments'].size > 0 else '否',
                '主要碎片离子': '; '.join([f'{f:.4f}' for f in matched_frags]) if matched_frags else '',
                '匹配碎片数': len(matched_frags),
                '诊断性离子个数': len(diag_ions),
                '诊断性离子': '; '.join([f'{f:.4f}' for f in diag_ions]) if diag_ions else '',
                '文献来源数': source_count,
                '文献来源': source_display,
                '评级': lvl_id,
                '评级名称': lvl_name,
                '置信度': confidence,
                '报告建议': suggestion,
                '综合得分': new_score
            }
            
            compound_key = (best_candidate['name_cn'], formula)
            
            # 收集药材名称
            if compound_key not in herbs_collection:
                herbs_collection[compound_key] = set()
            herbs_collection[compound_key].add(best_candidate['herb'])
            
            # 更新最佳记录（使用新得分作为择优依据）
            current_score = new_score
            if compound_key not in best_records or current_score > best_records[compound_key][1]:
                best_records[compound_key] = (record, current_score)
        
        # 合并药材名称并生成最终报告
        results = []
        for compound_key, (record, _) in best_records.items():
            herbs = herbs_collection[compound_key]
            # 合并药材名称（去重后按字母排序）
            record['药材名称'] = '; '.join(sorted(herbs))
            results.append(record)
        
        # 生成DataFrame并排序
        report_df = pd.DataFrame(results)
        
        if not report_df.empty:
            # 按评级和得分排序
            report_df = report_df.sort_values(by=['评级', '综合得分', 'ppm'], ascending=[True, False, True])
            report_df = report_df.reset_index(drop=True)
            report_df['序号'] = range(1, len(report_df) + 1)
        
        return report_df
    
    def save_report(self, report_df, output_path):
        """保存报告"""
        report_df.to_excel(output_path, index=False)
        print(f"\n报告已保存至: {output_path}")
    
    def print_summary(self, report_df, herb_name):
        """打印报告摘要（增强版：增加得分统计）"""
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
            score_high = (report_df['综合得分'] >= 80).sum()
            score_medium = ((report_df['综合得分'] >= 60) & (report_df['综合得分'] < 80)).sum()
            score_low = (report_df['综合得分'] < 60).sum()
            print(f"  - 高分 (≥80): {score_high} 个")
            print(f"  - 中分 (60-80): {score_medium} 个")
            print(f"  - 低分 (<60): {score_low} 个")
            
            print(f"\n【药材来源分布（前10）】")
            for herb, count in report_df['药材名称'].value_counts().head(10).items():
                print(f"  - {herb}: {count} 个")
        
        print("\n" + "="*100)
    
    def _print_initialization_info(self):
        """打印初始化信息"""
        print("\n" + "="*80)
        print("程序初始化完成（评分优化版）")
        print("="*80)
        print(f"  - 数据库记录数: {len(self.database)} 条")
        print(f"  - 正离子索引: {len(self.mz_values_pos)} 条")
        print(f"  - 负离子索引: {len(self.mz_values_neg)} 条")
        print(f"  - 正离子质谱: {len(self.ms_positive)} 条记录")
        print(f"  - 负离子质谱: {len(self.ms_negative)} 条记录")
        print(f"  - 诊断性离子库: {len(self.diagnostic_ions)} 类")
        print(f"  - 中性丢失数据: {len(self.neutral_losses)} 条")
        print(f"  - 并行处理: {'启用' if self.use_parallel else '禁用'}")
        print("="*80)


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


def deduplicate_report(input_file, output_file):
    """集成去重功能（保留，但增强版报告已内嵌去重）"""
    print('\n' + '='*80)
    print('【集成去重处理】')
    print('='*80)
    
    print(f'\n正在读取报告: {input_file}')
    df = pd.read_excel(input_file)
    print(f'  原始记录数: {len(df)}')
    
    # 标准化分子式
    print('\n步骤1: 标准化分子式（去除Unicode下标）')
    df['分子式_标准化'] = df['分子式'].apply(normalize_formula)
    changes = (df['分子式'] != df['分子式_标准化']).sum()
    print(f'  分子式标准化变化数: {changes}')
    
    # 合并同一化合物
    print('\n步骤2: 合并同一化合物（合并正负离子模式结果）')
    df['化合物ID'] = df['化合物中文名'].str.strip() + '_' + df['分子式_标准化'].str.strip()
    
    def select_best(group):
        sorted_group = group.sort_values(by=['评级', 'ppm'], ascending=[True, True])
        return sorted_group.iloc[0]
    
    def aggregate_modes(group):
        best = select_best(group)
        modes = group['离子化方式'].unique()
        best['离子化方式'] = '/'.join(modes)
        return best
    
    df_dedup = df.groupby('化合物ID', as_index=False).apply(aggregate_modes, include_groups=False)
    print(f'  合并后记录数: {len(df_dedup)}')
    
    # 最终去重
    print('\n步骤3: 最终去重（保留最佳评级和最小ppm）')
    df_dedup = df_dedup.sort_values(by=['评级', 'ppm'], ascending=[True, True])
    df_dedup = df_dedup.drop_duplicates(subset=['化合物ID'], keep='first')
    print(f'  最终记录数: {len(df_dedup)}')
    
    # 清理并保存
    print('\n步骤4: 生成最终报告')
    for col in ['化合物ID', '分子式_标准化']:
        if col in df_dedup.columns:
            df_dedup = df_dedup.drop(columns=[col])
    
    df_dedup = df_dedup.sort_values(by=['评级', 'ppm'], ascending=[True, True])
    df_dedup = df_dedup.reset_index(drop=True)
    df_dedup['序号'] = range(1, len(df_dedup) + 1)
    
    df_dedup.to_excel(output_file, index=False)
    print(f'  报告已保存至: {output_file}')
    
    # 统计报告
    print('\n' + '='*80)
    print('去重统计报告')
    print('='*80)
    print(f'原始记录数: {len(df)}')
    print(f'最终记录数: {len(df_dedup)}')
    print(f'去除重复: {len(df) - len(df_dedup)} 条 ({100*(len(df) - len(df_dedup))/len(df):.1f}%)')
    
    print('\n【评级分布】')
    for level in ['确证级', '高置信级', '推定级', '提示级', '参考级']:
        count = (df_dedup['评级名称'] == level).sum()
        print(f'  - {level}: {count} 个')
    
    print('\n【ppm误差分布】')
    ppm_10 = (df_dedup['ppm'] <= 10).sum()
    ppm_20 = ((df_dedup['ppm'] > 10) & (df_dedup['ppm'] <= 20)).sum()
    ppm_50 = ((df_dedup['ppm'] > 20) & (df_dedup['ppm'] <= 50)).sum()
    print(f'  - ≤10ppm: {ppm_10} 个 ({100*ppm_10/len(df_dedup):.1f}%)')
    print(f'  - 10-20ppm: {ppm_20} 个 ({100*ppm_20/len(df_dedup):.1f}%)')
    print(f'  - 20-50ppm: {ppm_50} 个 ({100*ppm_50/len(df_dedup):.1f}%)')
    
    print('\n【药材来源分布（前10）】')
    for herb, count in df_dedup['药材名称'].value_counts().head(10).items():
        print(f'  - {herb}: {count} 个')
    
    print('\n【确证级化合物列表（前15）】')
    confirmed = df_dedup[df_dedup['评级'] == 1].head(15)
    for _, row in confirmed.iterrows():
        print(f'  - {row["化合物中文名"]} ({row["分子式"]}) - ppm={row["ppm"]:.4f}')
    
    return df_dedup


# ============================================================================
# 云端优化：数据库加载函数（带缓存）
# ============================================================================

@st.cache_data
def load_database_cached(db_filename="TCM-SM-MS DB.xlsx"):
    """加载数据库（云端优化版本，使用缓存）"""
    # 可能的数据库文件名列表
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
                return df
            except Exception as e:
                continue
    
    return pd.DataFrame()


def find_database_path():
    """查找数据库文件路径"""
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


@st.cache_data
def load_diagnostic_ions_cached():
    """加载诊断离子数据库（带缓存）"""
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


def match_diagnostic_ions(user_mz_values, diagnostic_df, tolerance_ppm=10, ion_mode=None):
    """匹配诊断离子（用于单独筛查功能）"""
    if diagnostic_df.empty or not user_mz_values:
        return pd.DataFrame()
    
    # 过滤离子模式
    if ion_mode and ion_mode != "全部":
        filtered_df = diagnostic_df[diagnostic_df['离子模式'] == ion_mode].copy()
    else:
        filtered_df = diagnostic_df.copy()
    
    # 确保m/z列为数值类型
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

# 设置页面配置
st.set_page_config(
    page_title="中药化合物智能鉴定平台 v5.5（评分优化版）",
    page_icon="🌿",
    layout="wide",
    initial_sidebar_state="expanded"
)

# 设置中文字体支持
st.markdown("""
<style>
    @import url('https://fonts.googleapis.com/css2?family=Noto+Sans+SC:wght@300;400;500;700&display=swap');
    
    * {
        font-family: 'Noto Sans SC', 'Microsoft YaHei', sans-serif !important;
    }
    
    .main-header {
        background: linear-gradient(135deg, #2E7D32 0%, #1976D2 100%);
        padding: 2rem;
        border-radius: 10px;
        margin-bottom: 2rem;
        color: white;
    }
    
    .feature-card {
        background: white;
        border-radius: 10px;
        padding: 1.5rem;
        box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        margin-bottom: 1rem;
        border-left: 4px solid #2E7D32;
    }
    
    .stat-box {
        background: linear-gradient(135deg, #f5f7fa 0%, #e4e8eb 100%);
        border-radius: 8px;
        padding: 1rem;
        text-align: center;
        margin: 0.5rem 0;
    }
    
    .stat-number {
        font-size: 2rem;
        font-weight: 700;
        color: #2E7D32;
    }
    
    .stat-label {
        font-size: 0.9rem;
        color: #666;
    }
    
    .level-1 { background-color: #28a745; color: white; }
    .level-2 { background-color: #17a2b8; color: white; }
    .level-3 { background-color: #ffc107; color: black; }
    .level-4 { background-color: #fd7e14; color: white; }
    .level-5 { background-color: #6c757d; color: white; }
    
    .stDataFrame {
        border-radius: 8px;
    }
    
    div[data-testid="stExpander"] {
        border: 1px solid #e0e0e0;
        border-radius: 8px;
    }
</style>
""", unsafe_allow_html=True)


def load_css():
    """加载自定义CSS样式"""
    st.markdown("""
    <style>
        /* 整体应用样式 */
        .stApp {
            background: linear-gradient(180deg, #f8f9fa 0%, #ffffff 100%);
        }
        
        /* 标题样式 */
        h1, h2, h3 {
            color: #2E7D32 !important;
            font-weight: 600 !important;
        }
        
        /* 侧边栏样式 */
        [data-testid="stSidebar"] {
            background: linear-gradient(180deg, #ffffff 0%, #f5f7fa 100%);
            border-right: 1px solid #e0e0e0;
        }
        
        /* 按钮样式 */
        .stButton > button {
            background: linear-gradient(135deg, #2E7D32 0%, #1976D2 100%);
            color: white;
            border: none;
            border-radius: 8px;
            padding: 0.5rem 1rem;
            font-weight: 500;
            transition: all 0.3s ease;
        }
        
        .stButton > button:hover {
            transform: translateY(-2px);
            box-shadow: 0 4px 12px rgba(46, 125, 50, 0.4);
        }
        
        /* 进度条样式 */
        .stProgress > div > div {
            background: linear-gradient(90deg, #2E7D32 0%, #1976D2 100%);
        }
        
        /* 文件上传样式 */
        [data-testid="stFileUploader"] {
            border: 2px dashed #2E7D32;
            border-radius: 10px;
            padding: 2rem;
            background: rgba(46, 125, 50, 0.05);
        }
        
        /* 表格样式 */
        .dataframe {
            font-size: 0.9rem;
        }
        
        /* 指标卡片 */
        [data-testid="stMetricValue"] {
            font-size: 1.5rem !important;
            color: #2E7D32;
        }
        
        /* 警告和信息框 */
        .stAlert {
            border-radius: 8px;
        }
    </style>
    """, unsafe_allow_html=True)


def create_header():
    """创建应用头部"""
    st.markdown("""
    <div class="main-header">
        <h1 style="color: white !important; margin: 0;">🌿 中药化合物智能鉴定平台（评分优化版）</h1>
        <p style="margin: 0.5rem 0 0 0; opacity: 0.9;">基于评级的新评分系统 v5.5（云端版）</p>
    </div>
    """, unsafe_allow_html=True)


def create_sidebar():
    """创建侧边栏导航"""
    st.sidebar.markdown("""
    <div style="text-align: center; padding: 1rem 0;">
        <h2 style="color: #2E7D32; margin-bottom: 0.5rem;">🔬 TCM Identifier</h2>
        <p style="color: #666; font-size: 0.8rem;">中药化合物鉴定系统（评分优化版）</p>
    </div>
    """, unsafe_allow_html=True)
    
    page = st.sidebar.radio(
        "导航菜单",
        ["首页", "开始鉴定", "诊断离子筛查", "使用指南", "数据库预览", "结果分析"]
    )
    
    st.sidebar.markdown("---")
    
    st.sidebar.info("""
    **版本信息**
    - 程序版本：v5.5（评分优化版）
    - 数据库规模：35,828条化合物记录
    - 诊断离子：自动加载项目目录下的诊断离子.xlsx
    - 支持药材：291种
    - 核心特点：新评分系统（1级≥80，2级≥60）、药材名称合并、外部诊断离子
    """)
    
    st.sidebar.markdown("""
    <div style="text-align: center; color: #999; font-size: 0.7rem; padding: 1rem 0;">
        <p>© 2026 张永</p>
        <p>中药化合物智能鉴定平台</p>
    </div>
    """, unsafe_allow_html=True)
    
    return page


def show_home_page():
    """首页"""
    create_header()
    
    st.markdown("## 📊 系统概览")
    
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.markdown("""
        <div class="stat-box">
            <div class="stat-number">35,828</div>
            <div class="stat-label">数据库化合物数</div>
        </div>
        """, unsafe_allow_html=True)
    with col2:
        st.markdown("""
        <div class="stat-box">
            <div class="stat-number">291</div>
            <div class="stat-label">支持药材种类</div>
        </div>
        """, unsafe_allow_html=True)
    with col3:
        st.markdown("""
        <div class="stat-box">
            <div class="stat-number">6</div>
            <div class="stat-label">鉴定评级级别</div>
        </div>
        """, unsafe_allow_html=True)
    with col4:
        st.markdown("""
        <div class="stat-box">
            <div class="stat-number">400+</div>
            <div class="stat-label">化合物类型</div>
        </div>
        """, unsafe_allow_html=True)
    
    st.markdown("---")
    st.markdown("## ✨ 核心功能")
    
    col1, col2 = st.columns(2)
    with col1:
        st.markdown("""
        <div class="feature-card">
            <h4>🔍 智能化合物鉴定</h4>
            <p>基于高分辨质谱数据，在35,828条数据库记录中进行精准匹配，支持正负离子模式。</p>
        </div>
        <div class="feature-card">
            <h4>📈 六级评级标准</h4>
            <p>确证级、高置信级、推定级、提示级、参考级、排除级，科学评估鉴定结果可靠性。</p>
        </div>
        <div class="feature-card">
            <h4>🧪 外部诊断离子</h4>
            <p>支持上传自定义诊断离子文件，也可自动加载项目目录下的诊断离子.xlsx。</p>
        </div>
        """, unsafe_allow_html=True)
    with col2:
        st.markdown("""
        <div class="feature-card">
            <h4>⚡ 全新评分系统</h4>
            <p>基于评级计算综合得分，确保1级≥80分，2级≥60分，级别越小得分越高。</p>
        </div>
        <div class="feature-card">
            <h4>🌱 药材来源合并</h4>
            <p>同一化合物的所有药材来源自动合并显示，实现“所有的药材”。</p>
        </div>
        <div class="feature-card">
            <h4>🚀 云端优化</h4>
            <p>采用缓存和并行处理，支持Streamlit Community Cloud一键部署。</p>
        </div>
        """, unsafe_allow_html=True)
    
    st.markdown("---")
    st.markdown("## 🚀 快速开始")
    st.info("""
    ### 使用步骤
    1. **准备数据**：准备好质谱数据文件（Excel格式）和数据库文件（TCM-SM-MS DB.xlsx）
    2. **上传文件**：在"开始鉴定"页面上传正负离子质谱数据
    3. **可选**：上传自定义诊断离子文件（若不传，将自动加载项目目录下的诊断离子.xlsx）
    4. **配置参数**：设置ppm容差、保留时间容差等参数
    5. **运行鉴定**：点击"开始鉴定"按钮进行分析
    6. **查看结果**：在"结果分析"页面查看和导出鉴定报告
    """)
    
    if st.button("立即开始鉴定 →", type="primary"):
        st.session_state['page'] = '开始鉴定'
        st.rerun()


def show_analysis_page():
    """鉴定分析页面（增强版：增加外部诊断离子上传）"""
    create_header()
    
    st.markdown("## 📁 上传质谱数据")
    
    col1, col2 = st.columns(2)
    with col1:
        ms_positive_file = st.file_uploader("上传正离子模式质谱数据 (.xlsx)", type=['xlsx'], key='ms_positive')
        if ms_positive_file:
            st.success(f"✅ 已上传: {ms_positive_file.name}")
    with col2:
        ms_negative_file = st.file_uploader("上传负离子模式质谱数据 (.xlsx)", type=['xlsx'], key='ms_negative')
        if ms_negative_file:
            st.success(f"✅ 已上传: {ms_negative_file.name}")
    
    st.markdown("---")
    st.markdown("## 🧪 外部诊断离子（可选）")
    st.info("若不上传，将自动加载项目目录下的诊断离子.xlsx（如果存在）")
    diagnostic_file = st.file_uploader("上传自定义诊断离子文件 (.xlsx，格式：化合物类型、诊断碎片离子m/z)", type=['xlsx'], key='diagnostic')
    if diagnostic_file:
        st.success(f"✅ 已上传诊断离子文件: {diagnostic_file.name}")
    
    st.markdown("---")
    st.markdown("## ⚙️ 鉴定参数配置")
    
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        tolerance_ppm = st.number_input("ppm误差容限", min_value=10, max_value=100, value=50)
    with col2:
        max_candidates = st.number_input("最大候选数（调试用）", min_value=1, max_value=10, value=3)
    with col3:
        min_intensity = st.number_input("最小峰强度阈值", min_value=0, value=100)
    with col4:
        use_parallel = st.checkbox("启用多进程并行", value=True)
    
    col1, col2, col3 = st.columns(3)
    with col1:
        fragment_tolerance = st.number_input("碎片匹配容差 (Da)", min_value=0.01, max_value=1.0, value=0.05, step=0.01)
    with col2:
        rt_tolerance = st.number_input("保留时间容差 (min)", min_value=0.1, max_value=2.0, value=0.3, step=0.1)
    with col3:
        loss_tolerance = st.number_input("中性丢失容差 (Da)", min_value=0.01, max_value=0.5, value=0.02, step=0.01)
    
    herb_filter = st.selectbox("药材筛选模式", options=["使用全部数据库", "筛选特定药材"], index=0)
    if herb_filter == "筛选特定药材":
        herb_name = st.text_input("输入药材名称", placeholder="如：栀子、黄芩")
    else:
        herb_name = None
    
    st.markdown("---")
    
    if ms_positive_file and ms_negative_file:
        if st.button("🚀 开始化合物鉴定", type="primary", use_container_width=True):
            with st.spinner("正在初始化鉴定程序..."):
                try:
                    temp_dir = tempfile.gettempdir()
                    pos_path = os.path.join(temp_dir, ms_positive_file.name)
                    neg_path = os.path.join(temp_dir, ms_negative_file.name)
                    
                    with open(pos_path, 'wb') as f:
                        f.write(ms_positive_file.getbuffer())
                    with open(neg_path, 'wb') as f:
                        f.write(ms_negative_file.getbuffer())
                    
                    # 保存诊断离子文件（如果有）
                    diag_path = None
                    if diagnostic_file:
                        diag_path = os.path.join(temp_dir, diagnostic_file.name)
                        with open(diag_path, 'wb') as f:
                            f.write(diagnostic_file.getbuffer())
                    
                    db_path = find_database_path()
                    if not db_path:
                        st.error("未找到数据库文件！请将 TCM-SM-MS DB.xlsx 放在项目目录下。")
                        return
                    
                    st.info(f"✅ 已找到数据库文件: {db_path}")
                    
                    progress_bar = st.progress(0)
                    status_text = st.empty()
                    
                    status_text.text("【1/7】正在加载数据库...")
                    progress_bar.progress(10)
                    
                    identifier = UltimateGardeniaIdentifier(
                        database_path=db_path,
                        ms_positive_path=pos_path,
                        ms_negative_path=neg_path,
                        herb_name=herb_name,
                        use_parallel=use_parallel,
                        rt_tolerance=rt_tolerance,
                        loss_tolerance=loss_tolerance,
                        external_diagnostic_file=diag_path  # 若为None，将在内部自动查找默认文件
                    )
                    
                    status_text.text("【6/7】正在处理质谱数据...")
                    progress_bar.progress(40)
                    
                    report = identifier.generate_report('样品')
                    
                    progress_bar.progress(80)
                    status_text.text("【7/7】生成报告...")
                    
                    st.session_state['analysis_results'] = report
                    st.session_state['identifier'] = identifier
                    
                    progress_bar.progress(100)
                    status_text.text("鉴定完成！")
                    
                    st.success("✅ 化合物鉴定完成！")
                    
                    if st.button("查看鉴定结果 →"):
                        st.session_state['page'] = '结果分析'
                        st.rerun()
                    
                except Exception as e:
                    st.error(f"鉴定过程中出错：{str(e)}")
                    st.exception(e)
    else:
        st.warning("⚠️ 请先上传正负离子模式的质谱数据文件")


def show_diagnostic_ion_page():
    """诊断离子筛查页面（独立功能）"""
    create_header()
    
    st.markdown("## 🔬 诊断离子筛查")
    st.markdown("根据输入的m/z值，在诊断离子数据库中查找匹配的化合物特征离子，帮助快速识别化合物类别。")
    
    diagnostic_df = load_diagnostic_ions_cached()
    
    if diagnostic_df.empty:
        st.error("未找到诊断离子数据库文件（诊断离子.xlsx），请将文件放在项目目录下。")
        return
    
    st.success(f"✅ 已加载诊断离子数据库，包含 {len(diagnostic_df)} 条记录")
    
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        if '化合物类型' in diagnostic_df.columns:
            st.metric("化合物类型", diagnostic_df['化合物类型'].nunique())
    with col2:
        if '药材名' in diagnostic_df.columns:
            st.metric("药材种类", diagnostic_df['药材名'].nunique())
    with col3:
        if '诊断碎片离子m/z' in diagnostic_df.columns:
            st.metric("诊断离子数", diagnostic_df['诊断碎片离子m/z'].nunique())
    with col4:
        if '类特征性离子' in diagnostic_df.columns:
            char_ions = diagnostic_df[diagnostic_df['类特征性离子'] == True].shape[0]
            st.metric("类特征性离子", char_ions)
    
    st.markdown("---")
    st.markdown("### 📥 输入m/z值")
    
    col1, col2 = st.columns([3, 1])
    with col1:
        mz_input = st.text_area(
            "输入m/z值（每行一个值，或用逗号分隔）",
            placeholder="例如：\n151.003\n137.024, 121.029",
            height=150
        )
    with col2:
        st.markdown("#### 参数设置")
        tolerance_ppm = st.number_input("ppm容差", min_value=1, max_value=100, value=10)
        ion_mode = st.selectbox("离子模式", options=["全部", "正离子", "负离子"], index=0)
        show_only_class_specific = st.checkbox("仅显示类特征性离子", value=False)
    
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
            st.markdown("---")
            with st.spinner("正在匹配诊断离子..."):
                results_df = match_diagnostic_ions(user_mz_values, diagnostic_df, tolerance_ppm=tolerance_ppm, ion_mode=ion_mode)
            
            if show_only_class_specific and not results_df.empty:
                results_df = results_df[results_df['类特征性离子'] == True]
            
            st.markdown("### 📊 匹配结果统计")
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("输入离子数", len(user_mz_values))
            with col2:
                matched_ions = results_df['输入m/z'].nunique() if not results_df.empty else 0
                st.metric("匹配离子数", matched_ions)
            with col3:
                matched_compounds = results_df['化合物类型'].nunique() if not results_df.empty else 0
                st.metric("化合物类型数", matched_compounds)
            
            if not results_df.empty:
                st.markdown("#### 化合物类型分布")
                type_dist = results_df['化合物类型'].value_counts()
                st.bar_chart(type_dist)
                
                st.markdown("### 📋 匹配结果详情")
                display_cols = ['输入m/z', '匹配诊断离子m/z', '误差(ppm)', '化合物类型', '离子模式', '中文名称', '英文名称', '药材名', '类特征性离子']
                available_cols = [c for c in display_cols if c in results_df.columns]
                st.dataframe(results_df[available_cols], use_container_width=True, hide_index=True)
                
                col1, col2 = st.columns(2)
                with col1:
                    csv = results_df.to_csv(index=False, encoding='utf-8-sig')
                    st.download_button("📥 导出CSV", data=csv, file_name=f"诊断离子筛查_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv", mime="text/csv")
                with col2:
                    buffer = BytesIO()
                    with pd.ExcelWriter(buffer, engine='openpyxl') as writer:
                        results_df.to_excel(writer, index=False, sheet_name='诊断离子匹配结果')
                    st.download_button("📥 导出Excel", data=buffer.getvalue(), file_name=f"诊断离子筛查_{datetime.now().strftime('%Y%m%d_%H%M%S')}.xlsx", mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
            else:
                st.info("未找到匹配的诊断离子，请尝试增大ppm容差或更换离子模式。")


def show_guide_page():
    """使用指南页面"""
    create_header()
    st.markdown("## 📖 使用指南（评分优化版）")
    st.info("""
    ### 新增功能说明
    - **全新评分系统**：综合得分基于评级、ppm误差、匹配碎片数和诊断离子个数计算，确保1级化合物得分≥80，2级化合物得分≥60，级别越小得分越高。
    - **外部诊断离子**：在“开始鉴定”页面上传自定义诊断离子文件（Excel格式，需包含“化合物类型”和“诊断碎片离子m/z”列），程序将使用该库替代内置库进行诊断离子匹配。若不上传，程序会自动加载项目目录下的诊断离子.xlsx（如果存在）。
    - **择优输出**：每个母离子仅返回得分最高的候选，减少冗余。
    - **药材名称合并**：同一化合物的所有药材来源自动合并显示，实现“所有的药材”。
    
    ### 评级标准（同原六级标准）
    - 确证级、高置信级、推定级、提示级、参考级、排除级
    
    ### 使用步骤
    1. 上传正负离子质谱数据（Excel）
    2. （可选）上传自定义诊断离子文件
    3. 设置参数（ppm容差、保留时间容差等）
    4. 点击开始鉴定
    5. 在结果页面查看综合得分和评级
    """)


def show_database_page():
    """数据库预览页面"""
    create_header()
    st.markdown("## 🗃️ 数据库预览")
    
    db_path = find_database_path()
    if not db_path:
        st.error("未找到数据库文件！请将 TCM-SM-MS DB.xlsx 放在项目目录下。")
        return
    
    try:
        with st.spinner("正在加载数据库..."):
            df = load_database_cached()
        if df.empty:
            st.error("数据库文件为空或无法读取！")
            return
        
        st.success(f"✅ 成功加载数据库，包含 {len(df)} 条化合物记录")
        
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.metric("总记录数", len(df))
        with col2:
            herb_count = df['药材名称'].nunique() if '药材名称' in df.columns else 0
            st.metric("药材种类", herb_count)
        with col3:
            type_count = df['化合物类型'].nunique() if '化合物类型' in df.columns else 0
            st.metric("化合物类型", type_count)
        with col4:
            if '准分子离子（正）' in df.columns:
                valid_mz = df['准分子离子（正）'].notna().sum()
                st.metric("有效正离子记录", valid_mz)
        
        st.info(f"数据库文件路径: {db_path}")
        st.markdown("### 数据库预览（前10行）")
        st.dataframe(df.head(10), use_container_width=True)
        
    except Exception as e:
        st.error(f"加载数据库时出错：{str(e)}")
        st.exception(e)


def show_results_page():
    """结果分析页面（增强版）"""
    create_header()
    
    if 'analysis_results' not in st.session_state:
        st.warning("⚠️ 暂无鉴定结果，请先进行化合物鉴定。")
        if st.button("前往鉴定页面"):
            st.session_state['page'] = '开始鉴定'
            st.rerun()
        return
    
    report = st.session_state['analysis_results']
    
    st.markdown("## 📊 鉴定结果分析（评分优化版）")
    
    if report.empty:
        st.warning("鉴定结果为空，可能是因为没有匹配的化合物。")
        return
    
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.metric("鉴定化合物总数", len(report))
    with col2:
        confirmed = (report['评级名称'] == '确证级').sum()
        st.metric("确证级化合物", confirmed)
    with col3:
        high_score = (report['综合得分'] >= 80).sum()
        st.metric("高分化合物 (≥80)", high_score)
    with col4:
        avg_score = report['综合得分'].mean()
        st.metric("平均综合得分", f"{avg_score:.1f}")
    
    st.markdown("### 📊 评级分布")
    level_counts = report['评级名称'].value_counts()
    st.bar_chart(level_counts.reindex(['确证级', '高置信级', '推定级', '提示级', '参考级']).fillna(0))
    
    st.markdown("### 📈 综合得分分布")
    score_bins = pd.cut(report['综合得分'], bins=[0, 60, 70, 80, 90, 100], labels=['<60', '60-70', '70-80', '80-90', '90-100'])
    score_dist = score_bins.value_counts().sort_index()
    st.bar_chart(score_dist)
    
    st.markdown("### 📋 完整鉴定结果")
    
    all_columns = report.columns.tolist()
    default_cols = ['序号', '化合物中文名', '分子式', 'ppm', '评级名称', '药材名称', '综合得分']
    selected_cols = st.multiselect("选择显示的列", all_columns, default=[c for c in default_cols if c in all_columns])
    
    display_df = report[selected_cols] if selected_cols else report
    
    page_size = st.number_input("每页显示行数", min_value=10, max_value=100, value=20)
    page_num = st.number_input("当前页码", min_value=1, max_value=len(display_df)//page_size + 1, value=1)
    start_idx = (page_num - 1) * page_size
    end_idx = min(start_idx + page_size, len(display_df))
    
    st.dataframe(
        display_df.iloc[start_idx:end_idx],
        use_container_width=True,
        hide_index=True,
        column_config={
            'ppm': st.column_config.NumberColumn("ppm误差", format="%.4f"),
            '综合得分': st.column_config.NumberColumn("综合得分", format="%.2f")
        }
    )
    
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
    """主函数"""
    load_css()
    
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
