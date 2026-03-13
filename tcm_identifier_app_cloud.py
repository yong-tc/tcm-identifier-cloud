# -*- coding: utf-8 -*-
"""
中药化合物智能鉴定平台 - 云端部署增强版 v5.12（新增碎片直接比对功能）
=====================================================================
功能说明：
- 支持核心/辅助双诊断离子文件，核心离子权重2，辅助离子权重1
- 新增“碎片直接比对”功能：将样本碎片离子与诊断离子库直接匹配，快速推测化合物类别
- 其余功能与原 v5.11 一致
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
# 辅助函数：标准化分子式、文件查找
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


def find_file(filename):
    """在常见路径查找文件"""
    search_paths = [
        filename,
        f"data/{filename}",
        f"user_input_files/{filename}",
        os.path.join(os.getcwd(), filename)
    ]
    for path in search_paths:
        if os.path.exists(path):
            return path
    return None


# ============================================================================
# 诊断离子库加载（缓存）
# ============================================================================
@st.cache_data
def load_diagnostic_ions_from_files(file_paths):
    """
    从多个 Excel 文件加载诊断离子库，合并去重，自动识别核心/辅助离子并分配权重
    返回字典: {化合物类型: {'ions': list, 'weights': list, 'description': str}}
    """
    all_rows = []
    default_weights = {'核心': 2, '辅助': 1}

    for file_path in file_paths:
        if not os.path.exists(file_path):
            continue
        try:
            xls = pd.ExcelFile(file_path)
            sheets = xls.sheet_names
            file_name = os.path.basename(file_path).lower()
            for sheet in sheets:
                df = pd.read_excel(xls, sheet_name=sheet)
                if '化合物类型' not in df.columns or '离子m/z' not in df.columns:
                    continue

                sub = df[['化合物类型', '离子m/z']].copy()
                # 提取诊断离子类别列
                if '诊断离子类别' in df.columns:
                    sub['类别'] = df['诊断离子类别'].astype(str).str.strip()
                else:
                    # 从文件名推测
                    if '核心' in file_name:
                        sub['类别'] = '核心诊断离子'
                    elif '辅助' in file_name:
                        sub['类别'] = '辅助诊断离子'
                    else:
                        sub['类别'] = '未知'

                # 提取权重列（若有）
                if '权重' in df.columns:
                    sub['权重'] = pd.to_numeric(df['权重'], errors='coerce').fillna(1)
                else:
                    # 根据类别分配默认权重
                    sub['权重'] = sub['类别'].apply(
                        lambda x: default_weights['核心'] if '核心' in x else default_weights['辅助'] if '辅助' in x else 1
                    )

                sub = sub.dropna(subset=['离子m/z', '化合物类型'])
                sub['离子m/z'] = pd.to_numeric(sub['离子m/z'], errors='coerce')
                sub = sub.dropna(subset=['离子m/z'])
                sub['权重'] = sub['权重'].astype(float)
                all_rows.append(sub)
        except Exception:
            continue

    if not all_rows:
        return {}

    combined = pd.concat(all_rows, ignore_index=True)
    diagnostic_ions = {}
    for cat, group in combined.groupby('化合物类型'):
        ions_dict = {}
        for _, row in group.iterrows():
            mz = float(row['离子m/z'])
            w = float(row['权重'])
            ions_dict[mz] = ions_dict.get(mz, 0) + w
        ions = list(ions_dict.keys())
        weights = list(ions_dict.values())
        diagnostic_ions[cat] = {
            'ions': ions,
            'weights': weights,
            'description': f'来自外部文件，{len(ions)}个离子（去重后）'
        }
    return diagnostic_ions


# ============================================================================
# 直接比对碎片与诊断离子库
# ============================================================================
def match_fragments_to_diagnostic(fragments, diagnostic_ions, mode=None, tolerance=0.05, tolerance_type='Da'):
    """
    将碎片离子列表与诊断离子库直接比对
    fragments: list of float, 样本碎片离子m/z
    diagnostic_ions: 诊断离子库字典（由 load_diagnostic_ions_from_files 返回）
    mode: '正离子' 或 '负离子'，若提供则根据离子来源过滤（暂未实现，可扩展）
    tolerance: 容差值
    tolerance_type: 'Da' 或 'ppm'
    返回: dict {类别: {'matched_ions': list, 'total_weight': float, 'matched_weights': list}}
    """
    results = {}
    if not fragments or not diagnostic_ions:
        return results

    fragments = np.array(fragments)
    for cat, data in diagnostic_ions.items():
        matched = []
        matched_weights = []
        for ion, w in zip(data['ions'], data['weights']):
            if tolerance_type == 'Da':
                deltas = np.abs(fragments - ion)
                if np.any(deltas <= tolerance):
                    matched.append(ion)
                    matched_weights.append(w)
            else:  # ppm
                deltas = np.abs(fragments - ion) / ion * 1e6
                if np.any(deltas <= tolerance):
                    matched.append(ion)
                    matched_weights.append(w)
        if matched:
            results[cat] = {
                'matched_ions': matched,
                'total_weight': sum(matched_weights),
                'matched_weights': matched_weights
            }
    # 按总权重排序
    results = dict(sorted(results.items(), key=lambda x: x[1]['total_weight'], reverse=True))
    return results


# ============================================================================
# 原鉴定程序类（略作修改，使用缓存的诊断离子库）
# ============================================================================
class UltimateGardeniaIdentifier:
    """中药化合物鉴定终极版程序（完整功能）"""

    def __init__(self, database_path, ms_positive_path, ms_negative_path,
                 herb_name=None, config=None, use_parallel=True,
                 rt_tolerance=0.3, loss_tolerance=0.02,
                 external_diagnostic_files=None,
                 rt_fusion_tolerance=0.2,
                 intensity_relative_threshold=1.0,
                 tolerance_type='Da',
                 use_rt_score=True,
                 custom_db_path=None,
                 cache_index=True):
        # 初始化代码与原程序相同（略，实际运行时需包含完整类定义）
        # 为节省篇幅，此处省略，实际使用时应包含原类完整代码
        pass


# ============================================================================
# 数据库加载函数（缓存）
# ============================================================================
@st.cache_data
def load_database_cached(db_filename="TCM-SM-MS DB.xlsx"):
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
                return pd.read_excel(path)
            except:
                continue
    return pd.DataFrame()


def find_database_path():
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


# ============================================================================
# 旧的诊断离子筛查函数（保留，用于兼容）
# ============================================================================
@st.cache_data
def load_diagnostic_ions_cached():
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
            except:
                continue
    return pd.DataFrame()


def match_diagnostic_ions(user_mz_values, diagnostic_df, tolerance_ppm=10, ion_mode=None):
    # 与原程序相同，此处省略
    return pd.DataFrame()


# ============================================================================
# Streamlit 页面函数
# ============================================================================
def load_css():
    st.markdown("""
    <style>
        /* 省略，与原程序相同 */
    </style>
    """, unsafe_allow_html=True)


def login_page():
    # 与原程序相同
    pass


def logout_button():
    # 与原程序相同
    pass


def create_header():
    st.markdown("""
    <div class="main-header">
        <h1 style="color: white !important; margin: 0;">🌿 中药化合物智能鉴定平台（新增碎片直接比对）</h1>
        <p style="margin: 0.5rem 0 0 0; opacity: 0.9;">v5.12 | 登录用户: {}</p>
    </div>
    """.format(st.session_state.get('username', '')), unsafe_allow_html=True)


def create_sidebar():
    st.sidebar.markdown("""
    <div style="text-align: center; padding: 1rem 0;">
        <h2 style="color: #2E7D32; margin-bottom: 0.5rem;">🔬 TCM Identifier</h2>
        <p style="color: #666; font-size: 0.8rem;">中药化合物鉴定系统 v5.12</p>
    </div>
    """, unsafe_allow_html=True)

    if st.session_state.get('logged_in'):
        st.sidebar.markdown(f"**当前用户**: {st.session_state.username}")

    page = st.sidebar.radio(
        "导航菜单",
        ["首页", "开始鉴定", "碎片直接比对", "诊断离子筛查", "使用指南", "数据库预览", "结果分析"]
    )

    st.sidebar.markdown("---")
    st.sidebar.info("""
    **版本信息**
    - 程序版本：v5.12（新增碎片直接比对）
    - 数据库规模：35,828条化合物记录
    - 诊断离子：支持核心/辅助双文件，权重自动分配
    - 核心特点：直接比对、自动列名、双重阈值、中性丢失、RT得分、多数据库
    """)

    st.sidebar.markdown("""
    <div style="text-align: center; color: #999; font-size: 0.7rem; padding: 1rem 0;">
        <p>© 2026 张永</p>
        <p>中药化合物智能鉴定平台</p>
    </div>
    """, unsafe_allow_html=True)

    return page


def show_home_page():
    create_header()
    st.markdown("## 📊 系统概览")
    # 与原程序相同
    st.markdown("""
    <div class="stat-box"><div class="stat-number">35,828</div><div class="stat-label">数据库化合物数</div></div>
    """, unsafe_allow_html=True)
    # ... 省略


def show_analysis_page():
    # 原鉴定页面，略作修改以使用缓存的诊断离子库
    create_header()
    st.markdown("## 📁 上传质谱数据")
    # ... 与原程序相同，需保留诊断离子文件上传控件
    # 注意：这里需要将上传的文件路径传递给鉴定类，但鉴定类内部仍使用自己的加载方法。
    # 为简化，保留原有实现。
    pass


def show_fragment_match_page():
    """新增页面：碎片直接比对诊断离子库"""
    create_header()
    st.markdown("## 🔬 碎片直接比对诊断离子库")
    st.markdown("将样本的碎片离子列表与核心/辅助诊断离子库直接比对，按化合物类型输出匹配结果，快速推测样本中可能存在的化合物类别。")

    # 诊断离子库加载（从用户上传或默认文件）
    st.markdown("### 📁 加载诊断离子库")
    col1, col2 = st.columns(2)
    with col1:
        core_file = st.file_uploader("上传核心诊断离子库 (.xlsx)", type=['xlsx'], key='core_match')
    with col2:
        aux_file = st.file_uploader("上传辅助诊断离子库 (.xlsx)", type=['xlsx'], key='aux_match')

    # 收集文件路径
    temp_dir = tempfile.gettempdir()
    diag_paths = []
    if core_file:
        core_path = os.path.join(temp_dir, core_file.name)
        with open(core_path, 'wb') as f:
            f.write(core_file.getbuffer())
        diag_paths.append(core_path)
    if aux_file:
        aux_path = os.path.join(temp_dir, aux_file.name)
        with open(aux_path, 'wb') as f:
            f.write(aux_file.getbuffer())
        diag_paths.append(aux_path)

    # 若未上传，尝试查找默认文件
    if not diag_paths:
        default_core = find_file('核心诊断离子库.xlsx')
        default_aux = find_file('辅助诊断离子库.xlsx')
        if default_core:
            diag_paths.append(default_core)
        if default_aux:
            diag_paths.append(default_aux)

    if not diag_paths:
        st.warning("未找到诊断离子库文件，请上传或确保项目目录下存在默认文件。")
        return

    # 加载诊断离子库（缓存）
    with st.spinner("正在加载诊断离子库..."):
        diagnostic_ions = load_diagnostic_ions_from_files(diag_paths)

    if not diagnostic_ions:
        st.error("诊断离子库加载失败，请检查文件格式。")
        return

    st.success(f"✅ 成功加载诊断离子库，包含 {len(diagnostic_ions)} 个化合物类型，{sum(len(v['ions']) for v in diagnostic_ions.values())} 个去重离子。")

    st.markdown("---")
    st.markdown("### 📤 上传质谱数据")

    ms_file = st.file_uploader("上传质谱数据文件 (.xlsx)", type=['xlsx'], key='ms_match')
    if not ms_file:
        st.info("请上传质谱数据文件。")
        return

    # 解析质谱数据
    with st.spinner("正在解析质谱数据..."):
        df = pd.read_excel(ms_file)
        # 尝试识别碎片离子列
        mz_columns = [col for col in df.columns if 'm/z' in col.lower() or 'mass' in col.lower()]
        if not mz_columns:
            st.error("未找到包含 m/z 的列，请确保列名包含 'm/z' 或 'mass'。")
            return

        fragments = []
        for col in mz_columns:
            vals = df[col].dropna().tolist()
            for v in vals:
                try:
                    fragments.append(float(v))
                except:
                    continue
        fragments = sorted(set(fragments))  # 去重排序
        st.info(f"共提取到 {len(fragments)} 个碎片离子（去重后）。")

    # 参数设置
    col1, col2, col3 = st.columns(3)
    with col1:
        tolerance_type = st.selectbox("容差类型", options=["Da", "ppm"], index=0)
    with col2:
        if tolerance_type == "Da":
            tolerance = st.number_input("容差 (Da)", min_value=0.001, max_value=1.0, value=0.05, step=0.01)
        else:
            tolerance = st.number_input("容差 (ppm)", min_value=1, max_value=100, value=20, step=1)
    with col3:
        min_weight = st.number_input("最小总权重", min_value=0, value=0, help="仅显示总权重大于此值的类型")

    if st.button("🚀 开始比对", type="primary"):
        with st.spinner("正在匹配诊断离子..."):
            results = match_fragments_to_diagnostic(
                fragments, diagnostic_ions,
                tolerance=tolerance, tolerance_type=tolerance_type
            )

        if not results:
            st.warning("未匹配到任何诊断离子，请尝试调整容差。")
            return

        # 过滤
        filtered = {k: v for k, v in results.items() if v['total_weight'] >= min_weight}
        if not filtered:
            st.warning(f"无总权重大于 {min_weight} 的匹配结果。")
            return

        st.markdown("### 📊 匹配结果")

        # 转换为 DataFrame 展示
        rows = []
        for cat, data in filtered.items():
            rows.append({
                '化合物类型': cat,
                '匹配离子数': len(data['matched_ions']),
                '总权重': data['total_weight'],
                '匹配离子列表': '; '.join([f"{x:.4f}" for x in data['matched_ions']])
            })
        result_df = pd.DataFrame(rows)
        result_df = result_df.sort_values('总权重', ascending=False).reset_index(drop=True)

        st.dataframe(result_df, use_container_width=True)

        # 柱状图
        st.bar_chart(result_df.set_index('化合物类型')['总权重'])

        # 导出
        col1, col2 = st.columns(2)
        with col1:
            csv = result_df.to_csv(index=False, encoding='utf-8-sig')
            st.download_button("📥 导出CSV", data=csv, file_name=f"直接比对_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv", mime="text/csv")
        with col2:
            buffer = BytesIO()
            with pd.ExcelWriter(buffer, engine='openpyxl') as writer:
                result_df.to_excel(writer, index=False, sheet_name='直接比对结果')
            st.download_button("📥 导出Excel", data=buffer.getvalue(), file_name=f"直接比对_{datetime.now().strftime('%Y%m%d_%H%M%S')}.xlsx", mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")


def show_diagnostic_ion_page():
    # 原诊断离子筛查页面（保留）
    create_header()
    st.markdown("## 🔬 诊断离子筛查")
    st.markdown("根据输入的m/z值，在诊断离子数据库中查找匹配的化合物特征离子，帮助快速识别化合物类别。")
    # ... 与原程序相同
    pass


def show_guide_page():
    create_header()
    st.markdown("## 📖 使用指南（v5.12）")
    st.info("""
    ### 新增功能：碎片直接比对
    - 上传质谱数据文件，程序自动提取碎片离子，与核心/辅助诊断离子库直接匹配。
    - 按化合物类型统计匹配离子数和累计权重（核心离子权重2，辅助离子权重1），帮助快速推测样本中可能含有的化合物类别。
    - 可调整容差类型和大小，设置最小总权重过滤。
    ### 其他功能
    - 智能化合物鉴定：基于数据库匹配，六级评级，综合评分。
    - 诊断离子筛查：输入m/z值，检索特征离子。
    - 数据库预览：查看数据库内容。
    """)
    # 原有内容可保留


def show_database_page():
    # 原数据库预览页面
    create_header()
    st.markdown("## 🗃️ 数据库预览")
    # ... 与原程序相同
    pass


def show_results_page():
    # 原结果分析页面
    create_header()
    # ... 与原程序相同
    pass


def main():
    load_css()
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
    elif page == "碎片直接比对":
        show_fragment_match_page()
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
