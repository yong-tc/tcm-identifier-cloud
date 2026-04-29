# -*- coding: utf-8 -*-
"""
中药化合物鉴定报告筛选工具 v2.0
支持多条件组合筛选、智能排序、批量导出
"""

import pandas as pd
import os
import re
from datetime import datetime


def load_report(report_path):
    """加载鉴定报告"""
    if not os.path.exists(report_path):
        print(f"文件不存在: {report_path}")
        return None

    try:
        xlsx = pd.ExcelFile(report_path)
        print(f"发现 {len(xlsx.sheet_names)} 个工作表: {xlsx.sheet_names}")

        all_dfs = []
        for sheet in xlsx.sheet_names:
            df = pd.read_excel(xlsx, sheet_name=sheet)
            df['_来源Sheet'] = sheet
            all_dfs.append(df)

        combined_df = pd.concat(all_dfs, ignore_index=True)
        print(f"共加载 {len(combined_df)} 条记录")
        return combined_df
    except Exception as e:
        print(f"加载失败: {e}")
        return None


def quick_filter_presets(df):
    """快速筛选预设"""
    print("\n" + "="*60)
    print("快速筛选预设")
    print("="*60)
    presets = {
        '1': ('高置信模式', '确证级+高置信级, ppm≤20, 有碎片'),
        '2': ('中等置信', '推定级+, ppm≤50'),
        '3': ('宽松模式', 'ppm≤100, 有碎片'),
        '4': ('仅确证级', '确证级'),
        '5': ('仅黄酮类', 'ppm≤30, 黄酮'),
        '6': ('仅生物碱', 'ppm≤30, 生物碱'),
        '7': ('仅萜类', 'ppm≤30, 萜'),
        '8': ('确证+高置信无碎片', '确证级+高置信级, 无碎片'),
        '9': ('高质量无碎片', 'ppm≤10, 无碎片'),
        '10': ('综合评分高', '得分≥85'),
    }

    for num, (name, desc) in presets.items():
        print(f" {num}. {name}")
        print(f"    {desc}")
    print(" 0. 返回手动筛选")

    choice = input("\n选择预设 (0-10): ").strip()

    if choice == '1':
        filtered = df[(df['评级名称'].isin(['确证级', '高置信级'])) & (df['ppm'] <= 20) & (df['匹配碎片数'] > 0)].copy()
    elif choice == '2':
        filtered = df[(df['评级名称'].isin(['确证级', '高置信级', '推定级'])) & (df['ppm'] <= 50)].copy()
    elif choice == '3':
        filtered = df[(df['ppm'] <= 100) & (df['匹配碎片数'] > 0)].copy()
    elif choice == '4':
        filtered = df[df['评级名称'] == '确证级'].copy()
    elif choice == '5':
        filtered = df[(df['ppm'] <= 30) & (df['化合物中文名'].str.contains('黄酮|flavonoid', na=False, case=False))].copy()
    elif choice == '6':
        filtered = df[(df['ppm'] <= 30) & (df['化合物中文名'].str.contains('生物碱|alkaloid', na=False, case=False))].copy()
    elif choice == '7':
        filtered = df[(df['ppm'] <= 30) & (df['化合物中文名'].str.contains('萜|terpene', na=False, case=False))].copy()
    elif choice == '8':
        filtered = df[(df['评级名称'].isin(['确证级', '高置信级'])) & (df['匹配碎片数'] == 0)].copy()
    elif choice == '9':
        filtered = df[(df['ppm'] <= 10) & (df['匹配碎片数'] == 0)].copy()
    elif choice == '10':
        filtered = df[df['综合得分'] >= 85].copy()
    else:
        return None

    if choice in presets:
        print(f"\n应用: {presets[choice][0]}")
    return filtered


def create_filter_conditions(df):
    """创建筛选条件字典"""
    conditions = {}

    print("\n" + "="*60)
    print("多条件组合筛选器")
    print("="*60)

    # 1. 评级筛选
    if '评级名称' in df.columns:
        levels = df['评级名称'].dropna().unique().tolist()
        levels.sort()
        print("\n【1】评级筛选 (多选用逗号分隔，直接回车跳过):")
        print(f"   可选: {levels}")
        levels_input = input("   选择评级: ").strip()
        if levels_input:
            selected = [l.strip() for l in levels_input.split(',')]
            conditions['评级名称'] = selected

    # 2. ppm范围筛选
    if 'ppm' in df.columns:
        print("\n【2】ppm范围筛选 (直接回车跳过):")
        min_ppm = input("   最小ppm (如5): ").strip()
        max_ppm = input("   最大ppm (如50): ").strip()
        if min_ppm or max_ppm:
            conditions['ppm范围'] = (
                float(min_ppm) if min_ppm else -float('inf'),
                float(max_ppm) if max_ppm else float('inf')
            )

    # 3. 综合得分筛选
    if '综合得分' in df.columns:
        print("\n【3】综合得分筛选 (直接回车跳过):")
        min_score = input("   最低得分 (如70): ").strip()
        max_score = input("   最高得分 (如100): ").strip()
        if min_score or max_score:
            conditions['得分范围'] = (
                float(min_score) if min_score else 0,
                float(max_score) if max_score else 100
            )

    # 4. 匹配碎片数筛选
    if '匹配碎片数' in df.columns:
        print("\n【4】匹配碎片数筛选:")
        min_frags = input("   最少碎片数: ").strip()
        max_frags = input("   最多碎片数: ").strip()
        if min_frags or max_frags:
            conditions['碎片数范围'] = (
                int(min_frags) if min_frags else 0,
                int(max_frags) if max_frags else 999
            )
        print("   5. 仅显示有碎片")
        print("   6. 仅显示无碎片")
        frag_opt = input("   选择 (5/6跳过): ").strip()
        if frag_opt == '5':
            conditions['has_fragment'] = True
        elif frag_opt == '6':
            conditions['has_fragment'] = False

    # 5. 化合物名称关键词筛选
    print("\n【5】化合物名称关键词筛选:")
    include_kw = input("   包含关键词 (如: 苷、黄酮): ").strip()
    exclude_kw = input("   排除关键词 (如: 酸、酯): ").strip()
    if include_kw:
        conditions['包含关键词'] = [k.strip() for k in include_kw.split(',')]
    if exclude_kw:
        conditions['排除关键词'] = [k.strip() for k in exclude_kw.split(',')]

    # 6. 分子式筛选
    print("\n【6】分子式筛选:")
    formula_pattern = input("   包含元素 (如: C6H8O): ").strip()
    if formula_pattern:
        conditions['分子式包含'] = formula_pattern

    # 7. 化合物类型筛选
    if '化合物类型' in df.columns:
        types = df['化合物类型'].dropna().unique().tolist()
        types = [t for t in types if str(t) != 'nan']
        types.sort()
        print("\n【7】化合物类型筛选 (直接回车跳过):")
        for i, t in enumerate(types[:15], 1):
            print(f"   {i}. {t}")
        type_input = input(f"   选择类型编号 (1-{min(15,len(types))}): ").strip()
        if type_input:
            try:
                idx = int(type_input) - 1
                if 0 <= idx < len(types):
                    conditions['化合物类型'] = types[idx]
            except:
                conditions['化合物类型'] = type_input

    # 8. 分子量范围筛选
    if '匹配质量数' in df.columns:
        print("\n【8】分子量范围筛选:")
        min_mw = input("   最小分子量 (如100): ").strip()
        max_mw = input("   最大分子量 (如1000): ").strip()
        if min_mw or max_mw:
            conditions['分子量范围'] = (
                float(min_mw) if min_mw else 0,
                float(max_mw) if max_mw else 100000
            )

    # 9. RT保留时间筛选
    if '可能出峰时间' in df.columns:
        print("\n【9】保留时间范围筛选:")
        print("   提示: 可能出峰时间格式为逗号分隔的时间点")
        min_rt = input("   最早时间(分钟, 如5): ").strip()
        max_rt = input("   最晚时间(分钟, 如30): ").strip()
        if min_rt or max_rt:
            conditions['RT范围'] = (
                float(min_rt) if min_rt else 0,
                float(max_rt) if max_rt else 999
            )

    # 10. 数据库来源筛选
    if '数据来源' in df.columns:
        sources = df['数据来源'].dropna().unique().tolist()
        sources.sort()
        print("\n【10】数据库来源筛选 (直接回车跳过):")
        print(f"   可选: {sources}")
        source = input("   选择来源: ").strip()
        if source:
            conditions['数据来源'] = source

    # 11. 一级碎片列表筛选
    print("\n【11】母离子m/z筛选:")
    min_mz = input("   最小m/z (如100): ").strip()
    max_mz = input("   最大m/z (如1000): ").strip()
    if min_mz or max_mz:
        conditions['m/z范围'] = (
            float(min_mz) if min_mz else 0,
            float(max_mz) if max_mz else 100000
        )

    return conditions


def apply_filters(df, conditions):
    """应用筛选条件"""
    if df is None or df.empty:
        return df

    filtered = df.copy()
    total_before = len(filtered)

    print("\n" + "-"*60)
    print("筛选过程:")
    print("-"*60)

    # 评级筛选
    if '评级名称' in conditions:
        levels = conditions['评级名称']
        filtered = filtered[filtered['评级名称'].isin(levels)]
        print(f"  评级筛选 ({levels}): {len(filtered)} 条")

    # ppm范围筛选
    if 'ppm范围' in conditions:
        min_ppm, max_ppm = conditions['ppm范围']
        filtered = filtered[(filtered['ppm'] >= min_ppm) & (filtered['ppm'] <= max_ppm)]
        print(f"  ppm范围筛选 ({min_ppm}-{max_ppm}): {len(filtered)} 条")

    # 得分范围筛选
    if '得分范围' in conditions:
        min_score, max_score = conditions['得分范围']
        filtered = filtered[(filtered['综合得分'] >= min_score) & (filtered['综合得分'] <= max_score)]
        print(f"  得分范围筛选 ({min_score}-{max_score}): {len(filtered)} 条")

    # 碎片数范围筛选
    if '碎片数范围' in conditions:
        min_f, max_f = conditions['碎片数范围']
        filtered = filtered[(filtered['匹配碎片数'] >= min_f) & (filtered['匹配碎片数'] <= max_f)]
        print(f"  碎片数筛选 ({min_f}-{max_f}): {len(filtered)} 条")

    # 有无碎片筛选
    if 'has_fragment' in conditions:
        if conditions['has_fragment']:
            filtered = filtered[filtered['匹配碎片数'] > 0]
        else:
            filtered = filtered[filtered['匹配碎片数'] == 0]
        print(f"  碎片状态: {'有碎片' if conditions['has_fragment'] else '无碎片'}: {len(filtered)} 条")

    # 包含关键词筛选
    if '包含关键词' in conditions:
        for kw in conditions['包含关键词']:
            mask = filtered['化合物中文名'].str.contains(kw, na=False, case=False) | \
                   filtered['化合物英文名'].str.contains(kw, na=False, case=False) if '化合物英文名' in filtered.columns else \
                   filtered['化合物中文名'].str.contains(kw, na=False, case=False)
            filtered = filtered[mask]
        print(f"  包含关键词{conditions['包含关键词']}: {len(filtered)} 条")

    # 排除关键词筛选
    if '排除关键词' in conditions:
        for kw in conditions['排除关键词']:
            mask = ~filtered['化合物中文名'].str.contains(kw, na=False, case=False) & \
                   ~filtered['化合物英文名'].str.contains(kw, na=False, case=False) if '化合物英文名' in filtered.columns else \
                   ~filtered['化合物中文名'].str.contains(kw, na=False, case=False)
            filtered = filtered[mask]
        print(f"  排除关键词{conditions['排除关键词']}: {len(filtered)} 条")

    # 分子式筛选
    if '分子式包含' in conditions:
        pattern = conditions['分子式包含']
        filtered = filtered[filtered['分子式'].str.contains(pattern, na=False, case=False)]
        print(f"  分子式包含'{pattern}': {len(filtered)} 条")

    # 化合物类型筛选
    if '化合物类型' in conditions:
        compound_type = conditions['化合物类型']
        filtered = filtered[filtered['化合物类型'].str.contains(compound_type, na=False, case=False)]
        print(f"  类型包含'{compound_type}': {len(filtered)} 条")

    # 分子量范围筛选
    if '分子量范围' in conditions:
        min_mw, max_mw = conditions['分子量范围']
        filtered = filtered[(filtered['匹配质量数'] >= min_mw) & (filtered['匹配质量数'] <= max_mw)]
        print(f"  分子量范围 ({min_mw}-{max_mw}): {len(filtered)} 条")

    # RT范围筛选
    if 'RT范围' in conditions:
        min_rt, max_rt = conditions['RT范围']
        def in_rt_range(rt_str):
            if pd.isna(rt_str):
                return False
            try:
                rts = [float(r.strip()) for r in str(rt_str).split(',')]
                return any(min_rt <= r <= max_rt for r in rts)
            except:
                return False
        filtered = filtered[filtered['可能出峰时间'].apply(in_rt_range)]
        print(f"  RT范围 ({min_rt}-{max_rt}min): {len(filtered)} 条")

    # 数据来源筛选
    if '数据来源' in conditions:
        source = conditions['数据来源']
        filtered = filtered[filtered['数据来源'] == source]
        print(f"  来源='{source}': {len(filtered)} 条")

    # m/z范围筛选
    if 'm/z范围' in conditions:
        min_mz, max_mz = conditions['m/z范围']
        # 从一级碎片列提取m/z
        def extract_first_mz(frag_str):
            if pd.isna(frag_str):
                return None
            try:
                parts = str(frag_str).split(',')
                if parts:
                    mz = float(parts[0].strip())
                    return mz
            except:
                pass
            return None
        if '一级碎片' in filtered.columns:
            filtered['_tmp_mz'] = filtered['一级碎片'].apply(extract_first_mz)
            filtered = filtered[(filtered['_tmp_mz'] >= min_mz) & (filtered['_tmp_mz'] <= max_mz)]
            filtered = filtered.drop('_tmp_mz', axis=1)
            print(f"  m/z范围 ({min_mz}-{max_mz}): {len(filtered)} 条")

    print("-"*60)
    print(f"总计: {total_before} → {len(filtered)} 条 (过滤 {total_before - len(filtered)} 条)")

    return filtered


def print_filter_summary(df, original_count=0):
    """打印筛选结果摘要"""
    if df is None or df.empty:
        return

    print("\n" + "="*60)
    print("筛选结果摘要")
    print("="*60)

    print(f"\n共筛选出 {len(df)} 个化合物")
    if original_count > 0:
        pct = len(df) / original_count * 100
        print(f"占原始数据的 {pct:.1f}%")

    # 评级分布
    if '评级名称' in df.columns:
        print("\n【评级分布】")
        for level, count in df['评级名称'].value_counts().items():
            pct = count / len(df) * 100
            bar = '█' * int(pct / 5)
            print(f"  {level:10s}: {count:5d} ({pct:5.1f}%) {bar}")

    # 碎片匹配分布
    if '匹配碎片数' in df.columns:
        has_frag = (df['匹配碎片数'] > 0).sum()
        no_frag = (df['匹配碎片数'] == 0).sum()
        print(f"\n【碎片匹配】")
        print(f"  有碎片: {has_frag} ({has_frag/len(df)*100:.1f}%)")
        print(f"  无碎片: {no_frag} ({no_frag/len(df)*100:.1f}%)")

    # 得分统计
    if '综合得分' in df.columns:
        print(f"\n【综合得分】")
        print(f"  最高: {df['综合得分'].max():.1f}")
        print(f"  平均: {df['综合得分'].mean():.1f}")
        print(f"  中位数: {df['综合得分'].median():.1f}")

    # ppm统计
    if 'ppm' in df.columns:
        print(f"\n【ppm误差分布】")
        bins = [0, 10, 20, 30, 50, 100, 999]
        labels = ['0-10', '10-20', '20-30', '30-50', '50-100', '>100']
        df['_ppm_bin'] = pd.cut(df['ppm'], bins=bins, labels=labels)
        for label in labels:
            count = (df['_ppm_bin'] == label).sum()
            pct = count / len(df) * 100
            bar = '█' * int(pct / 5)
            print(f"  {label:>8s}: {count:5d} ({pct:5.1f}%) {bar}")
        df.drop('_ppm_bin', axis=1, inplace=True)

    # 化合物类型分布
    if '化合物类型' in df.columns:
        types = df['化合物类型'].value_counts().head(8)
        print(f"\n【化合物类型 (Top 8)】")
        for t, count in types.items():
            print(f"  {t}: {count}")

    # 分子量分布
    if '匹配质量数' in df.columns:
        print(f"\n【分子量统计】")
        print(f"  范围: {df['匹配质量数'].min():.1f} - {df['匹配质量数'].max():.1f}")
        print(f"  平均: {df['匹配质量数'].mean():.1f}")


def display_top_results(df, n=20):
    """显示前N条结果"""
    if df is None or df.empty:
        return

    print("\n" + "="*60)
    print(f"前 {min(n, len(df))} 条结果预览 (按综合得分排序):")
    print("="*60)

    # 选择关键列显示
    display_cols = []
    for col in ['序号', '化合物中文名', '分子式', 'ppm', '综合得分', '评级名称', '匹配碎片数']:
        if col in df.columns:
            display_cols.append(col)

    if display_cols:
        if '综合得分' in df.columns:
            display_df = df.sort_values('综合得分', ascending=False).head(n)
        else:
            display_df = df.head(n)

        # 设置pandas显示选项
        pd.set_option('display.max_columns', None)
        pd.set_option('display.width', None)
        pd.set_option('display.max_colwidth', 25)
        print(display_df[display_cols].to_string(index=False))
        pd.reset_option('display.max_columns')
        pd.reset_option('display.width')
        pd.reset_option('display.max_colwidth')


def sort_results(df):
    """排序结果"""
    if df is None or df.empty:
        return df

    print("\n排序方式:")
    print(" 1. 综合得分 (高→低)")
    print(" 2. ppm误差 (低→高)")
    print(" 3. 匹配碎片数 (多→少)")
    print(" 4. 分子量 (小→大)")
    print(" 5. 分子量 (大→小)")
    print(" 6. 评级 (确证→提示)")
    print(" 0. 保持原序")

    choice = input("\n选择排序方式: ").strip()

    if choice == '1':
        return df.sort_values('综合得分', ascending=False)
    elif choice == '2':
        return df.sort_values('ppm', ascending=True)
    elif choice == '3':
        return df.sort_values('匹配碎片数', ascending=False)
    elif choice == '4':
        return df.sort_values('匹配质量数', ascending=True) if '匹配质量数' in df.columns else df
    elif choice == '5':
        return df.sort_values('匹配质量数', ascending=False) if '匹配质量数' in df.columns else df
    elif choice == '6':
        level_order = {'确证级': 1, '高置信级': 2, '推定级': 3, '提示级': 4}
        if '评级名称' in df.columns:
            df['_sort_level'] = df['评级名称'].map(level_order).fillna(99)
            df = df.sort_values('_sort_level')
            df.drop('_sort_level', axis=1, inplace=True)
        return df
    else:
        return df


def save_filtered_report(df, output_path=None):
    """保存筛选后的报告"""
    if df is None or df.empty:
        print("没有数据可保存")
        return None

    if output_path is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_path = f"筛选结果_{timestamp}.xlsx"

    # 分割成两个sheet
    with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
        if '匹配碎片数' in df.columns:
            with_frag = df[df['匹配碎片数'] > 0].copy()
            no_frag = df[df['匹配碎片数'] == 0].copy()

            if not with_frag.empty:
                with_frag.to_excel(writer, sheet_name='有碎片匹配', index=False)
                print(f"  有碎片匹配: {len(with_frag)} 条")
            if not no_frag.empty:
                no_frag.to_excel(writer, sheet_name='无碎片匹配', index=False)
                print(f"  无碎片匹配: {len(no_frag)} 条")
        else:
            df.to_excel(writer, sheet_name='筛选结果', index=False)

    print(f"\n报告已保存: {output_path}")
    return output_path


def export_summary_report(df, output_path=None):
    """导出精简摘要报告"""
    if df is None or df.empty:
        return None

    if output_path is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_path = f"摘要报告_{timestamp}.xlsx"

    # 选择关键列
    summary_cols = ['序号', '化合物中文名', '分子式', '匹配质量数', 'ppm',
                    '综合得分', '评级名称', '匹配碎片数', '化合物类型']
    available_cols = [c for c in summary_cols if c in df.columns]

    summary_df = df[available_cols].copy()

    if '综合得分' in df.columns:
        summary_df = summary_df.sort_values('综合得分', ascending=False)

    summary_df.to_excel(output_path, index=False)
    print(f"摘要报告已保存: {output_path}")
    return output_path


def interactive_filter(df):
    """交互式筛选主函数"""
    original_count = len(df)
    current_df = df

    while True:
        print("\n" + "="*60)
        print("中药化合物鉴定报告筛选工具 v2.0")
        print("="*60)
        print(f"当前数据: {len(current_df)} 条 (原始: {original_count})")
        print("\n主菜单:")
        print(" 1. 快速预设筛选")
        print(" 2. 自定义多条件筛选")
        print(" 3. 对结果排序")
        print(" 4. 显示统计摘要")
        print(" 5. 预览前20条")
        print(" 6. 导出完整报告")
        print(" 7. 导出精简摘要")
        print(" 8. 重置(回到原始数据)")
        print(" 9. 保存当前筛选结果")
        print(" 0. 退出")

        choice = input("\n选择操作: ").strip()

        if choice == '1':
            filtered = quick_filter_presets(df)
            if filtered is not None:
                print_filter_summary(filtered, original_count)
                display_top_results(filtered)
                current_df = filtered

        elif choice == '2':
            conditions = create_filter_conditions(df)
            if conditions:
                filtered = apply_filters(df, conditions)
                if not filtered.empty:
                    print_filter_summary(filtered, original_count)
                    display_top_results(filtered)
                    current_df = filtered
                else:
                    print("\n筛选结果为空!")
            else:
                print("未设置任何筛选条件")

        elif choice == '3':
            current_df = sort_results(current_df)
            display_top_results(current_df)

        elif choice == '4':
            print_filter_summary(current_df, original_count)

        elif choice == '5':
            display_top_results(current_df)

        elif choice == '6':
            save_filtered_report(current_df)

        elif choice == '7':
            export_summary_report(current_df)

        elif choice == '8':
            current_df = df
            print(f"\n已重置，回到原始数据: {len(current_df)} 条")

        elif choice == '9':
            save_filtered_report(current_df)

        elif choice == '0':
            print("\n感谢使用，再见!")
            break

        else:
            print("\n无效选择，请重试")


def main():
    """主函数"""
    print("="*60)
    print("中药化合物鉴定报告筛选工具 v2.0")
    print("="*60)

    # 列出当前目录的Excel文件
    excel_files = [f for f in os.listdir('.') if f.endswith('.xlsx') and not f.startswith('~')]
    if excel_files:
        print("\n当前目录Excel文件:")
        for i, f in enumerate(excel_files[:10], 1):
            print(f"  {i}. {f}")

    # 输入报告文件路径
    default_path = "user_input_files/栀子全库鉴定报告_20260428_170519.xlsx"
    print(f"\n默认文件: {default_path}")
    report_path = input("输入报告文件路径 (直接回车使用默认): ").strip()

    if not report_path:
        # 尝试使用目录中第一个xlsx文件
        if excel_files:
            report_path = excel_files[0]
            print(f"使用: {report_path}")
        else:
            report_path = default_path

    # 加载报告
    df = load_report(report_path)
    if df is None:
        return

    # 进入交互式筛选
    interactive_filter(df)


if __name__ == "__main__":
    main()
