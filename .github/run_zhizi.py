# -*- coding: utf-8 -*-
"""
栀子化合物鉴定 - 使用全库数据比对
"""
import pandas as pd
import numpy as np
import os
from datetime import datetime
from tcm_identifier import UltimateGardeniaIdentifier, find_database_path

def main():
    print("="*80)
    print("栀子化合物鉴定分析（全库比对）")
    print("="*80)

    # 质谱数据路径
    ms_positive_path = "user_input_files/栀子-MS+.xlsx"
    ms_negative_path = "user_input_files/栀子-MS-.xlsx"

    # 查找数据库
    db_path = find_database_path()
    if not db_path:
        print("错误: 未找到主数据库文件!")
        return

    print(f"数据库路径: {db_path}")

    # 创建鉴定程序实例 - 使用全库数据（herb_name=None）
    identifier = UltimateGardeniaIdentifier(
        database_path=db_path,
        ms_positive_path=ms_positive_path,
        ms_negative_path=ms_negative_path,
        herb_name=None,  # 使用全库数据进行比对
        config={
            'min_intensity': 100,
            'fragment_tolerance': 0.15,  # 固定Da容差
            'fragment_tolerance_ppm': 20,
            'tolerance_ppm': 50,
            'max_candidates': 5,
            'max_ppm': 100,
            'ppm_tier1': 10,
            'ppm_tier2': 20,
            'ppm_tier3': 50,
        },
        use_parallel=True,
        rt_tolerance=0.3,
        tolerance_type='Da',
        use_rt_score=True,
        custom_db_path=None,
        english_db_path=None,
        standard_db_path=None,
        cache_index=True,
        strict_mode=False,
        enable_adduct_expansion=True,
        max_ppm=100,
        enable_diagnostic_logging=True
    )

    # 生成报告
    report = identifier.generate_report('栀子')

    # 保存报告
    output_path = f"user_input_files/栀子全库鉴定报告_{datetime.now().strftime('%Y%m%d_%H%M%S')}.xlsx"
    identifier.save_report(report, output_path)

    # 打印摘要
    identifier.print_summary(report, "栀子（全库）")

    return report

if __name__ == "__main__":
    main()
