表示碎片151.003出现在第0、1篇文献，碎片137.024出现在第1、2篇文献。
文献索引对应“文献来源”列中分号分隔的顺序（从0开始）。
若不指定索引（如`151.003`），则默认属于所有文献。

### 常见问题
- **鉴定结果为空**：检查ppm容差是否过小，或药材筛选是否正确。
- **缓存加载失败**：删除 `index_cache.pkl` 文件后重试。
- **文献匹配数为0**：检查数据库是否配置了带索引的碎片字段。
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
