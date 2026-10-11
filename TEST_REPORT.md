# 实际测试报告

运行环境：Linux，Python 3.12.14，RDKit 2026.9.1，Pillow 12.3.0。

- `python -m pytest -q`：21 项通过；测试覆盖三个指定反应式、PNG 像素与 DPI、SVG 分子路径、歧义拒绝、计量系数、离子电荷、盐的多组分分组、同位素不守恒、中文名称、箭头试剂标注、模拟 PubChem 单条/多条记录及立体信息保留。
- 三个指定示例实际通过命令行各生成 PNG 和 SVG，共六个文件，位于 `output/`。
- PNG 使用 Pillow 完整解码，写入 DPI 约 299.9994（PNG pHYs 整数换算造成误差，对应标准 300 DPI）。实际查看示例 1、2，真实骨架/键线清楚，白底，箭头及计量系数不重叠。
- 三个 SVG 均经过 XML 解析，并由独立 CairoSVG 渲染器成功转成 PNG；实际查看示例 1 的矢量渲染结果正常。CairoSVG 仅用于开发验证，不是运行依赖。
- Tkinter 在 Xvfb 虚拟显示中实际启动，自动调用生成按钮、检查预览 Canvas、分别保存 PNG/SVG、检查歧义错误后禁用保存且清空旧图。`DISPLAY=:99 python tests/gui_smoke.py` 通过。保存文件可打开。
- 当前云端 Tk 的字体列表不含中文字体，中文控件在此虚拟显示环境缺字；Windows 官方 Python 的中文界面显示尚未实机验证。Pillow 绘制中文条件使用系统 Noto CJK，已测试生成。
- PubChem 真实在线查询尝试 `acetone`，当前云端返回 ProxyError；因此不声称真实在线查询成功。查询处理逻辑用模拟 API 响应测试通过。受限网络应允许 `pubchem.ncbi.nlm.nih.gov`。
- PyPI 上 `rdkit==2026.9.1` 的 Python 3.12 Windows x64 wheel 已实际下载确认存在；没有在 Windows 安装或执行 `.bat`，双击流程未经 Windows 实机测试。
- 浏览器 Headless 尝试未完成，已停止；SVG 可打开的结论来自实际独立渲染器验证，不来自该浏览器尝试。

程序拒绝不守恒反应，不自动配平；准确结构及守恒不代表化学反应本身可行。
