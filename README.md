# ChemReactionDrawer 化学二维分子与反应图生成器

Python + RDKit + Tkinter。绘图来自真实分子连接结构，不使用 AI 文生图。支持中文界面、PNG（300 DPI 元数据）与纯矢量 SVG、反应条件、整数计量系数、离子、明确指定的多组分结构，以及元素/同位素/净电荷守恒检查。

## Windows 从零运行

1. 从 https://www.python.org/downloads/windows/ 安装 **Python 3.12 64 位**。安装时保留 `pip`、`tcl/tk and IDLE`、Python Launcher，并勾选 `Add python.exe to PATH`。
2. 下载完整项目 ZIP，并解压到可写目录，例如 `D:\ChemReactionDrawer`。不要直接从压缩包内部启动。
3. 双击 `start.bat`。首次启动自动创建 `.venv` 并安装依赖，需要联网；后续依赖完整时直接启动。若失败，窗口保留错误信息。
4. 输入单个分子名称、SMILES 或反应式，填写可选条件/备注，点击“生成图片”。预览成功后，选择 PNG/SVG，点击“保存图片”。生成一次即可分别保存两种格式。
5. SVG 可用浏览器、Inkscape 或 Illustrator 打开；PNG 可用 Windows 照片查看。图形界面预览会缩放，保存文件保持完整分辨率。

手动安装、运行（在项目目录打开终端）：

```bat
py -3.12 -m venv .venv
.venv\Scripts\python.exe -m pip install -r requirements.txt
.venv\Scripts\python.exe chem_drawer.py
```

`tkinter` 随官方 Windows Python 安装，不是 pip 包。RDKit 已确认有本版本 Python 3.12 / Windows x64 二进制安装包。Windows 的双击流程仍需在 Windows 上验证；云端测试运行于 Linux。

## 文件结构

```text
ChemReactionDrawer/
├── chem_drawer.py           完整程序，GUI + 命令行
├── start.bat                Windows 双击入口
├── requirements.txt         运行依赖
├── README.md                本说明
├── TEST_REPORT.md           实际测试记录及限制
├── tests/
│   └── test_reactions.py     自动化测试
└── output/                  示例图片，运行命令生成，Git 忽略
    ├── example1.png / example1.svg
    ├── example2.png / example2.svg
    └── example3.png / example3.svg
```

`.venv/` 首次启动自动创建，不包含在分发包中。

## 支持的输入

### 单个分子（无需反应箭头）

默认“自动识别”模式：没有箭头即作为一个完整分子名称解析；有箭头按反应式解析。也可显式选择“单个分子”或“反应式”。单个分子不做反应守恒检查，状态栏显示分子式和结构来源；备注放在结构下方。

直接输入以下任意一个即可生成独立二维结构图：

```text
Nα-(2,4-dinitro-5-fluorophenyl)-L-alaninamide
N-alpha-(2,4-dinitro-5-fluorophenyl)-L-alaninamide
乙醇
SMILES:CCO
SMILES:[Na+].[Cl-]
```

用户示例已经实际通过 PubChem 查到唯一记录 **CID 5486955**，分子式 C9H9FN4O5，IUPAC 名称 `(2S)-2-(5-fluoro-2,4-dinitroanilino)propanamide`。程序内置从该记录核验的含立体信息 SMILES，原始名称和 `N-alpha` 写法都可离线使用；其他准确名称仍按内置表/PubChem 查询，不进行通用 IUPAC 名称猜测。数据来源：https://pubchem.ncbi.nlm.nih.gov/compound/5486955

CLI 示例（Windows）：

```bat
.venv\Scripts\python.exe chem_drawer.py --offline --input "Nα-(2,4-dinitro-5-fluorophenyl)-L-alaninamide" --output output/single-alaninamide
```

`--input` 和原来的 `--equation` 为同一个参数的别名；可用 `--mode molecule` 强制单分子模式。包含配位键箭头的 SMILES（例如 `SMILES:N->[Cu+2]`）请使用“单个分子”模式，避免自动判为反应。

### 普通方程式 / 化学名称

```text
CH3COOH + C2H5OH -> CH3COOC2H5 + H2O
2 H2 + O2 -> 2 H2O
乙酸 + 乙醇 -> 乙酸乙酯 + 水
acetic acid + ethanol -> ethyl acetate + water
H+ + OH- -> H2O
Na+ + Cl- -> NaCl
SMILES:[Na+].[Cl-] -> NaCl
```

- 物种分隔加号**两侧必须有空格**；电荷加号不加空格，例如 `Na+ + Cl-`。
- 系数**后必须有空格**，例如 `2 H2`；只接受 1–1000 的正整数。分数系数请整体乘倍数。
- 箭头支持 `->`、`→`、`=>`；不处理可逆箭头、半反应电子 `e-` 或未确定的聚合物。
- 先查询内置明确结构表；未收录名称可通过 PubChem 查询，建议使用准确英文名称。中文名称依赖 PubChem 是否收录该同义词。
- 未收录分子式一律拒绝解析。例如 `C2H6O` 可能是乙醇或二甲醚，不能选择一个结构。请使用准确名称或 `SMILES:CCO`。
- 内置 `C2H5OH`、`CH3COOH` 是明确缩写；内置简单物种如 H2、O2、H2O 是经过明确指定的结构。未知缩写不自动猜测。
- PubChem 返回多个 CID 时拒绝生成，要求更明确结构。单条名称记录也不能替用户确定未指定的立体化学；需要精确异构体时提供含立体信息的 SMILES。
- 盐、溶剂化物、混合物应通过 `SMILES:` 指定完整点分结构。不根据 `CuSO4`、水合盐分子式猜测配位或离子结构。

### Reaction SMILES

```text
CC(=O)O.CCO>>CCOC(C)=O.O
CC(=O)O.CCO>O=S(=O)(O)O>CCOC(C)=O.O
([Na+].[Cl-])>>([Na+].[Cl-])
```

中间段为试剂/催化剂，不计入反应物和生成物守恒，并在箭头上标注。普通 Reaction SMILES 的点表示不同反应条目；要保留同一个盐的多组分单元，使用上述括号分组。不同输入模式不可混用；Reaction SMILES 中计量系数通过重复条目表达，不能写 `2 CCO`。原子映射编号仅用于解析，绘图时去除。

## 守恒与绘图行为

程序通过 RDKit 补足隐式氢后统计元素、同位素及形式电荷，再乘以计量系数。任何差异都会列出并停止生成，避免保存错误图。无法确定结构、空物种、非法 SMILES、通配原子或网络失败也会明确报错。守恒检查不能验证反应机理、产率、条件、平衡或实际可行性；程序不会自动配平。

分子用 RDKit 计算二维坐标，独立面板排列；箭头、加号、系数单独绘制。PNG 使用 RDKit Cairo 绘图加 Pillow 合成；SVG 保留 RDKit 矢量路径。白底、黑色键线，PNG 写入 300 DPI。物理尺寸由像素数和 DPI 共同决定，请勿在排版中无限放大 PNG；需要放大请选 SVG。

水等单原子骨架物种由 RDKit 标为 `H₂O`；氢气显示 H–H、氧气显示 O=O。它们仍来自真实分子数据，并非把整个方程式转成文字。

反应条件支持最多 240 字符、4 行；分子条目最多 30 个，超大画布会拒绝生成。复杂大分子和极长试剂标签可能需简化或拆分反应。Windows PNG 中文字体使用微软雅黑；Linux 优先使用 Noto CJK，无中文字体时请安装字体，否则中文条件可能显示缺字。SVG 中文依赖打开它的系统字体。

## 命令行 / 测试

```bat
.venv\Scripts\python.exe chem_drawer.py --offline --equation "CH3COOH + C2H5OH -> CH3COOC2H5 + H2O" --conditions "H2SO4 / 80 C" --output output/example1
.venv\Scripts\python.exe chem_drawer.py --offline --equation "2 H2 + O2 -> 2 H2O" --output output/example2
.venv\Scripts\python.exe chem_drawer.py --offline --equation "CC(=O)O.CCO>>CCOC(C)=O.O" --output output/example3
.venv\Scripts\python.exe -m pip install pytest
.venv\Scripts\python.exe -m pytest -q
```

命令行一次写出两种格式，解析失败退出码为 2。Linux/macOS 同样用 Python 虚拟环境安装 `requirements.txt` 后运行；Linux GUI 还需要系统 Tk 与图形显示环境。

## PubChem 网络要求

在线名称查询需要访问 `https://pubchem.ncbi.nlm.nih.gov`，请求保留 TLS 验证，超时会报错，不会静默替换结构。三个指定示例及本地名称无需网络。初次开发时云端请求被拒绝；本次更新已经成功实际查询并核验用户给出的单分子名称及结构。其他名称仍取决于数据库收录情况和网络可用性；如使用受限网络，请在环境设置允许该域名。网络不可用时不影响本地映射、SMILES、守恒检查与绘图。
