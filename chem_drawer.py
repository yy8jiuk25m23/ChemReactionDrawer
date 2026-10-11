"""真实分子数据驱动的二维分子/反应图生成器；GUI 或命令行均可使用。"""
from __future__ import annotations
import argparse
import io
import json
import math
import re
import threading
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from xml.sax.saxutils import escape

import requests
from PIL import Image, ImageDraw, ImageFont
from rdkit import Chem
from rdkit.Geometry import Point2D
from rdkit.Chem import rdDepictor, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdChemReactions

class ReactionError(ValueError):
    pass

# 此表仅含明确结构及常用、可确定的简单物种，分子式不进行在线“第一条结果”猜测。
STRUCTURES = {
    'CH3COOH': 'CC(=O)O', 'C2H5OH': 'CCO', 'CH3CH2OH': 'CCO',
    'CH3COOC2H5': 'CC(=O)OCC', 'H2O': 'O', 'H2': '[H][H]', 'O2': 'O=O',
    'N2': 'N#N', 'CO2': 'O=C=O', 'NH3': 'N', 'CH4': 'C',
    'H+': '[H+]', 'Na+': '[Na+]', 'K+': '[K+]', 'Cl-': '[Cl-]',
    'OH-': '[OH-]', 'Ca2+': '[Ca+2]', 'Mg2+': '[Mg+2]', 'Fe2+': '[Fe+2]',
    'Fe3+': '[Fe+3]', 'NH4+': '[NH4+]', 'SO4^2-': '[O-]S(=O)(=O)[O-]',
    'NaCl': '[Na+].[Cl-]', 'HCl': 'Cl', 'NaOH': '[Na+].[OH-]',
    '乙酸': 'CC(=O)O', '醋酸': 'CC(=O)O', '乙醇': 'CCO', '乙酸乙酯': 'CC(=O)OCC',
    '水': 'O', '氢气': '[H][H]', '氧气': 'O=O',
    'acetic acid': 'CC(=O)O', 'ethanol': 'CCO', 'ethyl acetate': 'CC(=O)OCC',
    'water': 'O', 'hydrogen': '[H][H]', 'oxygen': 'O=O',
}
FORMULA = re.compile(r'(?:[A-Z][a-z]?\d*|[()\[\]·.^\d+−-])+\Z')

# PUG REST 实际核验：CID 5486955，C9H9FN4O5，(2S) 构型。
# 来源：https://pubchem.ncbi.nlm.nih.gov/compound/5486955
VERIFIED_NAMES = {
    name.lower(): ('C[C@@H](C(=O)N)NC1=CC(=C(C=C1[N+](=O)[O-])[N+](=O)[O-])F', 5486955)
    for name in (
        'Nα-(2,4-dinitro-5-fluorophenyl)-L-alaninamide',
        'N-alpha-(2,4-dinitro-5-fluorophenyl)-L-alaninamide',
    )
}

@dataclass
class Species:
    mol: Chem.Mol
    label: str
    coefficient: int = 1
    source: str = '本地明确结构'

@dataclass
class Reaction:
    reactants: list[Species]
    products: list[Species]
    agents: list[Species]
    balance: str


def parse_mol(smiles: str) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or mol.GetNumAtoms() == 0:
        raise ReactionError(f'无效或空的 SMILES：{smiles}')
    if any(a.GetAtomicNum() == 0 for a in mol.GetAtoms()):
        raise ReactionError('含通配符或未确定原子的结构无法绘制、验证。')
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return mol


def resolve(text: str, online: bool = True) -> Species:
    key = text.strip()
    if key.lower() in VERIFIED_NAMES:
        smiles, cid = VERIFIED_NAMES[key.lower()]
        return Species(parse_mol(smiles), key, source=f'PubChem CID {cid}（已核验本地结构）')
    smiles = STRUCTURES.get(key, STRUCTURES.get(key.lower()))
    if smiles:
        return Species(parse_mol(smiles), key)
    if key.startswith('SMILES:'):
        return Species(parse_mol(key[7:].strip()), key, source='用户 SMILES')
    if FORMULA.fullmatch(key):
        raise ReactionError(f'“{key}”只有分子式或未支持的缩写，不能唯一确定连接结构、离子配位或异构体。请用准确名称、SMILES:结构或 Reaction SMILES；多组分盐也应明确给出离子结构。')
    if not online:
        raise ReactionError(f'本地未收录“{key}”；请开启 PubChem 查询或提供 SMILES。')
    # 多个 CID 不选择首项；名称不唯一时要求用户提供准确结构。
    from urllib.parse import quote
    base = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug'
    try:
        response = requests.get(f'{base}/compound/name/{quote(key, safe="")}/cids/JSON', timeout=20)
        if response.status_code == 404:
            raise ReactionError(f'PubChem 未找到“{key}”，请提供准确英文名称或 SMILES。')
        response.raise_for_status()
        cids = response.json()['IdentifierList']['CID']
        if len(cids) != 1:
            raise ReactionError(f'“{key}”对应 {len(cids)} 个 PubChem 记录，存在歧义；请提供准确 SMILES。')
        response = requests.get(f'{base}/compound/cid/{cids[0]}/property/IsomericSMILES/JSON', timeout=20)
        response.raise_for_status()
        props = response.json()['PropertyTable']['Properties'][0]
        smiles = props.get('SMILES') or props.get('IsomericSMILES')
        if not smiles:
            raise ReactionError('PubChem 未返回可用的分子连接结构。')
        return Species(parse_mol(smiles), key, source=f'PubChem CID {cids[0]}')
    except (requests.RequestException, KeyError, TypeError, json.JSONDecodeError) as exc:
        raise ReactionError(f'PubChem 查询失败：{type(exc).__name__}。请检查网络或改用明确 SMILES。') from exc


def equation_side(side: str, online: bool) -> list[Species]:
    # 加号分隔符必须两侧有空格，避免拆开 H+ / Fe3+ 和 SMILES [Na+]。
    parts = re.split(r'\s+\+\s+', side.strip())
    result = []
    for part in parts:
        match = re.fullmatch(r'(?:(\d+)\s+)?(.+)', part.strip())
        if not match:
            raise ReactionError('存在空物种。请使用“反应物 + 反应物 -> 生成物”。')
        count = int(match[1] or 1)
        if not 1 <= count <= 1000:
            raise ReactionError('化学计量系数必须是 1–1000 的整数，且与物种之间留空格。')
        name = match[2]
        if re.match(r'^\d', name):
            raise ReactionError('系数后必须留空格，例如“2 H2O”。')
        if '+' in name and not (name in STRUCTURES or name.startswith('SMILES:')):
            raise ReactionError('分隔物种的加号两侧必须留空格；离子电荷加号不加空格。')
        item = resolve(name, online)
        item.coefficient = count
        result.append(item)
    return result


def inventory(items: list[Species]) -> tuple[Counter, int]:
    atoms = Counter()
    charge = 0
    for item in items:
        mol = Chem.AddHs(item.mol)
        for atom in mol.GetAtoms():
            # 同位素也分别守恒，避免把重氢和氢混同。
            name = (str(atom.GetIsotope()) if atom.GetIsotope() else '') + atom.GetSymbol()
            atoms[name] += item.coefficient
            charge += atom.GetFormalCharge() * item.coefficient
    return atoms, charge


def parse_reaction(text: str, online: bool = True) -> Reaction:
    text = text.strip()
    if not text or len(text) > 10000:
        raise ReactionError('请输入反应式（最多 10000 字符）。')
    agents = []
    if '>' in text and '->' not in text and '=>' not in text and '→' not in text:
        if text.count('>') != 2:
            raise ReactionError('Reaction SMILES 格式应为 reactants>agents>products 或 reactants>>products。')
        try:
            rxn = rdChemReactions.ReactionFromSmarts(text, useSmiles=True)
            if rxn is None:
                raise ValueError('空反应')
            def convert(mols):
                return [Species(parse_mol(Chem.MolToSmiles(m)), '', source='Reaction SMILES') for m in mols]
            left, right = convert(rxn.GetReactants()), convert(rxn.GetProducts())
            agents = convert(rxn.GetAgents())
        except Exception as exc:
            raise ReactionError(f'Reaction SMILES 无法解析：{exc}') from exc
    else:
        parts = re.split(r'\s*(?:->|→|=>)\s*', text)
        if len(parts) != 2:
            raise ReactionError('需要唯一的反应箭头 -> 或 →；Reaction SMILES 使用 >>。')
        left, right = [equation_side(side, online) for side in parts]
    if not left or not right:
        raise ReactionError('反应物和生成物均不能为空。')
    if len(left) + len(right) > 30:
        raise ReactionError('最多支持 30 个反应物及生成物条目；请拆分复杂反应。')
    a, qa = inventory(left)
    b, qb = inventory(right)
    differences = [f'{e}: 左 {a[e]} / 右 {b[e]}' for e in sorted(a.keys() | b.keys()) if a[e] != b[e]]
    if qa != qb:
        differences.append(f'净电荷: 左 {qa} / 右 {qb}')
    if differences:
        raise ReactionError('反应不守恒，已停止生成：' + '；'.join(differences))
    balance = '元素、同位素及净电荷守恒检查通过（不代表反应机理或反应可行性已验证）。'
    return Reaction(left, right, agents, balance)


def parse_input(text: str, online: bool = True, mode: str = 'auto') -> Reaction | Species:
    """无箭头的自动输入按一个完整名称处理，不拆开名称中的连字符或盐组分。"""
    text = text.strip()
    if not text or len(text) > 10000:
        raise ReactionError('请输入分子名称、SMILES 或反应式（最多 10000 字符）。')
    if mode not in ('auto', 'molecule', 'reaction'):
        raise ReactionError('输入模式必须为 auto、molecule 或 reaction。')
    if mode == 'reaction' or (mode == 'auto' and re.search(r'[><→←⇌]', text)):
        return parse_reaction(text, online)
    if mode == 'auto' and not text.startswith('SMILES:') and re.search(r'\s+\+\s+', text):
        raise ReactionError('检测到多个物种，但没有反应箭头。单个分子请输入一个名称或完整 SMILES；反应式请补充 ->。')
    return resolve(text, online)


def font(size: int):
    paths = ['C:/Windows/Fonts/msyh.ttc', '/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc', '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf']
    for path in paths:
        if Path(path).exists():
            return ImageFont.truetype(path, size)
    return ImageFont.load_default(size=size)


def straight_fragment_coords(mol: Chem.Mol) -> bool:
    """树状、无指定立体构型的组分：主链共线，直属支链垂直展开。"""
    if mol.GetRingInfo().NumRings():
        return False
    if any(a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in mol.GetAtoms()):
        return False
    if any(b.GetStereo() != Chem.BondStereo.STEREONONE for b in mol.GetBonds()):
        return False
    from collections import deque
    def farthest(start):
        queue = deque([start])
        previous = {start: None}
        distance = {start: 0}
        while queue:
            node = queue.popleft()
            for neighbor in mol.GetAtomWithIdx(node).GetNeighbors():
                idx = neighbor.GetIdx()
                if idx not in previous:
                    previous[idx] = node
                    distance[idx] = distance[node] + 1
                    queue.append(idx)
        def endpoint_rank(idx):
            atom = mol.GetAtomWithIdx(idx)
            terminal_multiple = atom.GetDegree() == 1 and atom.GetBonds()[0].GetBondTypeAsDouble() > 1
            return distance[idx], not terminal_multiple, -idx
        end = max(distance, key=endpoint_rank)
        return end, previous
    start, _ = farthest(0)
    end, previous = farthest(start)
    path = [end]
    while previous[path[-1]] is not None:
        path.append(previous[path[-1]])
    if path[0] > path[-1]:
        path.reverse()
    coords = {idx: Point2D(i * 2.0, 0.0) for i, idx in enumerate(path)}
    trunk = set(path)
    for idx in path:
        branches = sorted(n.GetIdx() for n in mol.GetAtomWithIdx(idx).GetNeighbors() if n.GetIdx() not in trunk)
        if len(branches) > 2:
            return False  # 高配位物种保留 RDKit 的常规几何表达。
        for j, branch in enumerate(branches):
            coords[branch] = Point2D(coords[idx].x, 2.0 if j == 0 else -2.0)
    rdDepictor.Compute2DCoords(mol, canonOrient=False, coordMap=coords)
    return True


def orient_fragment_horizontal(mol: Chem.Mol):
    """复杂环系/立体结构只刚性旋转，不强制拉直，不镜像。"""
    rdDepictor.Compute2DCoords(mol, useRingTemplates=True)
    conf = mol.GetConformer()
    positions = [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]
    cx = sum(p.x for p in positions) / len(positions)
    cy = sum(p.y for p in positions) / len(positions)
    xx = sum((p.x - cx) ** 2 for p in positions)
    yy = sum((p.y - cy) ** 2 for p in positions)
    xy = sum((p.x - cx) * (p.y - cy) for p in positions)
    angle = 0.5 * math.atan2(2 * xy, xx - yy)
    cos_a, sin_a = math.cos(angle), math.sin(angle)
    for i, p in enumerate(positions):
        x, y = p.x - cx, p.y - cy
        conf.SetAtomPosition(i, (x * cos_a + y * sin_a, -x * sin_a + y * cos_a, 0))


def prepare_drawing_mol(item: Species) -> Chem.Mol:
    """直线主链 + 正交支链；多组分物种（含反离子）在同一基线上排布。"""
    mol = Chem.Mol(item.mol)
    mappings = []
    fragments = Chem.GetMolFrags(mol, asMols=True, fragsMolAtomMapping=mappings)
    conf = Chem.Conformer(mol.GetNumAtoms())
    conf.Set3D(False)
    offset = 0.0
    order = sorted(range(len(fragments)), key=lambda i: (-fragments[i].GetNumAtoms(), i))
    for i in order:
        fragment = fragments[i]
        straight = straight_fragment_coords(fragment)
        if not straight:
            orient_fragment_horizontal(fragment)
        positions = fragment.GetConformer().GetPositions()
        low = min(positions[:, 0])
        high = max(positions[:, 0])
        center_y = 0 if straight else (min(positions[:, 1]) + max(positions[:, 1])) / 2
        for local, original in enumerate(mappings[i]):
            x, y, z = positions[local]
            conf.SetAtomPosition(original, (x - low + offset, y - center_y, 0))
        offset += high - low + 3.0
    mol.RemoveAllConformers()
    mol.AddConformer(conf)
    return mol


def drawing_size(item: Species) -> tuple[int, int]:
    mol = prepare_drawing_mol(item)
    coords = mol.GetConformer().GetPositions()
    dx = max(coords[:, 0]) - min(coords[:, 0])
    dy = max(coords[:, 1]) - min(coords[:, 1])
    # 标签需要额外空间；画布按实际结构跨度设置，减少大片空白。
    return max(300, min(1800, int(dx * 40 + 160))), max(280, min(1200, int(dy * 40 + 120)))


def molecule_drawing(item: Species, width: int, height: int, svg: bool):
    mol = prepare_drawing_mol(item)
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height) if svg else rdMolDraw2D.MolDraw2DCairo(width, height)
    opts = drawer.drawOptions()
    opts.useBWAtomPalette()
    opts.clearBackground = True
    opts.padding = 0.08
    opts.bondLineWidth = 2.2
    opts.multipleBondOffset = 0.075  # 多重键平行线间距为默认值的一半。
    opts.fixedBondLength = 65
    opts.maxFontSize = 30
    # 半展开表示：保留真实键线，将隐式氢写为 CH3/CH2/CH/NH/NH2 等基团。
    # 不向分子添加氢原子，也不改变键、手性、形式电荷或同位素。
    # 特殊碳原子由 RDKit 原生标签显示，以保留其同位素/电荷/自由基。
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (6, 7, 8, 16):
            continue
        if atom.GetDegree() == 0:
            continue  # H2O、NH3 等独立物种使用 RDKit 的常规标签顺序。
        if atom.GetIsotope() or atom.GetFormalCharge() or atom.GetNumRadicalElectrons():
            continue
        hydrogens = atom.GetTotalNumHs()
        label = atom.GetSymbol()
        if hydrogens:
            label += 'H' + (f'<sub>{hydrogens}</sub>' if hydrogens > 1 else '')
        opts.atomLabels[atom.GetIdx()] = label
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def render(reaction: Reaction, conditions: str = '') -> tuple[bytes, bytes]:
    if len(conditions) > 240 or len(conditions.splitlines()) > 4:
        raise ReactionError('条件最多 240 字符、4 行；请简化标注。')
    labels = [x.label or Chem.MolToSmiles(x.mol) for x in reaction.agents]
    if labels:
        conditions = conditions + ('\n' if conditions else '') + '试剂: ' + ', '.join(labels)
    # 每个结构独占一个区域；复杂分子扩大区域，不压缩成不可读的小图。
    items = reaction.reactants + reaction.products
    sizes = [drawing_size(item) for item in items]
    widths = [w for w, h in sizes]
    arrow_width = max(360, max((len(line) for line in conditions.splitlines()), default=0) * 26 + 60)
    height = max(360, max(h for w, h in sizes) + 80)
    height = min(height, 1200)
    width = sum(widths) + 70 * (len(items)-2) + arrow_width + 80 + 60 * sum(i.coefficient != 1 for i in items)
    if width * height > 30000000:
        raise ReactionError('图像过大；请减少反应物或条件长度。')
    canvas = Image.new('RGB', (width, height), 'white')
    draw = ImageDraw.Draw(canvas)
    svgparts = [f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="{width}" height="{height}" viewBox="0 0 {width} {height}">', '<rect width="100%" height="100%" fill="white"/>']
    def write(text, x, y, size=34):
        draw.text((x, y), text, fill='black', font=font(size), anchor='mm')
        svgparts.append(f'<text x="{x}" y="{y}" dominant-baseline="middle" text-anchor="middle" font-family="Microsoft YaHei,Arial,sans-serif" font-size="{size}">{escape(text)}</text>')
    x, index = 40, 0
    for side_index, side in enumerate([reaction.reactants, reaction.products]):
        for j, item in enumerate(side):
            if j:
                write('+', x+35, height/2, 42)
                x += 70
            if item.coefficient != 1:
                write(str(item.coefficient), x+25, height/2, 38)
                x += 60
            w = widths[index]
            index += 1
            data = molecule_drawing(item, w, height-80, False)
            canvas.paste(Image.open(io.BytesIO(data)).convert('RGB'), (int(x), 40))
            raw = molecule_drawing(item, w, height-80, True)
            raw = raw[raw.index('<svg'):]
            # 根 SVG 移到独立面板，保留 RDKit 原始矢量内容。
            raw = raw.replace('<svg ', f'<svg x="{x}" y="40" ', 1)
            svgparts.append(raw)
            x += w
        if side_index == 0:
            y = height/2
            draw.line((x+25, y, x+arrow_width-25, y), fill='black', width=4)
            draw.polygon([(x+arrow_width-25,y),(x+arrow_width-47,y-11),(x+arrow_width-47,y+11)], fill='black')
            svgparts.append(f'<path d="M {x+25} {y} H {x+arrow_width-25} M {x+arrow_width-47} {y-11} L {x+arrow_width-25} {y} L {x+arrow_width-47} {y+11}" fill="none" stroke="black" stroke-width="4"/>')
            for row, line in enumerate(conditions.splitlines()):
                write(line, x+arrow_width/2, y-45-(len(conditions.splitlines())-1-row)*34, 26)
            x += arrow_width
    svgparts.append('</svg>')
    out = io.BytesIO()
    canvas.save(out, format='PNG', dpi=(300,300))
    return out.getvalue(), '\n'.join(svgparts).encode('utf-8')


def render_molecule(item: Species, notes: str = '') -> tuple[bytes, bytes]:
    """单个分子独立画布，不添加反应箭头、不执行反应守恒检查。"""
    if len(notes) > 240 or len(notes.splitlines()) > 4:
        raise ReactionError('备注最多 240 字符、4 行；请简化标注。')
    width, mol_height = drawing_size(item)
    width = max(500, width)
    width = max(width, max((len(line) for line in notes.splitlines()), default=0) * 26 + 80)
    height = mol_height + 34 * len(notes.splitlines()) + (30 if notes else 0)
    png_mol = molecule_drawing(item, width, mol_height, False)
    canvas = Image.new('RGB', (width, height), 'white')
    canvas.paste(Image.open(io.BytesIO(png_mol)).convert('RGB'), (0, 0))
    raw = molecule_drawing(item, width, mol_height, True)
    raw = raw[raw.index('<svg'):]
    svgparts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">', '<rect width="100%" height="100%" fill="white"/>', raw]
    draw = ImageDraw.Draw(canvas)
    for row, line in enumerate(notes.splitlines()):
        y = mol_height + 20 + row * 34
        draw.text((width / 2, y), line, fill='black', font=font(26), anchor='mm')
        svgparts.append(f'<text x="{width / 2}" y="{y}" dominant-baseline="middle" text-anchor="middle" font-family="Microsoft YaHei,Arial,sans-serif" font-size="26">{escape(line)}</text>')
    svgparts.append('</svg>')
    out = io.BytesIO()
    canvas.save(out, format='PNG', dpi=(300, 300))
    return out.getvalue(), '\n'.join(svgparts).encode('utf-8')


def generate_image(text: str, conditions: str = '', online: bool = True, mode: str = 'auto') -> tuple[bytes, bytes, str]:
    parsed = parse_input(text, online, mode)
    if isinstance(parsed, Species):
        png, svg = render_molecule(parsed, conditions)
        formula = rdMolDescriptors.CalcMolFormula(parsed.mol)
        return png, svg, f'单个分子：{formula}。数据来源：{parsed.source}。单分子不涉及反应守恒检查。'
    png, svg = render(parsed, conditions)
    sources = sorted(set(i.source for i in parsed.reactants + parsed.products + parsed.agents))
    return png, svg, parsed.balance + ' 数据来源：' + ', '.join(sources)


def launch_gui():
    import tkinter as tk
    from tkinter import ttk, filedialog
    from PIL import ImageTk
    root = tk.Tk()
    root.title('化学二维分子与反应图生成器')
    root.geometry('1100x760')
    root.minsize(760, 560)
    frame = ttk.Frame(root, padding=14)
    frame.pack(fill='both', expand=True)
    mode_frame = ttk.Frame(frame)
    mode_frame.pack(fill='x')
    mode = tk.StringVar(value='auto')
    for label, value in [('自动识别', 'auto'), ('单个分子', 'molecule'), ('反应式', 'reaction')]:
        ttk.Radiobutton(mode_frame, text=label, variable=mode, value=value).pack(side='left', padx=6)
    ttk.Label(frame, text='单个分子名称 / SMILES:结构 / 化学方程式 / Reaction SMILES').pack(anchor='w')
    ttk.Label(frame, text='直线半展开：主链水平，支链上下展开，反离子靠右；环系/立体结构保留必要几何').pack(anchor='w')
    entry = ttk.Entry(frame, font=('Microsoft YaHei', 12))
    entry.pack(fill='x', pady=6)
    entry.insert(0, 'CH3COOH + C2H5OH -> CH3COOC2H5 + H2O')
    ttk.Label(frame, text='反应条件 / 单分子备注（可分行，最多4行）').pack(anchor='w')
    conditions = tk.Text(frame, height=3, font=('Microsoft YaHei', 11))
    conditions.pack(fill='x', pady=6)
    controls = ttk.Frame(frame)
    controls.pack(fill='x')
    fmt = tk.StringVar(value='PNG')
    for value in ('PNG', 'SVG'):
        ttk.Radiobutton(controls, text=value, variable=fmt, value=value).pack(side='left', padx=6)
    online = tk.BooleanVar(value=True)
    ttk.Checkbutton(controls, text='允许 PubChem 名称查询（需网络）', variable=online).pack(side='left', padx=12)
    preview = tk.Canvas(frame, background='white', highlightthickness=1, highlightbackground='#ccc')
    preview.pack(fill='both', expand=True, pady=12)
    status = tk.StringVar(value='输入准确结构后点击生成。分子式存在歧义时将拒绝生成。')
    ttk.Label(frame, textvariable=status, wraplength=1000, foreground='#a02020').pack(fill='x')
    state = {'png': None, 'svg': None, 'photo': None}
    def update_preview(event=None):
        if state['png']:
            img = Image.open(io.BytesIO(state['png']))
            img.thumbnail((max(1, preview.winfo_width()-10), max(1, preview.winfo_height()-10)))
            state['photo'] = ImageTk.PhotoImage(img)
            preview.delete('all')
            preview.create_image(preview.winfo_width()/2, preview.winfo_height()/2, image=state['photo'])
    preview.bind('<Configure>', update_preview)
    def save():
        selected = fmt.get().lower()
        path = filedialog.asksaveasfilename(defaultextension='.'+selected, filetypes=[(selected.upper(), '*.'+selected)])
        if path:
            try:
                Path(path).write_bytes(state[selected])
                status.set('已保存：'+path)
            except OSError as exc:
                status.set('保存失败：'+str(exc))
    save_button = ttk.Button(controls, text='保存图片', command=save, state='disabled')
    save_button.pack(side='right')
    # 后台线程只做化学与网络工作，Tk 主线程定时取结果。
    from queue import Queue, Empty
    results = Queue()
    def poll():
        try:
            result = results.get_nowait()
        except Empty:
            root.after(100, poll)
            return
        generate_button.config(state='normal')
        if isinstance(result, Exception):
            status.set(str(result))
        else:
            png, svg, message = result
            state.update(png=png, svg=svg)
            save_button.config(state='normal')
            status.set(message)
            update_preview()
    def generate():
        equation, cond, allow_online, input_mode = entry.get(), conditions.get('1.0', 'end').strip(), online.get(), mode.get()
        state.update(png=None, svg=None)
        preview.delete('all')
        save_button.config(state='disabled')
        generate_button.config(state='disabled')
        status.set('正在解析结构并绘制（反应式会验证守恒）…')
        def worker():
            try:
                results.put(generate_image(equation, cond, allow_online, input_mode))
            except Exception as exc:
                results.put(exc)
        threading.Thread(target=worker, daemon=True).start()
        root.after(100, poll)
    generate_button = ttk.Button(controls, text='生成图片', command=generate)
    generate_button.pack(side='right', padx=12)
    root.mainloop()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--equation', '--input', dest='equation', help='分子名称、SMILES 或反应式；省略则打开中文图形界面')
    parser.add_argument('--mode', choices=['auto', 'molecule', 'reaction'], default='auto')
    parser.add_argument('--conditions', default='')
    parser.add_argument('--output', default='output/reaction', help='输出文件前缀，同时生成 PNG/SVG')
    parser.add_argument('--offline', action='store_true')
    args = parser.parse_args()
    if args.equation is None:
        launch_gui()
        return
    try:
        png, svg, message = generate_image(args.equation, args.conditions, not args.offline, args.mode)
        target = Path(args.output)
        target.parent.mkdir(parents=True, exist_ok=True)
        Path(str(target)+'.png').write_bytes(png)
        Path(str(target)+'.svg').write_bytes(svg)
        print(message)
        print('已生成：', str(target)+'.png', str(target)+'.svg')
    except (ReactionError, OSError) as exc:
        parser.exit(2, f'错误：{exc}\n')

if __name__ == '__main__':
    main()
