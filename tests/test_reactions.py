import io
import xml.etree.ElementTree as ET
import pytest
from PIL import Image
from chem_drawer import parse_reaction, render, ReactionError, inventory

EXAMPLES = ['CH3COOH + C2H5OH -> CH3COOC2H5 + H2O', '2 H2 + O2 -> 2 H2O', 'CC(=O)O.CCO>>CCOC(C)=O.O']

@pytest.mark.parametrize('equation', EXAMPLES)
def test_examples(equation, tmp_path):
    rxn = parse_reaction(equation, False)
    assert inventory(rxn.reactants) == inventory(rxn.products)
    png, svg = render(rxn, 'H2SO4\n80 °C / 乙醇')
    image = Image.open(io.BytesIO(png))
    image.load()
    assert image.info['dpi'][0] >= 299.99
    assert image.getpixel((0, 0)) == (255, 255, 255)
    assert image.getextrema()[0][0] < 100
    root = ET.fromstring(svg)
    assert len(root.findall('.//{http://www.w3.org/2000/svg}path')) > 5

@pytest.mark.parametrize('equation', ['C2H6O -> C2H6O', 'H2 + O2 -> H2O', 'H+ -> H2', '2H2 + O2 -> 2 H2O', 'C>>', '[*]>>[*]', 'H2+O2->H2O'])
def test_reject(equation):
    with pytest.raises(ReactionError):
        parse_reaction(equation, False)

@pytest.mark.parametrize('equation', ['H+ + OH- -> H2O', 'Na+ + Cl- -> NaCl', 'SMILES:[Na+].[Cl-] -> NaCl', '([Na+].[Cl-])>>([Na+].[Cl-])', '乙酸 + 乙醇 -> 乙酸乙酯 + 水'])
def test_ions_fragments_names(equation):
    rxn = parse_reaction(equation, False)
    assert inventory(rxn.reactants) == inventory(rxn.products)


def test_coefficients():
    rxn = parse_reaction(EXAMPLES[1], False)
    assert [x.coefficient for x in rxn.reactants] == [2, 1]
    assert rxn.products[0].coefficient == 2


def test_pubchem_unique_and_stereo(monkeypatch):
    from chem_drawer import resolve
    import requests
    class Response:
        status_code = 200
        def __init__(self, data):
            self.data = data
        def raise_for_status(self):
            pass
        def json(self):
            return self.data
    responses = iter([Response({'IdentifierList': {'CID': [123]}}), Response({'PropertyTable': {'Properties': [{'SMILES': 'C[C@H](O)C(=O)O'}]}})])
    monkeypatch.setattr(requests, 'get', lambda *a, **k: next(responses))
    item = resolve('specific stereoisomer')
    from rdkit import Chem
    assert '@' in Chem.MolToSmiles(item.mol)
    assert item.source == 'PubChem CID 123'


def test_pubchem_ambiguous(monkeypatch):
    from chem_drawer import resolve
    import requests
    class Response:
        status_code = 200
        def raise_for_status(self):
            pass
        def json(self):
            return {'IdentifierList': {'CID': [1, 2]}}
    monkeypatch.setattr(requests, 'get', lambda *a, **k: Response())
    with pytest.raises(ReactionError, match='歧义'):
        resolve('ambiguous name')


def test_group_not_split():
    rxn = parse_reaction('([Na+].[Cl-])>>([Na+].[Cl-])', False)
    assert len(rxn.reactants) == 1
    assert rxn.reactants[0].mol.GetNumAtoms() == 2


def test_isotope_imbalance():
    with pytest.raises(ReactionError, match='不守恒'):
        parse_reaction('[2H]O[2H]>>O', False)


def test_reaction_agents_drawn():
    rxn = parse_reaction('CC(=O)O.CCO>O=S(=O)(O)O>CCOC(C)=O.O', False)
    assert len(rxn.agents) == 1
    png, svg = render(rxn)
    assert '试剂: O=S(=O)(O)O' in svg.decode()


@pytest.mark.parametrize('text', ['乙醇', 'SMILES:CCO', 'H+', 'SMILES:[Na+].[Cl-]', 'Nα-(2,4-dinitro-5-fluorophenyl)-L-alaninamide'])
def test_single_molecule(text):
    from chem_drawer import generate_image, parse_input, Species
    item = parse_input(text, False)
    assert isinstance(item, Species)
    png, svg, message = generate_image(text, '单分子备注', False)
    im = Image.open(io.BytesIO(png)); im.load()
    assert im.info['dpi'][0] >= 299.99
    root = ET.fromstring(svg)
    assert not root.findall('./{http://www.w3.org/2000/svg}path')  # 外层无反应箭头
    assert '单个分子' in message and '反应守恒检查' in message


def test_named_chirality_and_formula():
    from chem_drawer import resolve
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
    item = resolve('Nα-(2,4-dinitro-5-fluorophenyl)-L-alaninamide', False)
    assert rdMolDescriptors.CalcMolFormula(item.mol) == 'C9H9FN4O5'
    assert [cip for index, cip in Chem.FindMolChiralCenters(item.mol)] == ['S']
    assert '5486955' in item.source
    # 芳环的 2、4 位硝基，5 位 F；SMILES 独立书写时仍为相同 S 结构。
    reference = Chem.MolFromSmiles('C[C@H](Nc1c([N+](=O)[O-])cc([N+](=O)[O-])c(F)c1)C(N)=O')
    assert Chem.MolToSmiles(item.mol) == Chem.MolToSmiles(reference)


@pytest.mark.parametrize('text', ['', 'C2H6O', 'H2 + O2', 'H2 ->', '乙醇 ⇌ 水'])
def test_auto_rejects_ambiguous_or_malformed(text):
    from chem_drawer import parse_input
    with pytest.raises(ReactionError):
        parse_input(text, False)


def test_modes_and_reaction_compatibility():
    from chem_drawer import parse_input, Reaction, Species
    assert isinstance(parse_input(EXAMPLES[0], False), Reaction)
    assert isinstance(parse_input(EXAMPLES[2], False), Reaction)
    assert isinstance(parse_input('SMILES:N->[Cu+2]', False, 'molecule'), Species)
    assert isinstance(parse_input('2 H2 + O2 => 2 H2O', False), Reaction)
    with pytest.raises(ReactionError):
        parse_input('乙醇', False, 'reaction')


@pytest.mark.parametrize('smiles', ['CCO', 'c1ccccc1', '[13CH3][C@H](O)C(=O)[O-]', '[CH2]CC'])
def test_mixed_formula_visible_carbons_and_preserved_structure(smiles):
    from chem_drawer import Species, parse_mol, molecule_drawing
    from rdkit import Chem
    item = Species(parse_mol(smiles), '')
    before = Chem.MolToSmiles(item.mol)
    svg = molecule_drawing(item, 900, 600, True)
    root = ET.fromstring(svg)
    # 实际输出的原子字形路径应覆盖全部碳，而非只存在键线中的端点编号。
    labelled = {p.attrib.get('class') for p in root.iter() if p.tag.endswith('path')}
    for atom in item.mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            assert f'atom-{atom.GetIdx()}' in labelled
    assert Chem.MolToSmiles(item.mol) == before


def test_horizontal_orientation_preserves_stereo():
    from chem_drawer import resolve, prepare_drawing_mol
    from rdkit import Chem
    item = resolve('Nα-(2,4-dinitro-5-fluorophenyl)-L-alaninamide', False)
    before = Chem.MolToSmiles(item.mol)
    drawn = prepare_drawing_mol(item)
    positions = drawn.GetConformer().GetPositions()
    assert positions[:, 0].var() >= positions[:, 1].var()
    assert Chem.MolToSmiles(drawn) == before
    assert [cip for i, cip in Chem.FindMolChiralCenters(drawn)] == ['S']
    assert item.mol.GetNumConformers() == 0


def test_straight_hydrazide_backbone_branches_and_counterion():
    from chem_drawer import Species, parse_mol, prepare_drawing_mol
    from rdkit import Chem
    item = Species(parse_mol('C[N+](C)(C)CC(=O)NN.[Cl-]'), '')
    before = Chem.MolToSmiles(item.mol)
    drawn = prepare_drawing_mol(item)
    coords = drawn.GetConformer().GetPositions()
    chain = [0, 1, 4, 5, 7, 8]
    assert max(coords[chain, 1]) - min(coords[chain, 1]) < 1e-8
    assert all(coords[b, 0] > coords[a, 0] for a, b in zip(chain, chain[1:]))
    assert coords[2, 1] > 0 and coords[3, 1] < 0
    assert coords[2, 0] == pytest.approx(coords[1, 0])
    assert coords[3, 0] == pytest.approx(coords[1, 0])
    assert coords[6, 1] > 0 and coords[6, 0] == pytest.approx(coords[5, 0])
    assert coords[9, 1] == pytest.approx(0)
    assert coords[9, 0] > max(coords[chain, 0])
    assert len(Chem.GetMolFrags(drawn)) == 2
    assert Chem.MolToSmiles(drawn) == before
    assert item.mol.GetNumConformers() == 0


@pytest.mark.parametrize('smiles', ['F/C=C/F', 'F/C=C\\F', 'C[C@H](O)C(=O)O', 'c1ccccc1'])
def test_stereo_and_rings_not_forced_collinear(smiles):
    from chem_drawer import Species, parse_mol, prepare_drawing_mol
    from rdkit import Chem
    item = Species(parse_mol(smiles), '')
    drawn = prepare_drawing_mol(item)
    coords = drawn.GetConformer().GetPositions()
    assert coords[:, 1].max() - coords[:, 1].min() > 0.1
    assert Chem.MolToSmiles(drawn) == Chem.MolToSmiles(item.mol)
    restored = Chem.MolFromMolBlock(Chem.MolToMolBlock(drawn))
    assert Chem.MolToSmiles(restored) == Chem.MolToSmiles(item.mol)


def test_custom_canvas_preserves_pixels():
    from chem_drawer import generate_image
    png, svg, _ = generate_image(EXAMPLES[0], online=False)
    original = Image.open(io.BytesIO(png)).convert('RGB')
    custom, vector, _ = generate_image(EXAMPLES[0], online=False, width=original.width+200, height=original.height+100)
    image = Image.open(io.BytesIO(custom)).convert('RGB')
    assert image.size == (original.width+200, original.height+100)
    assert image.crop((100, 50, original.width+100, original.height+50)).tobytes() == original.tobytes()
    root = ET.fromstring(vector)
    assert root.attrib['width'] == str(image.width)
    assert root.attrib['height'] == str(image.height)


def test_fit_preserves_aspect_and_output_size():
    from chem_drawer import generate_image
    from PIL import ImageChops
    png, _, _ = generate_image(EXAMPLES[0], online=False)
    original = Image.open(io.BytesIO(png)).convert('RGB')
    fitted, svg, _ = generate_image(EXAMPLES[0], online=False, width=600, height=400, size_mode='fit')
    image = Image.open(io.BytesIO(fitted)).convert('RGB')
    assert image.size == (600, 400)
    def aspect(im):
        ink = ImageChops.difference(im, Image.new('RGB', im.size, 'white')).convert('L').point(lambda v: 255 if v > 64 else 0)
        x1,y1,x2,y2 = ink.getbbox()
        return (x2-x1)/(y2-y1)
    assert aspect(image) == pytest.approx(aspect(original), rel=0.03)
    root = ET.fromstring(svg)
    group = root.find('{http://www.w3.org/2000/svg}g')
    # 矢量只用一个统一比例，横纵没有独立的 scale(x,y)。
    assert 'scale(' in group.attrib['transform']
    assert ',' not in group.attrib['transform']
    assert Image.open(io.BytesIO(fitted)).info['dpi'][0] >= 299.99


@pytest.mark.parametrize('kwargs', [{'width':64}, {'width':'abc'}, {'height':0}, {'width':10001}, {'width':10000,'height':10000}])
def test_invalid_or_clipping_size_rejected(kwargs):
    from chem_drawer import generate_image
    with pytest.raises(ReactionError):
        generate_image(EXAMPLES[0], online=False, **kwargs)


def test_compact_small_species_and_arrow():
    from chem_drawer import compact_panel, resolve, generate_image
    im, svg = compact_panel(resolve('H2O', False))
    assert im.width < 100
    png, svg, _ = generate_image(EXAMPLES[0], online=False)
    root = ET.fromstring(svg)
    arrow = root.find('./{http://www.w3.org/2000/svg}path')
    import re
    match = re.match(r'M (\S+) (\S+) H (\S+)', arrow.attrib['d'])
    assert float(match[3])-float(match[1]) <= 120
