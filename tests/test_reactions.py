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
