# tests/test_molecular_graph_builder.py

import pytest
from peplab.chemistry.builders.molecular_graph_builder import MolecularGraphBuilder
from peplab.core.graph.composites import MolecularGraph

@pytest.fixture
def builder():
    return MolecularGraphBuilder()

def test_from_smiles_valid(builder):
    smiles = "CCO"  # Ethanol
    graph = builder.from_smiles(smiles)
    assert isinstance(graph, MolecularGraph)
    assert len(graph.nodes) == 9  # 3 atoms and 6 implicit H atoms
    assert len(graph.edges) == 8  # 6 C-H bonds, 1 C-C bond, 1 C-O bond
    # Check atom types
    atom_elements = [node.element for node in graph.nodes]
    assert sorted(atom_elements) == ['C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'O']
    # Check bond types
    bond_types = [edge.properties['bond_type'] for edge in graph.edges]
    assert bond_types.count('SINGLE') == 8

def test_from_smiles_invalid(builder):
    smiles = "invalid_smiles"
    with pytest.raises(ValueError):
        builder.from_smiles(smiles)

def test_from_smiles_adamanthane(builder):
    # Adamanthane - non-natural amino acid
    smiles = "N[C@H](C(O)=O)CC12C[C@H]3C[C@@H](C2)C[C@@H](C1)C3"
    graph = builder.from_smiles(smiles)
    assert isinstance(graph, MolecularGraph)
    assert len(graph.nodes) == 37  # 1N + 13C + 2O + 21H
    assert len(graph.edges) == 39  # 38 single bonds + 1 double bond
    # Check atom types
    atom_elements = [node.element for node in graph.nodes]
    assert atom_elements.count('N') == 1   # N-terminus (NH3+)
    assert atom_elements.count('C') == 13  # Alpha carbon + carbonyl + adamantyl carbons
    assert atom_elements.count('O') == 2   # Carboxylic acid oxygens
    assert atom_elements.count('H') == 21  # NH3+ (3) + alpha (1) + adamantyl hydrogens (17)
    # Check bond types
    bond_types = [edge.properties['bond_type'] for edge in graph.edges]
    assert bond_types.count('SINGLE') == 38
    assert bond_types.count('DOUBLE') == 1  # Carboxylic acid C=O

def test_from_smiles_aromatic(builder):
    smiles = "c1ccccc1"  # Benzene
    graph = builder.from_smiles(smiles)
    assert isinstance(graph, MolecularGraph)
    assert len(graph.nodes) == 12  # 6 carbon atoms and 6 implicit H atoms
    assert len(graph.edges) == 12  # 6 C-H bonds and 6 aromatic C-C bonds
    # Check atom types
    atom_elements = [node.element for node in graph.nodes]
    assert sorted(atom_elements) == ['C'] * 6 + ['H'] * 6
    # Check bond types
    bond_types = [edge.properties['bond_type'] for edge in graph.edges]
    assert bond_types.count('AROMATIC') == 6
    assert bond_types.count('SINGLE') == 6

def test_from_smiles_triple_bond(builder):
    smiles = "C#CC"  # Propyne
    graph = builder.from_smiles(smiles)
    assert isinstance(graph, MolecularGraph)
    assert len(graph.nodes) == 7  # 3 carbon atoms and 4 implicit H atoms
    assert len(graph.edges) == 6  # 4 C-H bonds, 1 C#C bond, 1 C-C bond
    # Check atom types
    atom_elements = [node.element for node in graph.nodes]
    assert sorted(atom_elements) == ['C', 'C', 'C', 'H', 'H', 'H', 'H']
    # Check bond types
    bond_types = [edge.properties['bond_type'] for edge in graph.edges]
    assert bond_types.count('SINGLE') == 5
    assert bond_types.count('TRIPLE') == 1

def test_from_smiles_empty(builder):
    smiles = ""
    with pytest.raises(ValueError):
        builder.from_smiles(smiles)

if __name__ == "__main__":
    pytest.main()
