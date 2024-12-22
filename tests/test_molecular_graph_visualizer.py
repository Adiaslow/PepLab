import pytest
from peplab.visualization.graph.molecular_graph_visualizer import MolecularGraphVisualizer
import os

def test_show_molecule():
    smiles = "N[C@H](C(O)=O)CC12C[C@H]3C[C@@H](C2)C[C@@H](C1)C3"  # Adamanthane
    visualizer = MolecularGraphVisualizer(smiles)
    file_path = "test_molecule.png"
    visualizer.show_molecule(file_path)
    assert os.path.exists(file_path)
    os.remove(file_path)

def test_highlight_reactive_groups():
    smiles = "N[C@H](C(O)=O)CC12C[C@H]3C[C@@H](C2)C[C@@H](C1)C3"  # Adamanthane
    visualizer = MolecularGraphVisualizer(smiles)
    file_path = "test_highlighted_molecule.png"
    visualizer.highlight_reactive_groups(file_path)
    assert os.path.exists(file_path)
    os.remove(file_path)
