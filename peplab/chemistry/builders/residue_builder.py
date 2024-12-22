# peplab/chemistry/builders/residue_builder.py

"""
Module for building residues.

This module provides the ResidueBuilder class which allows for the construction of residue objects.
"""

from .molecular_graph_builder import MolecularGraphBuilder

class ResidueBuilder(MolecularGraphBuilder):
    """
    Class for building residues.
    """

    def __init__(self):
        super().__init__()
