# peplab/chemistry/builders/peptide_builder.py

"""
Module for building peptides.

This module provides the PeptideBuilder class which allows for the construction of peptide objects.
"""

from .molecular_graph_builder import MolecularGraphBuilder
from ...core.graph import MolecularGraph

class PeptideBuilder(MolecularGraphBuilder):
    """
    Class for building peptides.
    """

    def __init__(self):
        super().__init__()

    def add_residue(self, residue: MolecularGraph) -> 'PeptideBuilder':
        """
        Adds a residue to the peptide.

        Args:
            residue (MolecularGraph): The residue to add.

        Returns:
            PeptideBuilder: The builder instance for chaining.
        """
        max_index = max((node.index for node in self._graph.nodes), default=-1) + 1
        for node in residue.nodes:
            node.index += max_index
        for edge in residue.edges:
            edge.from_idx += max_index
            edge.to_idx += max_index

        self._graph.nodes.extend(residue.nodes)
        self._graph.edges.extend(residue.edges)
        return self
