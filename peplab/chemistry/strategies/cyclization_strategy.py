# peplab/chemistry/strategies/cyclization_strategy.py

"""
Module for cyclization strategies used in peptide construction.

This module defines the CyclizationStrategy protocol that all cyclization strategies must implement.
"""

from typing import Protocol
from ...core.graph.molecule_graph import MolecularGraph

class CyclizationStrategy(Protocol):
    """Protocol for cyclization strategies."""

    def cyclize(self, molecule: MolecularGraph) -> MolecularGraph:
        """
        Cyclize a molecular graph.

        Args:
            molecule (MolecularGraph): The molecular graph to cyclize.

        Returns:
            MolecularGraph: The cyclized molecular graph.
        """
        ...
