# peplab/chemistry/strategies/bonding_strategy.py

"""
Module for bonding strategies used in peptide construction.

This module defines the BondingStrategy protocol that all bonding strategies must implement.
"""

from typing import Protocol
from ...core.graph.molecule_graph import MolecularGraph

class BondingStrategy(Protocol):
    """Protocol for bonding strategies."""

    def form_bond(self, res1: MolecularGraph, res2: MolecularGraph) -> MolecularGraph:
        """
        Form a bond between two molecular graphs.

        Args:
            res1 (MolecularGraph): The first molecular graph.
            res2 (MolecularGraph): The second molecular graph.

        Returns:
            MolecularGraph: The molecular graph with the newly formed bond.
        """
        ...
