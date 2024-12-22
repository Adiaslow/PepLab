# peplab/chemistry/strategies/cyclization/triazole_cyclization.py

"""
Module for triazole cyclization strategy.

This module provides the TriazoleCyclization class which implements the CyclizationStrategy protocol
for cyclizing molecules by forming triazole bonds within themselves.
"""

from ..cyclization_strategy import CyclizationStrategy
from ..bonding import TriazoleBonding
from ....core.graph import MolecularGraph

class TriazoleCyclization(CyclizationStrategy):
    """Concrete strategy for triazole cyclization."""

    def __init__(self):
        self.bonding_strategy = TriazoleBonding()

    def cyclize(self, molecule: MolecularGraph) -> MolecularGraph:
        """
        Cyclize a molecular graph by forming a triazole bond within itself.

        Args:
            molecule (MolecularGraph): The molecular graph to cyclize.

        Returns:
            MolecularGraph: The cyclized molecular graph.
        """
        return self.bonding_strategy.form_bond(molecule, molecule)
