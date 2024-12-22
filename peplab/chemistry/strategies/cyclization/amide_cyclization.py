# peplab/chemistry/strategies/cyclization/amide_cyclization.py

"""
Module for amide cyclization strategy.

This module provides the AmideCyclization class which implements the CyclizationStrategy protocol
for cyclizing molecules by forming amide bonds within themselves.
"""

from ..cyclization_strategy import CyclizationStrategy
from ..bonding import AmideBonding
from ....core.graph import MolecularGraph

class AmideCyclization(CyclizationStrategy):
    """Concrete strategy for amide cyclization."""

    def __init__(self):
        self.bonding_strategy = AmideBonding()

    def cyclize(self, molecule: MolecularGraph) -> MolecularGraph:
        """
        Cyclize a molecular graph by forming an amide bond within itself.

        Args:
            molecule (MolecularGraph): The molecular graph to cyclize.

        Returns:
            MolecularGraph: The cyclized molecular graph.
        """
        return self.bonding_strategy.form_bond(molecule, molecule)
