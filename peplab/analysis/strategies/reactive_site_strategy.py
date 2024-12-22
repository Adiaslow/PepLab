# peplab/analysis/strategies/reactive_site_strategy.py

"""
Module defining the protocol for reactive site analyzers.

This module defines the ReactiveSiteStrategy protocol, which is implemented by various concrete strategies
for analyzing different types of reactive sites in a molecular graph.
"""

from typing import Protocol, List
from ...core.reaction import ReactivePattern
from ...core.graph.composites import MolecularGraph

class ReactiveSiteStrategy(Protocol):
    """
    Protocol for reactive site analysis strategies.

    Attributes:
        graph (MolecularGraph): The molecular graph to be analyzed.
    """

    graph: MolecularGraph

    def analyze(self) -> List[ReactivePattern]:
        """
        Analyzes the molecular graph for specific reactive sites.

        Returns:
            List[ReactivePattern]: List of identified reactive patterns.
        """
        ...
