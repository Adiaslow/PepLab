# peplab/analysis/strategies/factories/reactivity_analyzer_factory.py

"""
Factory module for creating reactivity-specific site analyzers.

This module defines a factory for creating instances of reactivity-specific site analyzers based on the specified type.
"""

from typing import Type
from ....core.graph.composites import MolecularGraph
from ..reactive_site_strategy import ReactiveSiteStrategy
from ..reactivity.amine_pattern_analyzer import AminePatternAnalyzer
from ..reactivity.carboxyl_pattern_analyzer import CarboxylPatternAnalyzer
from ..reactivity.azide_pattern_analyzer import AzidePatternAnalyzer
from ..reactivity.alkyne_pattern_analyzer import AlkynePatternAnalyzer

class ReactivityAnalyzerFactory:
    """
    Factory for creating reactivity-specific site analyzer instances.
    """

    @staticmethod
    def create_analyzer(analyzer_type: str, graph: MolecularGraph) -> ReactiveSiteStrategy:
        """
        Creates an instance of a reactivity-specific site analyzer based on the specified type.

        Args:
            analyzer_type (str): The type of analyzer to create.
            graph (MolecularGraph): The molecular graph to be analyzed.

        Returns:
            ReactiveSiteStrategy: An instance of the specified reactivity-specific site analyzer.
        """
        analyzers = {
            "amine": AminePatternAnalyzer,
            "carboxyl": CarboxylPatternAnalyzer,
            "azide": AzidePatternAnalyzer,
            "alkyne": AlkynePatternAnalyzer,
            # Add other analyzers here as needed
        }

        if analyzer_type in analyzers:
            return analyzers[analyzer_type](graph)
        else:
            raise ValueError(f"Unknown analyzer type: {analyzer_type}")
