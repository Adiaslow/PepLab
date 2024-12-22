# peplab/analysis/reactive_site_analyzer.py

"""
Context class for using reactive site analysis strategies.

This module defines the ReactiveSiteAnalyzer class, which uses different strategies for analyzing reactive sites in a molecular graph.
"""

from typing import List
from .strategies.reactive_site_strategy import ReactiveSiteStrategy
from ..core.reaction import ReactivePattern

class ReactiveSiteAnalyzer:
    """
    Context class for analyzing reactive sites using different strategies.

    Attributes:
        strategy (ReactiveSiteStrategy): The strategy used for analysis.
    """

    def __init__(self, strategy: ReactiveSiteStrategy):
        """
        Initializes the ReactiveSiteAnalyzer with a specific strategy.

        Args:
            strategy (ReactiveSiteStrategy): The strategy used for analysis.
        """
        self.strategy = strategy

    def analyze(self) -> List[ReactivePattern]:
        """
        Analyzes the molecular graph for reactive sites using the specified strategy.

        Returns:
            List[ReactivePattern]: List of identified reactive patterns.
        """
        return self.strategy.analyze()
