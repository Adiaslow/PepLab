# peplab/analysis/strategies/reactivity/azide_pattern_analyzer.py

"""
Concrete strategy for analyzing azide (N3) patterns in a molecular graph.
"""

from typing import List
from uuid import uuid4
from ....core.reaction import ReactivePattern
from ....core.reaction import ReactiveType
from ....core.graph.composites import MolecularGraph
from ..reactive_site_strategy import ReactiveSiteStrategy

class AzidePatternAnalyzer:
    """
    Concrete strategy for analyzing azide (N3) patterns in a molecular graph.
    """

    def __init__(self, graph: MolecularGraph):
        """
        Initializes the AzidePatternAnalyzer with a molecular graph.

        Args:
            graph (MolecularGraph): The molecular graph to be analyzed.
        """
        self.graph = graph

    def analyze(self) -> List[ReactivePattern]:
        """Find azide (N3) patterns."""
        reactive_patterns = []
        for node in self.graph.nodes:
            if node.element != 'N':
                continue

            neighbors = self.graph.get_neighbors(node.index)
            if len(neighbors) != 2:
                continue

            # Find middle nitrogen
            n2_pair = next(((n, e) for n, e in neighbors
                            if n.element == 'N' and e.properties.get('bond_type') == 'SINGLE'), None)
            if not n2_pair:
                continue

            n2, n1_n2_bond = n2_pair
            n2_neighbors = self.graph.get_neighbors(n2.index)
            if len(n2_neighbors) != 2:
                continue

            # Find terminal nitrogen
            n3_pair = next(((n, e) for n, e in n2_neighbors
                            if n.element == 'N' and n.index != node.index
                            and e.properties.get('bond_type') == 'TRIPLE'), None)
            if not n3_pair:
                continue

            n3, n2_n3_bond = n3_pair
            if len(self.graph.get_neighbors(n3.index)) == 1:
                node.properties['is_reactive_click'] = True
                reactive_patterns.append(ReactivePattern(
                    type=ReactiveType.AZIDE,
                    atoms=[node, n2, n3],
                    bonds=[n1_n2_bond, n2_n3_bond],
                    pattern_id=str(uuid4())
                ))
        return reactive_patterns
