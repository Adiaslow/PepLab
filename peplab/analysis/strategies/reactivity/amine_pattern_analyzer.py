# peplab/analysis/strategies/amine_pattern_analyzer.py

"""
Concrete strategy for analyzing primary and secondary amine patterns in a molecular graph.
"""

from typing import List
from uuid import uuid4
from ....core.reaction import ReactivePattern
from ....core.reaction import ReactiveType
from ....core.graph.composites import MolecularGraph
from ..reactive_site_strategy import ReactiveSiteStrategy

class AminePatternAnalyzer:
    """
    Concrete strategy for analyzing primary and secondary amine patterns in a molecular graph.
    """

    def __init__(self, graph: MolecularGraph):
        """
        Initializes the AminePatternAnalyzer with a molecular graph.

        Args:
            graph (MolecularGraph): The molecular graph to be analyzed.
        """
        self.graph = graph

    def analyze(self) -> List[ReactivePattern]:
        """Find primary and secondary amine patterns."""
        reactive_patterns = []
        for node in self.graph.nodes:
            if node.element != 'N':
                continue

            neighbors = self.graph.get_neighbors(node.index)
            h_neighbors = [(n, e) for n, e in neighbors if n.element == 'H']
            heavy_neighbors = [(n, e) for n, e in neighbors if n.element != 'H']

            # Check for NH2 pattern
            if len(h_neighbors) == 2 and len(heavy_neighbors) == 1:
                heavy_atom, heavy_bond = heavy_neighbors[0]
                if heavy_atom.element == 'C' and heavy_bond.properties.get('bond_type') == 'SINGLE':
                    node.properties['is_reactive_nuc'] = True
                    reactive_patterns.append(ReactivePattern(
                        type=ReactiveType.NH2,
                        atoms=[node, *[n for n, _ in h_neighbors], heavy_atom],
                        bonds=[e for _, e in neighbors],
                        pattern_id=str(uuid4())
                    ))

            # Check for NH pattern
            elif len(h_neighbors) == 1 and len(heavy_neighbors) == 2:
                if all(n.element == 'C' and e.properties.get('bond_type') == 'SINGLE'
                      for n, e in heavy_neighbors):
                    node.properties['is_reactive_nuc'] = True
                    reactive_patterns.append(ReactivePattern(
                        type=ReactiveType.NH,
                        atoms=[node, h_neighbors[0][0], *[n for n, _ in heavy_neighbors]],
                        bonds=[e for _, e in neighbors],
                        pattern_id=str(uuid4())
                    ))
        return reactive_patterns
