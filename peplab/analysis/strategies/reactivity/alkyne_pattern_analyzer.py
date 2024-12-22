# peplab/analysis/strategies/reactivity/alkyne_pattern_analyzer.py

"""
Concrete strategy for analyzing terminal alkyne (C≡C) patterns in a molecular graph.
"""

from typing import List
from uuid import uuid4
from ....core.reaction import ReactivePattern
from ....core.reaction import ReactiveType
from ....core.graph.composites import MolecularGraph
from ..reactive_site_strategy import ReactiveSiteStrategy

class AlkynePatternAnalyzer:
    """
    Concrete strategy for analyzing terminal alkyne (C≡C) patterns in a molecular graph.
    """

    def __init__(self, graph: MolecularGraph):
        """
        Initializes the AlkynePatternAnalyzer with a molecular graph.

        Args:
            graph (MolecularGraph): The molecular graph to be analyzed.
        """
        self.graph = graph

    def analyze(self) -> List[ReactivePattern]:
        """Find terminal alkyne (C≡C) patterns."""
        reactive_patterns = []
        for node in self.graph.nodes:
            if node.element != 'C':
                continue

            neighbors = self.graph.get_neighbors(node.index)
            triple_bond_pair = next(((n, e) for n, e in neighbors
                                     if n.element == 'C' and e.properties.get('bond_type') == 'TRIPLE'), None)
            if not triple_bond_pair:
                continue

            other_c, triple_bond = triple_bond_pair

            # Get non-triple bond connections
            current_other_bonds = [(n, e) for n, e in neighbors if e != triple_bond]
            other_c_neighbors = self.graph.get_neighbors(other_c.index)
            other_c_other_bonds = [(n, e) for n, e in other_c_neighbors if e != triple_bond]

            # Check for terminal alkyne
            is_current_terminal = len(current_other_bonds) <= 1
            is_other_terminal = len(other_c_other_bonds) <= 1

            if ((is_current_terminal and len(other_c_other_bonds) <= 2) or
                (is_other_terminal and len(current_other_bonds) <= 2)):
                node.properties['is_reactive_click'] = True
                reactive_patterns.append(ReactivePattern(
                    type=ReactiveType.ALKYNE,
                    atoms=[node, other_c],
                    bonds=[triple_bond],
                    pattern_id=str(uuid4())
                ))
        return reactive_patterns
