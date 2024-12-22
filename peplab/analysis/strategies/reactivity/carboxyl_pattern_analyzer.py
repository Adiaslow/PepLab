# peplab/analysis/strategies/reactivity/carboxyl_pattern_analyzer.py

"""
Concrete strategy for analyzing carboxyl (COOH) patterns in a molecular graph.
"""

from typing import List
from uuid import uuid4
from ....core.reaction import ReactivePattern
from ....core.reaction import ReactiveType
from ....core.graph.composites import MolecularGraph
from ..reactive_site_strategy import ReactiveSiteStrategy

class CarboxylPatternAnalyzer:
    """
    Concrete strategy for analyzing carboxyl (COOH) patterns in a molecular graph.
    """

    def __init__(self, graph: MolecularGraph):
        """
        Initializes the CarboxylPatternAnalyzer with a molecular graph.

        Args:
            graph (MolecularGraph): The molecular graph to be analyzed.
        """
        self.graph = graph

    def analyze(self) -> List[ReactivePattern]:
        """Find carboxyl (COOH) patterns."""
        reactive_patterns = []
        for node in self.graph.nodes:
            if node.element != 'C':
                continue

            neighbors = self.graph.get_neighbors(node.index)
            o_neighbors = [(n, e) for n, e in neighbors if n.element == 'O']
            heavy_neighbors = [(n, e) for n, e in neighbors if n.element not in {'O', 'H'}]

            if len(o_neighbors) != 2 or not heavy_neighbors:
                continue

            # Look for C(=O)OH pattern
            double_o = None
            single_o = None
            for o_atom, o_bond in o_neighbors:
                if o_bond.properties.get('bond_type') == 'DOUBLE':
                    double_o = (o_atom, o_bond)
                elif o_bond.properties.get('bond_type') == 'SINGLE':
                    o_neighbors2 = self.graph.get_neighbors(o_atom.index)
                    h_neighbors = [(n, e) for n, e in o_neighbors2 if n.element == 'H']
                    if len(h_neighbors) == 1:
                        single_o = (o_atom, o_bond)

            if double_o and single_o and len(heavy_neighbors) == 1:
                heavy_atom, heavy_bond = heavy_neighbors[0]
                if heavy_bond.properties.get('bond_type') == 'SINGLE':
                    node.properties['is_reactive_elec'] = True
                    reactive_patterns.append(ReactivePattern(
                        type=ReactiveType.COOH,
                        atoms=[node, double_o[0], single_o[0],
                               next(n for n, _ in self.graph.get_neighbors(single_o[0].index)
                                    if n.element == 'H')],
                        bonds=[double_o[1], single_o[1],
                               next(e for _, e in self.graph.get_neighbors(single_o[0].index)
                                    if next(n for n in self.graph.nodes
                                            if n.index == (e.to_idx if e.from_idx == single_o[0].index
                                                           else e.from_idx)).element == 'H')],
                        pattern_id=str(uuid4())
                    ))
        return reactive_patterns
