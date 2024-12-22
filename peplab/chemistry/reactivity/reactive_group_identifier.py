# peplab/chemistry/reactivity/reactive_group_identifier.py

from typing import List, Dict, Set, Optional
from uuid import uuid4

from peplab.core.graph.composites import MolecularGraph
from peplab.core.graph.leaves import Node, Edge
from .reactive_group import ReactiveGroup
from .reactive_group_pattern_registry import pattern_registry
from peplab.utils.logging import Logger

class ReactiveGroupIdentifier:
    """Identifies reactive groups in a molecular graph."""

    def __init__(self, graph: MolecularGraph):
        self.graph = graph
        self.patterns = pattern_registry.get_all_patterns()
        self.logger = Logger().get_logger()

    def identify_reactive_groups(self) -> List[ReactiveGroup]:
        """Identifies all reactive groups in the molecular graph."""
        reactive_groups = []

        # First identify aromatic rings to mark aromatic nodes
        aromatic_rings = self._identify_aromatic_rings()
        self.logger.debug(f"Aromatic rings found: {len(aromatic_rings)}")
        for ring in aromatic_rings:
            for node_idx in ring:
                for node in self.graph.nodes:
                    if node.index == node_idx:
                        node.properties['is_aromatic'] = True
                        self.logger.debug(f"Node {node.index} marked as aromatic")
                        break

        # Then identify all reactive groups
        for node in self.graph.nodes:
            patterns = pattern_registry.get_patterns_for_element(node.element)
            self.logger.debug(f"Node {node.index} ({node.element}) has {len(patterns)} patterns")

            for pattern in patterns:
                for criteria in pattern.node_criteria:
                    if self._node_matches_criteria(node, criteria):
                        match = self._build_match_from_node(node, pattern)
                        if match:
                            group = ReactiveGroup(
                                name=pattern.name,
                                nodes=match['nodes'],
                                edges=match['edges'],
                                group_id=f"{pattern.name}_{uuid4().hex[:8]}"
                            )
                            reactive_groups.append(group)
                            self.logger.debug(f"Reactive group {group.name} identified")
                            break

        self.logger.info(f"Total reactive groups identified: {len(reactive_groups)}")
        return reactive_groups

    def _node_matches_criteria(self, node: Node, criteria: Dict) -> bool:
        """Checks if a node matches given criteria."""
        self.logger.debug(f"Checking node {node.index} against criteria: {criteria}")

        # Check element
        if node.element != criteria['element']:
            return False

        # Check neighbors count
        if 'neighbors' in criteria:
            neighbors = self._get_node_neighbors(node)
            if len(neighbors) != criteria['neighbors']:
                return False

        # Check hydrogen count
        if 'hydrogen_count' in criteria:
            h_count = sum(1 for n, _ in self._get_node_neighbors(node) if n.element == 'H')
            if h_count != criteria['hydrogen_count']:
                return False

        # Check aromatic system
        if criteria.get('check_aromatic_system', False):
            if not node.properties.get('is_aromatic', False):
                return False

        # Check for specific bond types to elements
        if 'single_bond' in criteria:
            if not any(n.element == criteria['single_bond'] and e.properties.get('bond_type') == 'SINGLE'
                       for n, e in self._get_node_neighbors(node)):
                return False

        if 'double_bond' in criteria:
            if not any(n.element == criteria['double_bond'] and e.properties.get('bond_type') == 'DOUBLE'
                       for n, e in self._get_node_neighbors(node)):
                return False

        if 'triple_bond' in criteria:
            if not any(n.element == criteria['triple_bond'] and e.properties.get('bond_type') == 'TRIPLE'
                       for n, e in self._get_node_neighbors(node)):
                return False

        return True

    def _build_match_from_node(self, start_node: Node, pattern) -> Optional[Dict]:
        """Builds a complete match starting from a node."""
        self.logger.debug(f"Building match from node {start_node.index} for pattern: {pattern.name}")
        matched_nodes = [start_node]
        matched_edges = []

        if pattern.edge_criteria:
            edges = self._find_matching_edges(start_node, pattern.edge_criteria)
            if edges:
                # Add connected nodes for appropriate patterns
                for edge in edges:
                    other_idx = edge.to_idx if edge.from_idx == start_node.index else edge.from_idx
                    other_node = next(n for n in self.graph.nodes if n.index == other_idx)
                    if other_node.element != 'H':  # Only include non-hydrogen atoms
                        matched_nodes.append(other_node)
                matched_edges.extend(edges)
                return {'nodes': matched_nodes, 'edges': matched_edges}

        return None

    def _find_matching_edges(self, node: Node, edge_criteria: List[Dict]) -> List[Edge]:
        """Finds edges that match given criteria."""
        neighbors = self._get_node_neighbors(node)
        matching_edges = []

        # Convert single edge criteria to list for uniform handling
        if not isinstance(edge_criteria, list):
            edge_criteria = [edge_criteria]

        for criteria in edge_criteria:
            # Handle multiple allowed bond types
            if isinstance(criteria['bond_type'], list):
                allowed_types = criteria['bond_type']
                required_count = criteria.get('count', 1)
                edges = []
                for neighbor, edge in neighbors:
                    if edge.properties.get('bond_type') in allowed_types:
                        edges.append(edge)
                if len(edges) == required_count:
                    matching_edges.extend(edges)
                else:
                    return []
            else:
                # Handle single bond type with element requirements
                found = False
                for neighbor, edge in neighbors:
                    if (edge.properties.get('bond_type') == criteria['bond_type'] and
                        ((edge.from_idx == node.index and
                          neighbor.element == criteria.get('to_element')) or
                         (edge.to_idx == node.index and
                          neighbor.element == criteria.get('from_element')))):
                        matching_edges.append(edge)
                        found = True
                        break
                if not found:
                    return []

        return matching_edges

    def _get_node_neighbors(self, node: Node) -> List[tuple[Node, Edge]]:
        """Gets all neighboring nodes and their connecting edges."""
        neighbors = []
        for edge in self.graph.edges:
            if edge.from_idx == node.index:
                neighbor = next(n for n in self.graph.nodes if n.index == edge.to_idx)
                neighbors.append((neighbor, edge))
            elif edge.to_idx == node.index:
                neighbor = next(n for n in self.graph.nodes if n.index == edge.from_idx)
                neighbors.append((neighbor, edge))
        return neighbors

    def _identify_aromatic_rings(self) -> List[Set[int]]:
        """Identifies all aromatic rings in the molecule."""
        rings = []
        visited = set()

        for node in self.graph.nodes:
            if node.index not in visited and node.element == 'C' and node.properties.get('is_aromatic', False):
                ring = self._find_aromatic_ring(node.index, visited)
                if ring:
                    rings.append(ring)
                    visited.update(ring)

        return rings

    def _find_aromatic_ring(self, start_idx: int, visited: Set[int], path: List[int] = None) -> Optional[Set[int]]:
        """Finds an aromatic ring starting from a given carbon atom."""
        if path is None:
            path = []

        path = path + [start_idx]
        visited.add(start_idx)

        if len(path) == 6:
            if start_idx == path[0] and self._check_alternating_bonds(path):
                return set(path)
            return None

        if len(path) > 6:
            return None

        neighbors = []
        for edge in self.graph.edges:
            other_idx = None
            if edge.from_idx == start_idx and edge.to_idx not in visited:
                other_idx = edge.to_idx
            elif edge.to_idx == start_idx and edge.from_idx not in visited:
                other_idx = edge.from_idx

            if other_idx is not None:
                other_node = next(n for n in self.graph.nodes if n.index == other_idx)
                if (other_node.element == 'C' and
                    edge.properties.get('bond_type') in ['SINGLE', 'DOUBLE']):
                    neighbors.append(other_idx)

        for next_idx in neighbors:
            result = self._find_aromatic_ring(next_idx, visited, path)
            if result:
                return result

        return None

    def _check_alternating_bonds(self, path: List[int]) -> bool:
        """Checks if a ring has alternating single/double bonds."""
        bonds = []
        for i in range(len(path)):
            j = (i + 1) % len(path)
            edge = self._find_edge(path[i], path[j])
            if not edge:
                return False
            bonds.append(edge.properties.get('bond_type'))

        expected_patterns = [
            ['SINGLE', 'DOUBLE', 'SINGLE', 'DOUBLE', 'SINGLE', 'DOUBLE'],
            ['DOUBLE', 'SINGLE', 'DOUBLE', 'SINGLE', 'DOUBLE', 'SINGLE']
        ]
        return bonds in expected_patterns

    def _find_edge(self, from_idx: int, to_idx: int) -> Optional[Edge]:
        """Finds an edge between two nodes."""
        for edge in self.graph.edges:
            if ((edge.from_idx == from_idx and edge.to_idx == to_idx) or
                (edge.from_idx == to_idx and edge.to_idx == from_idx)):
                return edge
        return None
