# peplab/core/graph/molecule_graph.py

"""
Module for representing and manipulating a molecular graph.

This module defines the MolecularGraph class, which represents a molecular graph consisting of nodes (atoms)
and edges (bonds). It provides methods to get neighbors of a node and convert the graph to a dictionary format.
"""

import logging
from typing import List, Tuple, Dict
from ..leaves import Edge, Node


class MolecularGraph:
    """
    Represents a molecular graph consisting of nodes (atoms) and edges (bonds).

    Attributes:
        nodes (List[Node]): List of nodes (atoms) in the molecular graph.
        edges (List[Edge]): List of edges (bonds) in the molecular graph.
        logger (logging.Logger): Logger for logging information and errors.
    """

    def __init__(self):
        """
        Initializes the MolecularGraph with empty lists of nodes and edges, and sets up logging.
        """
        self.nodes: List[Node] = []
        self.edges: List[Edge] = []
        self.logger = logging.getLogger(self.__class__.__name__)

    def add_node(self, node: Node) -> None:
        """
        Adds a node to the molecular graph.

        Args:
            node (Node): The node to add.
        """
        self.nodes.append(node)

    def add_edge(self, edge: Edge) -> None:
        """
        Adds an edge to the molecular graph.

        Args:
            edge (Edge): The edge to add.
        """
        self.edges.append(edge)

    def get_neighbors(self, node_index: int) -> List[Tuple[Node, Edge]]:
        """
        Gets all neighbors of a node along with the connecting edges using node index.

        Args:
            node_index (int): The index of the node whose neighbors are to be found.

        Returns:
            List[Tuple[Node, Edge]]: List of tuples containing neighboring nodes and the connecting edges.
        """
        return [(n, e) for e in self.edges
                if (e.from_idx == node_index and (n := next(n for n in self.nodes if n.index == e.to_idx)))
                or (e.to_idx == node_index and (n := next(n for n in self.nodes if n.index == e.from_idx)))]

    def to_dict(self) -> Dict:
        """
        Converts the molecular graph to a dictionary format.

        Returns:
            Dict: Dictionary representation of the molecular graph.
        """
        return {
            'nodes': [n.to_dict() for n in self.nodes],
            'edges': [e.to_dict() for e in self.edges]
        }

    def __str__(self) -> str:
        """
        Returns a string representation of the molecular graph.

        Returns:
            str: String representation of the molecular graph.
        """
        node_str = ', '.join([str(node) for node in self.nodes])
        edge_str = ', '.join([str(edge) for edge in self.edges])
        return f'MolecularGraph(nodes=[{node_str}], edges=[{edge_str}])'

    def __repr__(self) -> str:
        """
        Returns a detailed string representation of the molecular graph for debugging.

        Returns:
            str: Detailed string representation of the molecular graph.
        """
        return f'MolecularGraph(nodes={self.nodes!r}, edges={self.edges!r})'
