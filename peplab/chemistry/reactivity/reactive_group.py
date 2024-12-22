# peplab/chemistry/strategies/reactivity/reactive_group.py

"""
Module defining the ReactiveGroup class for representing reactive groups within a molecular graph.

This module defines the ReactiveGroup class, which provides common functionality for all reactive groups.
"""

from typing import List, Dict
from ...core.graph.leaves import Edge, Node

class ReactiveGroup:
    """
    Represents a reactive group within a molecular graph.

    Attributes:
        name (str): The name of the reactive group.
        nodes (List[Node]): The nodes that make up the reactive group.
        edges (List[Edge]): The edges that make up the reactive group.
        group_id (str): Unique identifier for the reactive group.
    """

    def __init__(self, name: str, nodes: List[Node], edges: List[Edge], group_id: str):
        """
        Initializes a ReactiveGroup with the given name, nodes, edges, and group identifier.

        Args:
            name (str): The name of the reactive group.
            nodes (List[Node]): The nodes that make up the reactive group.
            edges (List[Edge]): The edges that make up the reactive group.
            group_id (str): Unique identifier for the reactive group.
        """
        self.name = name
        self.nodes = nodes
        self.edges = edges
        self.group_id = group_id

    def to_dict(self) -> Dict:
        """
        Converts the reactive group to a dictionary format.

        Returns:
            Dict: Dictionary representation of the reactive group.
        """
        return {
            'name': self.name,
            'nodes': [node.to_dict() for node in self.nodes],
            'edges': [edge.to_dict() for edge in self.edges],
            'group_id': self.group_id
        }
