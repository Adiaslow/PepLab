# peplab/chemistry/strategies/factories/reactive_group_factory.py

"""
Module defining the ReactiveGroupFactory class for creating reactive groups.

This module defines the ReactiveGroupFactory class, which provides a method for creating reactive groups based on specified criteria.
"""

from typing import List
from ...core.graph.leaves import Edge, Node
from ..reactivity.reactive_group import ReactiveGroup

class ReactiveGroupFactory:
    """
    Factory class for creating reactive groups.
    """

    @staticmethod
    def create_reactive_group(group_type: str, nodes: List[Node], edges: List[Edge], group_id: str) -> ReactiveGroup:
        """
        Creates a reactive group based on the specified type.

        Args:
            group_type (str): The type of reactive group to create.
            nodes (List[Node]): The nodes that make up the reactive group.
            edges (List[Edge]): The edges that make up the reactive group.
            group_id (str): Unique identifier for the reactive group.

        Returns:
            ReactiveGroup: An instance of ReactiveGroup with the specified type.
        """
        return ReactiveGroup(group_type, nodes, edges, group_id)
