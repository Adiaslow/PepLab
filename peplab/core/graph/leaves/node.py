# peplab/core/graph/node.py

"""
Module defining the Node class for representing atoms in a molecular graph.

This module defines the Node class, which represents an atom within a molecular graph.
"""

from typing import Dict
from uuid import uuid4

class Node:
    """
    Represents an atom within a molecular graph.

    Attributes:
        index (int): The index of the atom in the molecular graph.
        element (str): The chemical element of the atom.
        properties (Dict): A dictionary of properties associated with the atom.
        node_id (str): Unique identifier for the node.
    """

    def __init__(self, index: int, element: str, properties: Dict):
        """
        Initializes a Node with the given index, element, and properties.

        Args:
            index (int): The index of the atom in the molecular graph.
            element (str): The chemical element of the atom.
            properties (Dict): A dictionary of properties associated with the atom.
        """
        self.index = index
        self.element = element
        self.properties = properties
        self.node_id = str(uuid4())

    def to_dict(self) -> Dict:
        """
        Converts the node to a dictionary format.

        Returns:
            Dict: Dictionary representation of the node.
        """
        return {
            'index': self.index,
            'element': self.element,
            'properties': self.properties,
            'node_id': self.node_id
        }

    def __str__(self) -> str:
        """
        Returns a string representation of the node.

        Returns:
            str: String representation of the node.
        """
        return f'Node(index={self.index}, element={self.element}, properties={self.properties})'

    def __repr__(self) -> str:
        """
        Returns a detailed string representation of the node for debugging.

        Returns:
            str: Detailed string representation of the node.
        """
        return f'Node(index={self.index}, element={self.element}, properties={self.properties!r}, node_id={self.node_id!r})'
