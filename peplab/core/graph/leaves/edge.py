# peplab/core/graph/edge.py

"""
Module defining the Edge class for representing bonds in a molecular graph.

This module defines the Edge class, which represents a bond within a molecular graph.
"""

from typing import Dict
from uuid import uuid4

class Edge:
    """
    Represents a bond within a molecular graph.

    Attributes:
        from_idx (int): The index of the starting atom of the bond.
        to_idx (int): The index of the ending atom of the bond.
        properties (Dict): A dictionary of properties associated with the bond.
        edge_id (str): Unique identifier for the edge.
    """

    def __init__(self, from_idx: int, to_idx: int, properties: Dict):
        """
        Initializes an Edge with the given starting and ending atom indices and properties.

        Args:
            from_idx (int): The index of the starting atom of the bond.
            to_idx (int): The index of the ending atom of the bond.
            properties (Dict): A dictionary of properties associated with the bond.
        """
        self.from_idx = from_idx
        self.to_idx = to_idx
        self.properties = properties
        self.edge_id = str(uuid4())

    def to_dict(self) -> Dict:
        """
        Converts the edge to a dictionary format.

        Returns:
            Dict: Dictionary representation of the edge.
        """
        return {
            'from_idx': self.from_idx,
            'to_idx': self.to_idx,
            'properties': self.properties,
            'edge_id': self.edge_id
        }

    def __str__(self) -> str:
        """
        Returns a string representation of the edge.

        Returns:
            str: String representation of the edge.
        """
        return f'Edge(from_idx={self.from_idx}, to_idx={self.to_idx}, properties={self.properties})'

    def __repr__(self) -> str:
        """
        Returns a detailed string representation of the edge for debugging.

        Returns:
            str: Detailed string representation of the edge.
        """
        return f'Edge(from_idx={self.from_idx}, to_idx={self.to_idx}, properties={self.properties!r}, edge_id={self.edge_id!r})'
