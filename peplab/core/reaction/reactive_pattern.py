# peplab/core/reaction/reactive_pattern.py

"""
Module defining the ReactivePattern class.

This module defines the ReactivePattern class, which represents a reactive pattern identified in a molecular graph.
"""

from typing import List
from uuid import uuid4
from ..molecule.atom import GraphNode
from ..molecule.bond import GraphEdge
from .reactive_type import ReactiveType

class ReactivePattern:
    """
    Represents a reactive pattern identified in a molecular graph.

    Attributes:
        pattern_id (str): Unique identifier for the reactive pattern.
        type (ReactiveType): Type of the reactive pattern.
        atoms (List[GraphNode]): List of atoms involved in the reactive pattern.
        bonds (List[GraphEdge]): List of bonds involved in the reactive pattern.
    """

    def __init__(self, type: ReactiveType, atoms: List[GraphNode], bonds: List[GraphEdge], pattern_id: str = None):
        """
        Initializes a ReactivePattern with specified type, atoms, and bonds.

        Args:
            type (ReactiveType): Type of the reactive pattern.
            atoms (List[GraphNode]): List of atoms involved in the reactive pattern.
            bonds (List[GraphEdge]): List of bonds involved in the reactive pattern.
            pattern_id (str, optional): Unique identifier for the reactive pattern. Defaults to a new UUID.
        """
        self.pattern_id = pattern_id or str(uuid4())
        self.type = type
        self.atoms = atoms
        self.bonds = bonds

    def to_dict(self) -> dict:
        """
        Converts the reactive pattern to a dictionary format.

        Returns:
            dict: Dictionary representation of the reactive pattern.
        """
        return {
            'pattern_id': self.pattern_id,
            'type': self.type.name,
            'atoms': [atom.to_dict() for atom in self.atoms],
            'bonds': [bond.to_dict() for bond in self.bonds]
        }
