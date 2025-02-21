# peplab/backend/src/domain/models/molecular_graph.py
"""
The molecular graph model.

Classes:
    MolecularGraph: The molecular graph model.
"""

# Standard Library Imports
from typing import Any, Dict, List

# Third Party Imports
from networkx.classes.graph import Graph

# Internal Imports


class MolecularGraph(Graph):
    """The molecular graph model."""

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.G: Graph = Graph(*args, **kwargs)
