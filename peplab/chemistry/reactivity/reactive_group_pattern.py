# peplab/chemistry/reactivity/reactive_group_pattern.py

"""
Module defining the ReactiveGroupPattern class for reactive group patterns.

This module defines the ReactiveGroupPattern class, which provides a structure for defining reactive group patterns.
"""

from typing import List, Dict, Optional

class ReactiveGroupPattern:
    def __init__(self, name: str, node_criteria: List[Dict], edge_criteria: Optional[List[Dict]] = None):
        """
        Initializes a ReactiveGroupPattern.

        Args:
            name (str): The name of the reactive group pattern.
            node_criteria (List[Dict]): A list of dictionaries defining the criteria for nodes.
            edge_criteria (List[Dict], optional): A list of dictionaries defining the criteria for edges. Defaults to None.
        """
        self.name = name
        self.node_criteria = node_criteria
        self.edge_criteria = edge_criteria or []

    def __repr__(self):
        return f"ReactiveGroupPattern(name={self.name}, node_criteria={self.node_criteria}, edge_criteria={self.edge_criteria})"
