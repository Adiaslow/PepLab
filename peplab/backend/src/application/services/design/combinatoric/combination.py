# peplab/backend/src/application/services/design/combinatoric/combination.py
"""
This module defines the Combination class for generating combinations of items.

Classes:
    Combination: A class for generating combinations of items.
"""

# Standard library imports
from itertools import combinations

# External imports
from typing import List, Any
import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

# Internal imports
from peplab.backend.src.application.interfaces.design.composition import Composition

class Combination:
    @staticmethod
    def generate_composition(items: List[Any], length: int) -> List[List[Any]]:
        """
        Generates all unique combinations of items with a specified length, where order does not matter.

        Parameters:
        - items: List of items to combine (single set)
        - length: Desired length of each combination sequence.

        Returns:
        - List of lists, each representing a unique combination.
        """
        all_combinations = list(combinations(items, length))
        return [list(comb) for comb in all_combinations]


"""
EXAMPLE USAGE
if __name__ == "__main__":
    items = ['a', 'b', 'c']
    combs = Combination.generate_combinations(items, length=2)
    print("Combinations (length 2):", combs)
"""
