# peplab/backend/src/application/services/design/combinatoric/permutation.py
"""
This module defines the Permutation class for generating permutations of items.

Classes:
    Permutation: A class for generating permutations of items.
"""

from itertools import permutations
from typing import List, Any


class Permutation:
    @staticmethod
    def generate_composition(items: List[Any], length: int = None, **kwargs) -> List[List[Any]]:
        """
        Generates all possible permutations of items, where order matters.

        Parameters:
        - items: List of items to permute (single set)
        - length: The length of permutations to generate. Defaults to all items if None.

        Returns:
        - List of lists, each representing a unique permutation.
        """
        if length is not None:
            all_permutations = list(permutations(items, r=length))
        else:
            all_permutations = list(permutations(items))
        return [list(perm) for perm in all_permutations]


"""
EXAMPLE USAGE
if __name__ == "__main__":
    items = ['a', 'b', 'c']
    perms = Permutation.generate_permutations(items)
    print("Permutations:", perms)
"""
