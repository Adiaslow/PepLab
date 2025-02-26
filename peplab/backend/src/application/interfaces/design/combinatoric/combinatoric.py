# peplab/backend/src/application/interfaces/design/combinatoric/combinatoric.py
"""
This module defines the CombinatorialComposition class for combinatorial compositions.

Classes:
    CombinatorialComposition: Abstract base class for combinatorial compositions.
"""

# Standard library imports
from itertools import combinations, permutations, product

# External imports
from typing import List, Any

# Internal imports
from peplab.backend.src.application.interfaces.design.composition import Composition


class CombinatorialComposition(Composition):
    """
    Abstract base class for combinatorial compositions.
    """

    def _generate_elements(self, *args, **kwargs):
        raise NotImplementedError("This method must be implemented in a subclass")


class CombinationComposition(CombinatorialComposition):
    """
    Generates combinations of a given length from a list of items.
    """

    def _generate_elements(self, items: List[Any], r: int) -> List[List[Any]]:
        """
        Implements the abstract method to generate combinations.
        """
        return [list(comb) for comb in combinations(items, r)]


class PermutationComposition(CombinatorialComposition):
    """
    Generates permutations of a given length from a list of items.
    """

    def _generate_elements(self, items: List[Any], r: int = None) -> List[List[Any]]:
        """
        Implements the abstract method to generate permutations.
        """
        r = r or len(items)
        return [list(perm) for perm in permutations(items, r)]


class CartesianProductComposition(CombinatorialComposition):
    """
    Generates the Cartesian product of multiple lists.
    """

    def _generate_elements(self, *item_lists: List[List[Any]]) -> List[List[Any]]:
        """
        Implements the abstract method to generate Cartesian products.
        """
        return [list(prod) for prod in product(*item_lists)]


class KFoldCartesianProductComposition(CombinatorialComposition):
    """
    Generates k-fold Cartesian products of a single list.
    """

    def _generate_elements(self, items: List[Any], k: int) -> List[List[Any]]:
        """
        Implements the abstract method to generate k-fold Cartesian products.
        """
        return [list(prod) for prod in product(items, repeat=k)]


from itertools import combinations, permutations, product
from typing import List, Any
from peplab.design.library_design.composition import Composition
from abc import ABC, abstractmethod


class GroupTheoreticComposition(Composition):
    """
    Abstract base class for group-theoretic compositions.
    """

    @abstractmethod
    def _generate_elements(self, *args, **kwargs) -> List[Any]:
        """Abstract method to generate elements for group-theoretic compositions."""
        raise NotImplementedError("This method msut be implemented in a subclass")


class CyclicPermutationComposition(GroupTheoreticComposition):
    """
    Generates all cyclic permutations of a list of items.
    """

    def _generate_elements(self, items: List[Any]) -> List[List[Any]]:
        """
        Generates all cyclic permutations of the input list.
        Each element is shifted to the front in turn.

        Parameters:
        - items (List[Any]): The input list.

        Returns:
        - List[List[Any]]: All cyclic permutations of the input list.
        """
        n = len(items)
        return [items[i:] + items[:i] for i in range(n)]


class DihedralPermutationComposition(GroupTheoreticComposition):
    """
    Generates all dihedral permutations (rotations and reflections) of a list of items.
    """

    def _generate_elements(self, items: List[Any]) -> List[List[Any]]:
        """
        Generates all dihedral permutations of the input list.

        Parameters:
        - items (List[Any]): The input list.

        Returns:
        - List[List[Any]]: All dihedral permutations of the input list.
        """
        n = len(items)
        permutations = []

        # Handle the special case for single-element input
        if n == 1:
            return [items, items[::-1]]

        # Generate rotations
        for i in range(n):
            rotated = items[i:] + items[:i]
            permutations.append(rotated)

            # Add reflection of each rotation
            reflected = rotated[::-1]
            permutations.append(reflected)

        # Remove duplicates
        unique_permutations = [list(x) for x in set(tuple(p) for p in permutations)]
        return unique_permutations
