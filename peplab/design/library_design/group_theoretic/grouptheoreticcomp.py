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
