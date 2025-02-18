from itertools import combinations, permutations, product
from typing import List, Any
from peplab.design.library_design.composition import Composition
from abc import ABC, abstractmethod


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
