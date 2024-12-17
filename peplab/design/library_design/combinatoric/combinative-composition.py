from itertools import combinations, permutations, product
from typing import List, Any
from composition import Composition
from abc import ABC, abstractmethod

class CombinatorialComposition(Composition):
    """
    Abstract base class for combinatorial compositions.
    """

    @abstractmethod
    def generate_composition(self, *args, **kwargs):
        raise NotImplementedError("This method must be implemented in a subclass")


class CombinationComposition(CombinatorialComposition):
    """
    Generates combinations of a given length from a list of items.
    """

    def generate_composition(self, items: List[Any], r: int) -> List[List[Any]]:
        """
        Generates all combinations of length r from the input list.

        Parameters:
        - items (List[Any]): The input list.
        - r (int): The length of the combinations.

        Returns:
        - List[List[Any]]: All combinations of length r.
        """
        return [list(comb) for comb in combinations(items, r)]


class PermutationComposition(CombinatorialComposition):
    """
    Generates permutations of a given length from a list of items.
    """

    def generate_composition(self, items: List[Any], r: int = None) -> List[List[Any]]:
        """
        Generates all permutations of length r from the input list.
        If r is None, permutations will be of the same length as the input list.

        Parameters:
        - items (List[Any]): The input list.
        - r (int, optional): The length of the permutations. Defaults to the length of items.

        Returns:
        - List[List[Any]]: All permutations of the given length.
        """
        r = r or len(items)
        return [list(perm) for perm in permutations(items, r)]


class CartesianProductComposition(CombinatorialComposition):
    """
    Generates the Cartesian product of multiple lists.
    """

    def generate_composition(self, *item_lists: List[List[Any]]) -> List[List[Any]]:
        """
        Generates the Cartesian product of the input lists.

        Parameters:
        - item_lists (List[List[Any]]): A variable number of input lists.

        Returns:
        - List[List[Any]]: The Cartesian product of the input lists.
        """
        return [list(prod) for prod in product(*item_lists)]


class KFoldCartesianProductComposition(CombinatorialComposition):
    """
    Generates k-fold Cartesian products of a single list.
    """

    def generate_composition(self, items: List[Any], k: int) -> List[List[Any]]:
        """
        Generates k-fold Cartesian products of the input list.

        Parameters:
        - items (List[Any]): The input list.
        - k (int): The number of repetitions for the Cartesian product.

        Returns:
        - List[List[Any]]: The k-fold Cartesian product.
        """
        return [list(prod) for prod in product(items, repeat=k)]
