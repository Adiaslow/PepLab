import sys
import os
import unittest
# Dynamically add the root directory to sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../PepLab")))

from peplab.design.library_design.combinatoric.combinative_composition import (
    CombinationComposition,
    PermutationComposition,
    CartesianProductComposition,
    KFoldCartesianProductComposition,
)


class TestCombinatorialComposition(unittest.TestCase):
    """Unit tests for CombinatorialComposition subclasses."""

    def test_combination_composition(self):
        """Test CombinationComposition."""
        comb = CombinationComposition()

        # Test standard combinations
        result = comb.generate_composition([1, 2, 3], r=2)
        expected = [[1, 2], [1, 3], [2, 3]]
        self.assertEqual(result, expected)

        # Test edge case: r=0
        result = comb.generate_composition([1, 2, 3], r=0)
        expected = [[]]  # Only one combination: the empty combination
        self.assertEqual(result, expected)

        # Test edge case: r > len(items)
        result = comb.generate_composition([1, 2], r=3)
        expected = []  # No combinations possible
        self.assertEqual(result, expected)

    def test_permutation_composition(self):
        """Test PermutationComposition."""
        perm = PermutationComposition()

        # Test standard permutations
        result = perm.generate_composition([1, 2, 3], r=2)
        expected = [
            [1, 2],
            [1, 3],
            [2, 1],
            [2, 3],
            [3, 1],
            [3, 2],
        ]
        self.assertEqual(result, expected)

        # Test edge case: r=None (full-length permutations)
        result = perm.generate_composition([1, 2, 3])
        expected = [
            [1, 2, 3],
            [1, 3, 2],
            [2, 1, 3],
            [2, 3, 1],
            [3, 1, 2],
            [3, 2, 1],
        ]
        self.assertEqual(result, expected)

        # Test edge case: r=1
        result = perm.generate_composition([1, 2, 3], r=1)
        expected = [[1], [2], [3]]
        self.assertEqual(result, expected)

    def test_cartesian_product_composition(self):
        """Test CartesianProductComposition."""
        cart = CartesianProductComposition()

        # Test standard Cartesian product
        result = cart.generate_composition([1, 2], ['a', 'b'])
        expected = [[1, 'a'], [1, 'b'], [2, 'a'], [2, 'b']]
        self.assertEqual(result, expected)

        # Test edge case: One empty list
        result = cart.generate_composition([], ['a', 'b'])
        expected = []  # Cartesian product with an empty list is empty
        self.assertEqual(result, expected)

        # Test edge case: Multiple empty lists
        result = cart.generate_composition([], [])
        expected = []  # No combinations possible
        self.assertEqual(result, expected)

    def test_k_fold_cartesian_product_composition(self):
        """Test KFoldCartesianProductComposition."""
        kfold = KFoldCartesianProductComposition()

        # Test standard k-fold Cartesian product
        result = kfold.generate_composition([1, 2], k=2)
        expected = [[1, 1], [1, 2], [2, 1], [2, 2]]
        self.assertEqual(result, expected)

        # Test edge case: k=0
        result = kfold.generate_composition([1, 2], k=0)
        expected = [[]]  # Only one combination: the empty combination
        self.assertEqual(result, expected)

        # Test edge case: Empty list
        result = kfold.generate_composition([], k=2)
        expected = []  # No combinations possible
        self.assertEqual(result, expected)


if __name__ == "__main__":
    unittest.main()
