import sys
import os
import unittest
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../..")))

from peplab.design.library_design.group_thoeretic.grouptheoreticcomp import (
CyclicPermutationComposition, 
DihedralPermutationComposition
) 


class TestGroupTheoreticComposition(unittest.TestCase):
    """Unit tests for GroupTheoreticComposition subclasses."""

    def setUp(self):
        """Set up instances for testing."""
        self.cyclic_composer = CyclicPermutationComposition()
        self.dihedral_composer = DihedralPermutationComposition()

    def test_cyclic_permutations(self):
        """Test cyclic permutations."""
        result = self.cyclic_composer.generate_composition(['a', 'b', 'c'])
        expected = [['a', 'b', 'c'], ['b', 'c', 'a'], ['c', 'a', 'b']]
        self.assertEqual(result, expected)

    def test_cyclic_single_element(self):
        """Test cyclic permutations with a single element."""
        result = self.cyclic_composer.generate_composition(['a'])
        expected = [['a']]
        self.assertEqual(result, expected)

    def test_cyclic_empty_list(self):
        """Test cyclic permutations with an empty list."""
        result = self.cyclic_composer.generate_composition([])
        expected = []
        self.assertEqual(result, expected)

    def test_dihedral_permutations(self):
        """Test dihedral permutations."""
        result = self.dihedral_composer.generate_composition(['a', 'b', 'c'])
        expected = [
            ['a', 'b', 'c'],
            ['c', 'b', 'a'],
            ['b', 'c', 'a'],
            ['a', 'c', 'b'],
            ['c', 'a', 'b'],
            ['b', 'a', 'c']
        ]
        self.assertCountEqual(result, expected)  # Order doesn't matter for dihedral permutations

    def test_dihedral_single_element(self):
        """Test dihedral permutations with a single element."""
        result = self.dihedral_composer.generate_composition(['a'])
        expected = [['a'], ['a']]  # Single element rotated and reflected
        self.assertEqual(result, expected)

    def test_dihedral_empty_list(self):
        """Test dihedral permutations with an empty list."""
        result = self.dihedral_composer.generate_composition([])
        expected = []
        self.assertEqual(result, expected)


if __name__ == "__main__":
    unittest.main()
