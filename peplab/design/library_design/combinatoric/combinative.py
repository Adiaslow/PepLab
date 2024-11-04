from itertools import combinations
from typing import List, Any

class Combination:
    @staticmethod
    def generate_combinations(items: List[Any], length: int) -> List[List[Any]]:
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

'''
EXAMPLE USAGE
if __name__ == "__main__":
    items = ['a', 'b', 'c']
    combs = Combination.generate_combinations(items, length=2)
    print("Combinations (length 2):", combs)
'''
