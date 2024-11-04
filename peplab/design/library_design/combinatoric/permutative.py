from itertools import permutations
from typing import List, Any

class Permutation:
    @staticmethod
    def generate_permutations(items: List[Any]) -> List[List[Any]]:
        """
        Generates all possible permutations of items, where order matters.

        Parameters:
        - items: List of items to permute.

        Returns:
        - List of lists, each representing a unique permutation.
        """
        all_permutations = list(permutations(items))
        return [list(perm) for perm in all_permutations]

'''
EXAMPLE USAGE
if __name__ == "__main__":
    items = ['a', 'b', 'c']
    perms = Permutation.generate_permutations(items)
    print("Permutations:", perms)
'''