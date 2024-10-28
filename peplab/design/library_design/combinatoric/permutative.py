from itertools import permutations
from typing import List, Callable, Any, Optional, Dict


# One feature that can be included is Filter Functionality: Allows custom filtering with filter_func, 
# so you can limit permutations based on specific criteria.
class Permutative:
    def __init__(self, items: List[Any]):
        """
        Initializes the Permutative class with a list of items to be permuted.

        Parameters:
        - items: A list containing items to permute.
        """
        self.items = items

    def generate_permutations(self, length: Optional[int] = None) -> List[List[Any]]:
        """
        Generates all possible permutations of items.

        Parameters:
        - length: Optional integer specifying the length of each permutation. 
                  If not provided, the length of `self.items` will be used.

        Returns:
        - List of lists, where each sublist is a valid permutation of items.
        """
        # Default length to the length of items if not specified
        length = length or len(self.items)
        
        # Generate all permutations of the specified length
        all_permutations = list(permutations(self.items, length))
        
        # Return all permutations without filtering
        return [list(perm) for perm in all_permutations]

    def to_dict(self, permutations: List[List[Any]]) -> Dict[int, List[Any]]:
        """
        Converts a list of permutations into a dictionary where keys are indices and values are permutations.

        Parameters:
        - permutations: List of permutations to convert.

        Returns:
        - Dictionary with index keys and permutation values.
        """
        return {i: perm for i, perm in enumerate(permutations)}


'''
EXAMPLE USAGE
perm = Permutative(['a', 'b', 'c'])
all_perms = perm.generate_permutations(length=2)
print("All Permutations (length 2):", all_perms)
perm_dict = perm.to_dict(all_perms)
print("Permutations Dictionary:", perm_dict
'''
