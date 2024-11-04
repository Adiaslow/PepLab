class CyclicPermutative:
    @staticmethod
    def generate_cyclic_permutations(items):
        """
        Generates all cyclic permutations of a list of items.

        Parameters:
        - items: List of items to permute in a cyclic manner.

        Returns:
        - List of lists, each representing a unique cyclic permutation.
        """
        n = len(items)
        return [items[i:] + items[:i] for i in range(n)]
'''
EXAMPLE USAGE
if __name__ == "__main__":
    items = ['a', 'b', 'c']
    cyclic_perms = CyclicPermutative.generate_cyclic_permutations(items)
    print("Cyclic Permutations:", cyclic_perms)
'''