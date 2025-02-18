class CyclicPermutative:
    @staticmethod
    def generate_cyclic_permutations(items, deduplicate=False, reorder=False):
        """
        Generates all cyclic permutations of a list of items by moving each element through every position.

        Parameters:
        - items: List of items to permute.
        - deduplicate (bool): If True, removes any duplicate permutations.
        - reorder (bool): If True, sorts the permutations in lexicographic order.

        Returns:
        - List of lists, each representing a unique cyclic permutation with each element moving through each position.
        """
        n = len(items)
        cyclic_permutations = []
        
        # Generate all cyclic permutations by moving each item through each position
        for moving_index in range(n):
            moving_item = items[moving_index]
            fixed_items = items[:moving_index] + items[moving_index + 1:]
            
            for i in range(n):
                perm = fixed_items[:i] + [moving_item] + fixed_items[i:]
                cyclic_permutations.append(perm)

        # Optional: Deduplicate permutations if deduplicate is True
        if deduplicate:
            cyclic_permutations = CyclicPermutative._deduplicate(cyclic_permutations)

        # Optional: Reorder permutations if reorder is True
        if reorder:
            cyclic_permutations = sorted(cyclic_permutations)

        return cyclic_permutations

    @staticmethod
    def _deduplicate(permutations):
        """
        Helper method to remove duplicate permutations.
        """
        unique_permutations = []
        for perm in permutations:
            if perm not in unique_permutations:
                unique_permutations.append(perm)
        return unique_permutations
'''
EXAMPLE USAGE
if __name__ == "__main__":
    items = ['a', 'b', 'c']
    
    # Generate cyclic permutations without deduplication or reordering
    cyclic_perms = CyclicPermutative.generate_cyclic_permutations(items)
    print("Cyclic Permutations:", cyclic_perms)

    # Generate with deduplication
    cyclic_perms_dedup = CyclicPermutative.generate_cyclic_permutations(items, deduplicate=True)
    print("Deduplicated Cyclic Permutations:", cyclic_perms_dedup)

    # Generate with deduplication and reordering
    cyclic_perms_dedup_sorted = CyclicPermutative.generate_cyclic_permutations(items, deduplicate=True, reorder=True)
    print("Deduplicated and Sorted Cyclic Permutations:", cyclic_perms_dedup_sorted)
'''