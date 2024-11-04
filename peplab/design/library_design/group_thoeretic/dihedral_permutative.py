class DihedralPermutative:
    @staticmethod
    def generate_dihedral_permutations(items):
        """
        Generates all dihedral permutations of a list of items (both rotations and reflections).

        Parameters:
        - items: List of items to permute with dihedral symmetry.

        Returns:
        - List of lists, each representing a unique dihedral permutation.
        """
        n = len(items)
        permutations = []
        
        # Generate rotations and add them to permutations
        for i in range(n):
            rotated = items[i:] + items[:i]
            permutations.append(rotated)
            
            # Add the reflection of each rotation
            reflected = rotated[::-1]
            permutations.append(reflected)
        
        # Remove duplicates by converting to set of tuples, then back to list of lists
        unique_permutations = [list(x) for x in set(tuple(p) for p in permutations)]
        return unique_permutations
'''
EXAMPLE USAGE
if __name__ == "__main__":
    items = ['a', 'b', 'c']
    dihedral_perms = DihedralPermutative.generate_dihedral_permutations(items)
    print("Dihedral Permutations:", dihedral_perms)
'''