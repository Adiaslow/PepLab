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
        
        # Generate all rotations
        for i in range(n):
            rotated = items[i:] + items[:i]
            permutations.append(rotated)
            
            # Generate reflection of the rotated list
            reflected = rotated[::-1]
            permutations.append(reflected)
        
        return permutations

# Example usage
if __name__ == "__main__":
    items = [1, 2, 3, 4]
    dihedral_perms = DihedralPermutative.generate_dihedral_permutations(items)
    print("Dihedral Permutations:", dihedral_perms)
