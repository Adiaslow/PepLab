import sys
import os

# Add the project root to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../Peplab")))

from peplab.design.library_design.group_theoretic.grouptheoreticcomp import DihedralPermutationComposition
from peplab.design.library_design.composer import Composer

def main():
    # Input list
    items = ['a', 'b', 'c', 'd']

    # Initialize DihedralPermutationComposition strategy
    dihedral_strategy = DihedralPermutationComposition()
    composer = Composer(dihedral_strategy)

    # Generate dihedral permutations
    dihedral_permutations = composer.generate_library(items)

    # Display the library
    print("Generated Dihedral Permutations:")
    for permutation in dihedral_permutations:
        print("_".join(permutation))

    # Print total permutations
    print(f"Total Dihedral Permutations: {len(dihedral_permutations)}")

    # Export to CSV
    composer.export_to_csv(dihedral_permutations, "dihedral_permutations.csv")


if __name__ == "__main__":
    main()
