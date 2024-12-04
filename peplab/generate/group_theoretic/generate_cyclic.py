import sys
import os

# Add the project root to the Python path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from peplab.design.library_design.group_theoretic.grouptheoreticcomp import CyclicPermutationComposition
from peplab.design.library_design.composer import Composer

def main():
    # Input list
    items = ['a', 'b', 'c', 'd']

    # Initialize CyclicPermutationComposition strategy
    cyclic_strategy = CyclicPermutationComposition()
    composer = Composer(cyclic_strategy)

    # Generate cyclic permutations
    cyclic_permutations = composer.generate_library(items)

    # Display the library
    print("Generated Cyclic Permutations:")
    for permutation in cyclic_permutations:
        print("_".join(permutation))

    # Print total permutations
    print(f"Total Cyclic Permutations: {len(cyclic_permutations)}")

    # Export to CSV
    composer.export_to_csv(cyclic_permutations, "cyclic_permutations.csv")


if __name__ == "__main__":
    main()
