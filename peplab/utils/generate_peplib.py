from peplab.design.library_design.combinatoric.combinative_composition import CombinationComposition
from peplab.design.library_design.composer import Composer

# Define input items (e.g., amino acids)
amino_acids = ['A', 'R', 'N', 'D']

# Initialize Composer with CombinationComposition
combination_strategy = CombinationComposition()
composer = Composer(combination_strategy)

# Generate a peptide library of length 3
peptide_library = composer.generate_library(amino_acids, r=3)

# Display the library
print("Generated Peptide Library:")
for peptide in peptide_library:
    print(''.join(peptide))

# Export the library to a CSV file
composer.export_to_csv(peptide_library, "peptide_library.csv")
