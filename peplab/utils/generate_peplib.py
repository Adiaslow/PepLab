from peplab.design.library_design.combinatoric.combinative_composition import CartesianProductComposition
from peplab.design.library_design.composer import Composer

def main():
    # Define input items (e.g., amino acids)
    amino_acids = [
        ['Leu', 'Phe', 'Val', 'Ala'],
        ['DLeu', 'DPhe', 'DVal', 'DAla'],
        ['LeuMe', 'PheMe', 'ValMe', 'AlaMe'],
        ['DLeuMe', 'DPheMe', 'DValMe', 'DAlaMe']
    ]

    # Initialize Composer with CartesianProductComposition
    cartesian_strategy = CartesianProductComposition()
    composer = Composer(cartesian_strategy)

    # Generate the peptide library
    peptide_library = composer.generate_library(*amino_acids)

    # Display the library
    print("Generated Peptide Library:")
    for peptide in peptide_library[:10]:  # Limit output for readability
        print("_".join(peptide))

    # Print total combinations
    print(f"Total Combinations: {len(peptide_library)}")

    # Optional: Export the library to a CSV file
    composer.export_to_csv(peptide_library, "peptide_library.csv")


if __name__ == "__main__":
    main()
