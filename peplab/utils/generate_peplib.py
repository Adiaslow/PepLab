from peplab.design.library_design.combinatoric.combinative_composition import CartesianProductComposition
from peplab.design.library_design.composer import Composer

def main():
    amino_acids = [
        ['Leu', 'Phe', 'Val', 'Ala'],
        ['DLeu', 'DPhe', 'DVal', 'DAla'],
        ['LeuMe', 'PheMe', 'ValMe', 'AlaMe'],
        ['DLeuMe', 'DPheMe', 'DValMe', 'DAlaMe']
    ]

    cartesian_strategy = CartesianProductComposition()
    composer = Composer(cartesian_strategy)

    peptide_library = composer.generate_library(*amino_acids)

    print("Generated Peptide Library:")
    for peptide in peptide_library[:5]:
        print(peptide)
    print(f"Total Combinations: {len(peptide_library)}")
    composer.export_to_csv(peptide_library, "peptide_library.csv")


if __name__ == "__main__":
    main()
