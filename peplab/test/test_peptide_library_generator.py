from pathlib import Path

from peplab.core.library.peptide_library_generator import PeptideLibraryGenerator

def run_test(input_path: str, output_dir: Path) -> int:
    """Example usage of the peptide library generator."""
    # Set up paths
    current_dir = Path.cwd()

    # Initialize generator with corrected visualization config
    generator = PeptideLibraryGenerator(
        input_path=input_path,
        output_dir=output_dir,
        property_config={
            'alogp': True,
            'exact_mass': True,
            'rotatable_bonds': True,
            'hbd_count': True,
            'hba_count': True
        },
        visualization_config={
            'size': (400, 400),  # Changed from mol_size to size
            'fmt': 'svg',
            'graph_size': (10, 10)  # This will be handled via kwargs
        }
    )

    try:
        # Generate cyclic peptides
        peptides = generator.generate_peptides(
            cyclize=True,
            save_intermediates=True
        )

        # Analyze peptides
        results = generator.analyze_peptides(peptides)

        # Generate visualizations
        generator.visualize_peptides(peptides)

        # Print summary
        print("\nAnalysis Summary:")
        print(results['summary'])

    except Exception as e:
        print(f"Error: {str(e)}")
        return 1

    return 0
