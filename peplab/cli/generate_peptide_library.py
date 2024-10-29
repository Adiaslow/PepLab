"""Command-line interface for peptide library generation."""

import argparse

from ..core.library.peptide_library_generator import PeptideLibraryGenerator

parser = argparse.ArgumentParser(
    description="Generate and analyze peptide libraries"
)

parser.add_argument(
    'input_path',
    type=str,
    help="Path to input library file (CSV or JSON)"
)

parser.add_argument(
    '-o', '--output-dir',
    type=str,
    default='output',
    help="Directory for output files"
)

parser.add_argument(
    '--cyclize',
    action='store_true',
    help="Generate cyclic peptides"
)

parser.add_argument(
    '--no-visualization',
    action='store_true',
    help="Skip visualization generation"
)

parser.add_argument(
    '--no-analysis',
    action='store_true',
    help="Skip property analysis"
)

parser.add_argument(
    '--save-intermediates',
    action='store_true',
    help="Save intermediate structures"
)

args = parser.parse_args()

try:
    # Initialize generator
    generator = PeptideLibraryGenerator(
        input_path=args.input_path,
        output_dir=args.output_dir
    )

    # Generate peptides
    peptides = generator.generate_peptides(
        cyclize=args.cyclize,
        save_intermediates=args.save_intermediates
    )

    # Analyze peptides
    if not args.no_analysis:
        generator.analyze_peptides(peptides)

    # Generate visualizations
    if not args.no_visualization:
        generator.visualize_peptides(peptides)

except Exception as e:
    print(f"Error: {str(e)}")
