#!/usr/bin/env python3

import argparse
import json
import logging
import sys
from pathlib import Path

import pandas as pd
from peplab.core.library.peptide_library_generator import PeptideLibraryGenerator

def setup_logging(level=logging.INFO):
    """Configure logging for the script."""
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler('peptide_library_generation.log')
        ]
    )
    return logging.getLogger('peptide-library-generator')

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Generate combinatorial peptide libraries')
    parser.add_argument(
        '--input', '-i',
        required=True,
        type=Path,
        help='Input library file path (CSV or JSON)'
    )
    parser.add_argument(
        '--config', '-c',
        required=True,
        type=Path,
        help='Configuration JSON file path'
    )
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose logging'
    )
    return parser.parse_args()

def load_config(config_path):
    """Load configuration from JSON file."""
    with open(config_path) as f:
        return json.load(f)

def main():
    args = parse_args()
    logger = setup_logging(level=logging.DEBUG if args.verbose else logging.INFO)
    
    try:
        # Load configuration
        logger.info(f"Loading configuration from {args.config}")
        config = load_config(args.config)
        
        # Configure property calculations
        property_config = {
            'alogp': config.get('calc_alogp', True),
            'exact_mass': config.get('calc_exact_mass', True),
            'rotatable_bonds': config.get('calc_rotatable_bonds', True),
            'hbd_count': config.get('calc_hbd_count', True),
            'hba_count': config.get('calc_hba_count', True)
        }
        
        # Initialize generator
        logger.info(f"Initializing peptide library generator with input file: {args.input}")
        generator = PeptideLibraryGenerator(
            input_path=str(args.input),
            property_config=property_config
        )
        
        # Generate peptides
        logger.info("Generating peptide library...")
        peptides = generator.generate_peptides(
            cyclize=(config.get('peptide_type', 'linear') == 'cyclic')
        )
        
        # Analyze peptides
        logger.info("Analyzing peptide properties...")
        results = generator.analyze_peptides(peptides)
        
        # Save results
        output_file = 'library.csv'
        logger.info(f"Saving peptide library to {output_file}")
        results['dataframe'].to_csv(output_file, index=False)
        
        # Log summary
        logger.info("Peptide Library Generation Summary:")
        logger.info(results['summary'])
        
    except Exception as e:
        logger.error(f"Error during peptide library generation: {str(e)}", exc_info=True)
        sys.exit(1)

if __name__ == '__main__':
    main()