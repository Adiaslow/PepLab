from peplab.core.library.peptide_library_generator import PeptideLibraryGenerator

input_file = "peplab/test/ppLibInputAzideTest1.csv" # @param {type:"string"}
#@markdown Select peptide type:
peptide_type = "cyclic" #@param ["linear", "cyclic"]

#@markdown Configure property calculations:
calc_alogp = True #@param {type:"boolean"}
calc_exact_mass = True #@param {type:"boolean"}
calc_rotatable_bonds = True #@param {type:"boolean"}
calc_hbd_count = True #@param {type:"boolean"}
calc_hba_count = True #@param {type:"boolean"}

#@markdown Output filename:
output_file = "peptide_library.csv" #@param {type:"string"}

# Configure property calculations
property_config = {
    'alogp': calc_alogp,
    'exact_mass': calc_exact_mass,
    'rotatable_bonds': calc_rotatable_bonds,
    'hbd_count': calc_hbd_count,
    'hba_count': calc_hba_count
}

try:
    # Initialize generator
    generator = PeptideLibraryGenerator(
        input_path=input_file,
        property_config=property_config
    )

    # Generate peptides
    peptides = generator.generate_peptides(
        cyclize=(peptide_type == "cyclic")
    )

    # Analyze peptides
    results = generator.analyze_peptides(peptides)

    # Save results
    results['dataframe'].to_csv(output_file, index=False)

    # Display summary
    print("\nAnalysis Summary:")
    print(results['summary'])

except Exception as e:
    print(f"Error: {str(e)}")
