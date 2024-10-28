# combinative.py

from rdkit import Chem
from rdkit.Chem import Descriptors
import csv

# Step 1: Define Core Structures and Functional Groups
def get_core_structures():
    return [
        "C1=CC=CC=C1",  # Benzene
        "C1=CC=NC=C1"   # Pyridine
    ]

def get_functional_groups():
    return [
        "C",    # Methyl
        "O",    # Hydroxyl
        "N",    # Amine
        "C=O"   # Carbonyl
    ]

# Generate Combinations, Calculate Properties, and Store in Dictionary
def generate_combinations():
    cores = get_core_structures()
    functional_groups = get_functional_groups()
    
    # Dictionary to store combinations with metadata
    molecules = {}
    
    for core in cores:
        molecules[core] = []  # Initialize a list for each core structure

        for group in functional_groups:
            # Combine core and functional group (this is a simplified combination)
            smiles = f"{core}.{group}"
            
            # Convert SMILES to RDKit molecule
            mol = Chem.MolFromSmiles(smiles)
            
            # Calculate molecular weight
            mol_weight = Descriptors.MolWt(mol) if mol else None
            
            # Store each combination as a dictionary with metadata
            molecule_data = {
                "smiles": smiles,
                "functional_group": group,
                "molecular_weight": mol_weight
            }
            
            # Add to list of molecules for this core structure
            molecules[core].append(molecule_data)
    
    return molecules

# Output Results to CSV
def save_to_csv(molecule_dict, filename="output_molecules.csv"):
    """
    Saves the generated combinations with metadata to a CSV file.
    Each row in the CSV represents a molecule with its core, functional group, SMILES, and properties.
    """
    with open(filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        # Write the header row
        writer.writerow(["Core Structure", "Functional Group", "SMILES", "Molecular Weight"])

        # Write each molecule's data
        for core, molecules in molecule_dict.items():
            for molecule in molecules:
                writer.writerow([
                    core,
                    molecule["functional_group"],
                    molecule["smiles"],
                    molecule["molecular_weight"]
                ])

# Main Function to Run the Entire Process
def main():
    # Generate the dictionary with combinations and metadata
    molecule_dict = generate_combinations()
    
    # Save results to CSV
    save_to_csv(molecule_dict, filename="output_molecules.csv")
    print("Combinations generated and saved to output_molecules.csv")

if __name__ == "__main__":
    main() # Tentative driver code to be adapted as the framework continues
