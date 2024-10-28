# combinative.py

from rdkit import Chem
from rdkit.Chem import Descriptors
import csv

# 1. Define Building Blocks
def get_core_structures():
    """
    Returns a list of core structures (molecular scaffolds).
    These will serve as the base to attach functional groups.
    """
    cores = [
        '''
        [ insert sample here ]
        '''
    ]
    return cores

def get_functional_groups():
    """
    Returns a list of functional groups or substituents that can be attached to core structures.
    """
    functional_groups = [
        "C",     # Methyl
        "O",     # Hydroxyl
        "N",     # Amine
        "C=O"    # Carbonyl
    ]
    return functional_groups

# 2. Generate Combinations
def combine_fragments():
    """
    Generates combinations of core structures with functional groups.
    """
    cores = get_core_structures()
    functional_groups = get_functional_groups()

    # Using a dictionary to store combinations
    molecules = {}
    for core in cores:
        molecules[core] = []
        for group in functional_groups:
            smiles = f"{core}.{group}"
            mol = Chem.MolFromSmiles(smiles)
            
            # Calculate properties
            mol_weight = Descriptors.MolWt(mol) if mol else None
            
            # Store as a dictionary with additional information
            molecule_data = {
                "smiles": smiles,
                "functional_group": group,
                "molecular_weight": mol_weight
            }
            
            molecules[core].append(molecule_data)


    return molecules

# 3. Filtering Functionality
def filter_by_molecular_weight(molecule_smiles, max_weight=500):
    """
    Filters a molecule by its molecular weight.
    """
    mol = Chem.MolFromSmiles(molecule_smiles)
    if mol:
        mol_weight = Descriptors.MolWt(mol)
        return mol_weight <= max_weight
    return False

def filter_molecules(molecule_list, max_weight=500):
    """
    Filters a list of molecules based on molecular weight.
    """
    filtered_molecules = []
    for smiles in molecule_list:
        if filter_by_molecular_weight(smiles, max_weight):
            filtered_molecules.append(smiles)
    return filtered_molecules

# 4. Outputting Results
def save_to_csv(molecule_list, filename="output_molecules.csv"):
    """
    Saves a list of molecules (SMILES strings) to a CSV file.
    """
    with open(filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["SMILES"])
        for molecule in molecule_list:
            writer.writerow([molecule])

# Main Function to Bring It All Together
def main():
    # Step 1: Generate combinations
    molecules = combine_fragments()

    # Step 2: Filter molecules by molecular weight
    filtered_molecules = filter_molecules(molecules, max_weight=500)

    # Step 3: Save the filtered molecules to a CSV file
    save_to_csv(filtered_molecules, filename="filtered_molecules.csv")
    print(f"Generated and saved {len(filtered_molecules)} molecules.")

if __name__ == "__main__":
    main()
    # driver code, tentative - will be updated as we continue rolling the framework
