# Might not have the best understanding of dihedral permutations - trying to learn :)

from rdkit import Chem
from itertools import permutations

def get_core_structure():
    """
    Returns a core structure with defined attachment points
    (e.g., cyclobutane ring) suitable for dihedral symmetry operations.
    """
    return "C1C2CC2C1"  # Example: Cyclobutane ring for D4 symmetry

def get_functional_groups():
    return ["C", "O", "N", "C=O"]  # Methyl, Hydroxyl, Amine, Carbonyl

def apply_dihedral_permutation(core, functional_groups, symmetry_order):
    """
    Generates unique configurations of functional groups around a core 
    structure based on dihedral symmetry (Dₙ).
    
    Parameters:
    - core: SMILES string of the core structure.
    - functional_groups: List of functional groups to attach.
    - symmetry_order: Order of dihedral symmetry (e.g., 4 for D4 symmetry).

    Returns:
    - List of unique SMILES strings with dihedral symmetry applied.
    """
    unique_configurations = set()
    
    # Generate permutations with dihedral symmetry operations
    for perm in permutations(functional_groups):
        for i in range(symmetry_order):
            # Apply rotation symmetry
            rotated_perm = perm[i:] + perm[:i]
            # Add reflection symmetry
            reflected_perm = rotated_perm[::-1]
            
            # Generate SMILES strings and add to unique set
            smiles_rotated = f"{core}." + ".".join(rotated_perm)
            smiles_reflected = f"{core}." + ".".join(reflected_perm)
            
            unique_configurations.add(smiles_rotated)
            unique_configurations.add(smiles_reflected)
    
    return list(unique_configurations)

def main():
    core = get_core_structure()
    functional_groups = get_functional_groups()
    symmetry_order = 4  # Example: D4 symmetry
    
    dihedral_molecules = apply_dihedral_permutation(core, functional_groups, symmetry_order)
    print("Generated dihedral molecules:", dihedral_molecules)

if __name__ == "__main__":
    main()
