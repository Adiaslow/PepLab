# cyclic_permutative.py

from rdkit import Chem
from itertools import permutations

def get_core_structure():
    """
    Returns a core structure with defined attachment points
    (e.g., benzene ring) suitable for cyclic symmetry operations.
    """
    return "C1=CC=CC=C1"  # Example: Benzene ring for C6 symmetry

def get_functional_groups():
    return ["C", "O", "N"]  # Methyl, Hydroxyl, Amine

def apply_cyclic_permutation(core, functional_groups, symmetry_order):
    """
    Generates unique configurations of functional groups around a core 
    structure based on a cyclic symmetry (Cₙ).
    
    Parameters:
    - core: SMILES string of the core structure.
    - functional_groups: List of functional groups to attach.
    - symmetry_order: Order of cyclic symmetry (e.g., 3 for C3 symmetry).

    Returns:
    - List of unique SMILES strings with cyclic symmetry applied.
    """
    unique_configurations = set()
    
    # Generate permutations
    for perm in permutations(functional_groups):
        # Generate cyclic configurations
        for i in range(symmetry_order):
            rotated_perm = perm[i:] + perm[:i]
            smiles = f"{core}." + ".".join(rotated_perm)
            unique_configurations.add(smiles)
    
    return list(unique_configurations)

def main():
    core = get_core_structure()
    functional_groups = get_functional_groups()
    symmetry_order = 3  # Example: C3 symmetry
    
    cyclic_molecules = apply_cyclic_permutation(core, functional_groups, symmetry_order)
    print("Generated cyclic molecules:", cyclic_molecules)

if __name__ == "__main__":
    main()
    # More driver code to be adapted to framework as we progress