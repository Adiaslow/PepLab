# backend/core/properties/property_calculator.py

from rdkit import Chem
from rdkit.Chem import Descriptors

def calculate_properties(sequence):
    """Calculates molecular properties using RDKit."""
    mol = Chem.MolFromSequence(sequence)
    return {
        "mass": Descriptors.ExactMolWt(mol),
        "logP": Descriptors.MolLogP(mol),
        "hbd": Descriptors.NumHDonors(mol),
        "hba": Descriptors.NumHAcceptors(mol),
    }

class RDKitCalculator:
    """RDKit-based property calculator."""
    def calculate(self, peptide):
        return calculate_properties(peptide.sequence)
