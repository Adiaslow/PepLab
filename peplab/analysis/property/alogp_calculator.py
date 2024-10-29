from rdkit import Chem
from rdkit.Chem import Descriptors

from .property_calculator import PropertyCalculator

class ALogPCalculator(PropertyCalculator):
    """Calculates atomic logP."""

    def calculate(self, mol: Chem.Mol) -> float:
        return Descriptors.MolLogP(mol)
