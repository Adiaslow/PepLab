from rdkit import Chem
from rdkit.Chem import Descriptors

from .property_calculator import PropertyCalculator

class RotatableBondsCalculator(PropertyCalculator):
    """Counts rotatable bonds."""

    def calculate(self, mol: Chem.Mol) -> int:
        return Descriptors.NumRotatableBonds(mol)
