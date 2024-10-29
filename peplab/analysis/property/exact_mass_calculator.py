from rdkit import Chem
from rdkit.Chem import Descriptors

from property_calculator import PropertyCalculator

class ExactMassCalculator(PropertyCalculator):
    """Calculates exact molecular mass."""

    def calculate(self, mol: Chem.Mol) -> float:
        return Descriptors.ExactMolWt(mol)
