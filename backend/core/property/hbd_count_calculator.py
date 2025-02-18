from rdkit import Chem
from rdkit.Chem import Descriptors

from .property_calculator import PropertyCalculator

class HBDCountCalculator(PropertyCalculator):
    """Counts hydrogen bond donors."""

    def calculate(self, mol: Chem.Mol) -> int:
        return Descriptors.NumHDonors(mol)
