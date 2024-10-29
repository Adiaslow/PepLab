from rdkit import Chem
from rdkit.Chem import Descriptors

from property_calculator import PropertyCalculator

class HBACountCalculator(PropertyCalculator):
    """Counts hydrogen bond acceptors."""

    def calculate(self, mol: Chem.Mol) -> int:
        return Descriptors.NumHAcceptors(mol)
