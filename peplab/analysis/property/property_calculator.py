from abc import ABC, abstractmethod
from rdkit import Chem

class PropertyCalculator(ABC):
    """Abstract base class for property calculators."""

    @abstractmethod
    def calculate(self, mol: Chem.Mol) -> float:
        """Calculates property for a molecule."""
        pass
