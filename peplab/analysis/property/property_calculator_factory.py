from .property_calculator import PropertyCalculator

from .alogp_calculator import ALogPCalculator
from .exact_mass_calculator import ExactMassCalculator
from .rotatable_bonds_calculator import RotatableBondsCalculator
from .hbd_count_calculator import HBDCountCalculator
from .hba_count_calculator import HBACountCalculator

class PropertyCalculatorFactory:
    """Factory for creating property calculators."""

    _calculators = {
        'alogp': ALogPCalculator,
        'exact_mass': ExactMassCalculator,
        'rotatable_bonds': RotatableBondsCalculator,
        'hbd_count': HBDCountCalculator,
        'hba_count': HBACountCalculator
    }

    @classmethod
    def create(cls, property_name: str) -> PropertyCalculator:
        """Creates appropriate calculator for property.

        Args:
            property_name: Name of the property to calculate.

        Returns:
            PropertyCalculator instance.

        Raises:
            ValueError: If property is not supported.
        """
        calculator_class = cls._calculators.get(property_name)
        if calculator_class is None:
            raise ValueError(f"Unsupported property: {property_name}")
        return calculator_class()
