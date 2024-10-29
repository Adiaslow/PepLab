from typing import Any, Dict, Optional
from dataclasses import dataclass

class PropertyStore:
    """Flexible key-value store for molecular properties with type hints."""

    def __init__(self):
        self._properties: Dict[str, Any] = {}
        self._property_types: Dict[str, type] = {}

    def set_property(self, key: str, value: Any, property_type: Optional[type] = None) -> None:
        if property_type is not None and not isinstance(value, property_type):
            raise TypeError(f"Value {value} is not of type {property_type}")
        self._properties[key] = value
        if property_type:
            self._property_types[key] = property_type

    def get_property(self, key: str) -> Any:
        return self._properties.get(key)

    def __hash__(self):
        return hash(frozenset(self._properties.items()))

@dataclass
class PropertyConfig:
    """Configuration for property calculations.

    Attributes:
        alogp: Calculate atomic logP.
        exact_mass: Calculate exact molecular mass.
        rotatable_bonds: Count rotatable bonds.
        hbd_count: Count hydrogen bond donors.
        hba_count: Count hydrogen bond acceptors.
    """
    alogp: bool = True
    exact_mass: bool = True
    rotatable_bonds: bool = True
    hbd_count: bool = True
    hba_count: bool = True

    @classmethod
    def from_dict(cls, config_dict: Dict[str, bool]) -> 'PropertyConfig':
        """Creates PropertyConfig from dictionary."""
        return cls(**{
            k: v for k, v in config_dict.items()
            if k in cls.__dataclass_fields__
        })
