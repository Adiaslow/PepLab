from typing import Any, Dict, Optional

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
