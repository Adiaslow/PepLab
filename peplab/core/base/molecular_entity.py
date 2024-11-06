from abc import ABC, abstractmethod
from typing import Optional
import uuid

from .property_store import PropertyStore

class MolecularEntity(ABC):
    """Base class for all molecular entities."""

    def __init__(self, id: Optional[str] = None):
        self.id = id if id is not None else str(uuid.uuid4())
        self.properties = PropertyStore()

    @abstractmethod
    def get_type(self) -> str:
        pass

    def __hash__(self):
        return hash(self.id)
