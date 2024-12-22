from dataclasses import dataclass
from typing import FrozenSet

from .reactive_type import ReactiveType

@dataclass(frozen=True)
class ReactiveSite:
    """Represents a reactive site within a molecule."""
    atoms: FrozenSet[int]  # atom indices
    reactivity_type: ReactiveType

    def __hash__(self):
        return hash((frozenset(self.atoms), self.reactivity_type))
