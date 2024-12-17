from dataclasses import dataclass
from typing import FrozenSet

from ...analysis.property.reactivity_profile import ReactivityProfile
from .reactive_type import ReactiveType

@dataclass(frozen=True)
class ReactiveSite:
    """Represents a reactive site within a molecule."""
    atoms: FrozenSet[int]  # atom indices
    reactivity_type: ReactiveType
    profile: ReactivityProfile

    def __hash__(self):
        return hash((frozenset(self.atoms), self.reactivity_type))
