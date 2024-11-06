from typing import Dict, Optional
import math

from .thermodynamics import ThermodynamicParameters
from ...core.reaction.reaction_conditions import ReactionConditions

class ReactivityProfile:
    """Defines chemical reactivity characteristics."""

    def __init__(self,
                 thermodynamics: ThermodynamicParameters,
                 electron_density: float,
                 orbital_energies: Optional[Dict[str, float]] = None):
        self.thermodynamics = thermodynamics
        self.electron_density = electron_density
        self.orbital_energies = orbital_energies or {}

    def get_reactivity_score(self, conditions: 'ReactionConditions') -> float:
        rate_constant = self._calculate_rate_constant(conditions)
        free_energy = self.thermodynamics.get_free_energy(conditions.temperature)
        return rate_constant * self.electron_density * math.exp(-free_energy / (8.314 * conditions.temperature))

    def _calculate_rate_constant(self, conditions) -> float:
        return 0.0
