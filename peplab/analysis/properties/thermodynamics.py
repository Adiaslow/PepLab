from dataclasses import dataclass
from typing import Optional

@dataclass(frozen=True)
class ThermodynamicParameters:
    """Thermodynamic and kinetic parameters for reactions."""
    activation_energy: float
    enthalpy: float
    entropy: float
    pKa: Optional[float] = None

    def get_free_energy(self, temperature: float) -> float:
        return self.enthalpy - (temperature * self.entropy / 1000)
