from dataclasses import dataclass
from typing import Dict

from .residue import ResidueInfo

@dataclass
class PositionInfo:
    """Information about a position in the peptide sequence.

    Attributes:
        residues: Dictionary mapping residue keys to ResidueInfo objects.
        temperature: Optional temperature parameter for the position.
    """
    residues: Dict[str, ResidueInfo]
    temperature: float = 20.0
