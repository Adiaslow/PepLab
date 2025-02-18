from dataclasses import dataclass
from typing import Dict

@dataclass
class GraphEdge:
    """Represents a bond in the molecular graph."""
    from_idx: int
    to_idx: int
    bond_type: str
    is_aromatic: bool = False
    is_conjugated: bool = False
    in_ring: bool = False
    stereo: str = 'NONE'

    def to_dict(self) -> Dict:
        """Convert edge to dictionary representation."""
        return {
            'from_idx': self.from_idx,
            'to_idx': self.to_idx,
            'bond_type': self.bond_type,
            'is_aromatic': self.is_aromatic,
            'is_conjugated': self.is_conjugated,
            'in_ring': self.in_ring,
            'stereo': self.stereo
        }
