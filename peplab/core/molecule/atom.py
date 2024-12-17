from dataclasses import dataclass
from typing import Dict

@dataclass
class GraphNode:
    """Represents an atom in the molecular graph."""
    id: str
    index: int
    element: str
    atomic_num: int
    formal_charge: int
    implicit_valence: int
    explicit_valence: int
    is_aromatic: bool
    hybridization: str
    num_explicit_hs: int
    num_implicit_hs: int
    total_num_hs: int
    degree: int
    in_ring: bool
    chiral: bool = False
    chiral_tag: str = 'CHI_UNSPECIFIED'
    is_reactive_nuc: bool = False
    is_reactive_elec: bool = False
    is_reactive_click: bool = False

    def to_dict(self) -> Dict:
        """Convert node to dictionary representation."""
        return {
            'id': self.id,
            'index': self.index,
            'element': self.element,
            'atomic_num': self.atomic_num,
            'formal_charge': self.formal_charge,
            'implicit_valence': self.implicit_valence,
            'explicit_valence': self.explicit_valence,
            'is_aromatic': self.is_aromatic,
            'hybridization': self.hybridization,
            'num_explicit_hs': self.num_explicit_hs,
            'num_implicit_hs': self.num_implicit_hs,
            'total_num_hs': self.total_num_hs,
            'degree': self.degree,
            'in_ring': self.in_ring,
            'chiral': self.chiral,
            'chiral_tag': self.chiral_tag,
            'is_reactive_nuc': self.is_reactive_nuc,
            'is_reactive_elec': self.is_reactive_elec,
            'is_reactive_click': self.is_reactive_click
        }

    def __hash__(self):
        return hash(self.id)
