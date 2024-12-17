# peplab/molecule/residue.py
from dataclasses import dataclass

@dataclass
class ResidueInfo:
    """Information about a single residue."""
    id: str
    name: str
    smiles: str
    nucleophile: str
    electrophile: str

    @classmethod
    def from_dict(cls, data: dict) -> 'ResidueInfo':
        return cls(
            id=['id'],
            name=data['name'],  # Use 'id' from dict as name
            smiles=data['smiles'],
            nucleophile=data['nuc'],
            electrophile=data['elec']
        )

    def to_dict(self) -> dict:
        """Convert to dictionary format matching expected structure."""
        return {
            'id': self.name,  # Use name as 'id' in dict
            'smiles': self.smiles,
            'nuc': self.nucleophile,
            'elec': self.electrophile
        }
