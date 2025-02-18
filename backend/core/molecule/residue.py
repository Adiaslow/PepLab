from dataclasses import dataclass

@dataclass
class ResidueInfo:
    """Information about a single residue."""
    name: str
    smiles: str
    nucleophile: str
    electrophile: str
