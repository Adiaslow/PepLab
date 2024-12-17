from dataclasses import dataclass
from typing import Dict, List, Optional

@dataclass
class PeptideInfo:
    """Information about a generated peptide."""
    graph: Dict
    id: str
    sequence: List[str]
    peptide_type: str
    smiles: Optional[str] = None
