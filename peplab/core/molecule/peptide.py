@dataclass
class PeptideInfo:
    """Information about a generated peptide."""
    graph: Dict
    sequence: List[str]
    peptide_type: str
    smiles: Optional[str] = None
