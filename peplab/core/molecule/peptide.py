from typing import List, Optional, Union

from ..base.molecular_entity import MolecularEntity
from ..graph.molecule_graph import MoleculeGraph
from .residue import Residue

class Peptide(MolecularEntity):
    """Represents a peptide molecule, either built from residues or loaded directly."""

    def __init__(self,
                 source: Union[str, List['Residue']],
                 name: Optional[str] = None):
        super().__init__()
        self.name = name or f"Peptide_{self.id[:8]}"

        if isinstance(source, str):
            # Create from SMILES
            self.molecule = MoleculeGraph.from_smiles(source)
            self.sequence = None  # Could be determined through analysis
            self._residues = None
        else:
            # Build from residues
            self._residues = tuple(source)
            self.sequence = ''.join(res.code for res in source)
            # Combine residue graphs into single molecule
            self.molecule = self._combine_residues()

        self.properties.set_property("sequence", self.sequence, str)

    def _combine_residues(self) -> MoleculeGraph:
        # Implementation for combining residue graphs
        pass

    @classmethod
    def from_smiles(cls, smiles: str, name: Optional[str] = None) -> 'Peptide':
        """Create peptide directly from SMILES string."""
        return cls(source=smiles, name=name)

    def get_type(self) -> str:
        return "Peptide"
