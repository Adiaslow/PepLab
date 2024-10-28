from typing import Dict, Set
from rdkit.Chem import Descriptors

from ..base.property_store import PropertyStore
from ..molecule import Peptide

class PeptideLibrary:
    """Represents a collection of peptides with analysis capabilities."""

    def __init__(self, name: str):
        self.name = name
        self._peptides: Dict[str, Peptide] = {}
        self.properties = PropertyStore()

    @classmethod
    def from_smiles_dict(cls, name: str, smiles_dict: Dict[str, str]) -> 'PeptideLibrary':
        """Create library directly from dictionary of SMILES strings."""
        library = cls(name)
        for peptide_name, smiles in smiles_dict.items():
            peptide = Peptide.from_smiles(smiles, name=peptide_name)
            library.add_peptide(peptide)
        return library

    @classmethod
    def from_residue_library(cls, name: str, residue_library: Dict) -> 'PeptideLibrary':
        """Create library by combining residues according to specifications."""
        library = cls(name)
        # Implementation for building peptides from residues
        return library

    def add_peptide(self, peptide: Peptide) -> None:
        self._peptides[peptide.name] = peptide

    def analyze_properties(self, property_set: Set[str]) -> Dict[str, Dict[str, float]]:
        """Calculate specified properties for all peptides."""
        results = {}
        for name, peptide in self._peptides.items():
            results[name] = self._calculate_properties(peptide, property_set)
        return results

    def _calculate_properties(self, peptide: Peptide,
                            property_set: Set[str]) -> Dict[str, float]:
        """Calculate molecular properties for a peptide."""
        results = {}
        mol = peptide.molecule.mol

        for prop in property_set:
            if prop == "exact_mass":
                results[prop] = Descriptors.ExactMolWt(mol)
            elif prop == "alogp":
                results[prop] = Descriptors.MolLogP(mol)
            elif prop == "rotatable_bonds":
                results[prop] = Descriptors.NumRotatableBonds(mol)
            elif prop == "hbd_count":
                results[prop] = Descriptors.NumHDonors(mol)
            elif prop == "hba_count":
                results[prop] = Descriptors.NumHAcceptors(mol)

        return results
