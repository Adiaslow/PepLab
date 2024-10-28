from typing import Dict
from peptide_library import PeptideLibrary

class LibraryParser:
    """Parser for different types of library input formats."""

    @staticmethod
    def parse_library(input_data: Dict) -> PeptideLibrary:
        """Parse library from input dictionary."""
        if "smiles_peptides" in input_data:
            # Direct SMILES peptide library
            return PeptideLibrary.from_smiles_dict(
                input_data["name"],
                input_data["smiles_peptides"]
            )
        else:
            # Residue-based library
            return PeptideLibrary.from_residue_library(
                input_data["name"],
                input_data
            )
