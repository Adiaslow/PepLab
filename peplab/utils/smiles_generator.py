from typing import Dict, Optional, Tuple
from rdkit import Chem

from peplab.utils.rdkit_utils import RDKitUtils

class SMILESGenerator:
    """Generates SMILES strings from molecular graphs."""

    @staticmethod
    def generate(
        graph_dict: Dict,
        canonical: bool = True
    ) -> Tuple[Optional[str], bool]:
        """Generates SMILES string from graph dictionary.

        Args:
            graph_dict: Dictionary containing molecular graph data.
            canonical: Whether to return canonical SMILES.

        Returns:
            Tuple of (SMILES string or None, success boolean).
        """

        try:
            mol = RDKitUtils.graph_dict_to_mol(graph_dict)
            if mol is None:
                return None, False

            # Generate SMILES
            smiles = Chem.MolToSmiles(mol, canonical=canonical)

            # Verify SMILES is valid
            test_mol = Chem.MolFromSmiles(smiles)
            if test_mol is None:
                return None, False

            return smiles, True

        except Exception as e:
            print(f"Error generating SMILES: {str(e)}")
            return None, False
