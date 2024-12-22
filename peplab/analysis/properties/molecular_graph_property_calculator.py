# peplab/analysis/property/molecular_graph_property_calculator.py

"""
Module for calculating properties of atoms and bonds from RDKit molecules.
"""

from rdkit import Chem

class MolecularGraphPropertyCalculator:
    """
    Class for calculating properties of atoms and bonds from RDKit molecules.
    """

    @staticmethod
    def get_atom_properties(atom: Chem.Atom) -> dict:
        """
        Extracts properties of an atom.

        Args:
            atom (Chem.Atom): RDKit atom object.

        Returns:
            dict: Dictionary of atom properties.
        """
        return {
            'atomic_num': atom.GetAtomicNum(),
            'formal_charge': atom.GetFormalCharge(),
            'implicit_valence': atom.GetImplicitValence(),
            'explicit_valence': atom.GetExplicitValence(),
            'is_aromatic': atom.GetIsAromatic(),
            'hybridization': str(atom.GetHybridization()),
            'num_explicit_hs': atom.GetNumExplicitHs(),
            'num_implicit_hs': atom.GetNumImplicitHs(),
            'total_num_hs': atom.GetTotalNumHs(),
            'degree': atom.GetDegree(),
            'in_ring': atom.IsInRing(),
            'chiral': atom.GetChiralTag() != Chem.rdchem.ChiralType.CHI_UNSPECIFIED,
            'chiral_tag': str(atom.GetChiralTag())
        }

    @staticmethod
    def get_bond_properties(bond: Chem.Bond) -> dict:
        """
        Extracts properties of a bond.

        Args:
            bond (Chem.Bond): RDKit bond object.

        Returns:
            dict: Dictionary of bond properties.
        """
        return {
            'bond_type': str(bond.GetBondType()),
            'is_aromatic': bond.GetIsAromatic(),
            'is_conjugated': bond.GetIsConjugated(),
            'in_ring': bond.IsInRing(),
            'stereo': str(bond.GetStereo())
        }
