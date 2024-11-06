# @title RDKit Utilities

"""RDKit utility functions for molecular operations.

This module provides helper functions for working with RDKit molecules
and extracting chemical information.
"""

from typing import Dict
from rdkit import Chem

class RDKitUtils:
    @staticmethod
    def get_atom_info(atom: Chem.Atom) -> Dict:
        """Extracts detailed information about an RDKit atom.

        Args:
            atom: RDKit Atom object.

        Returns:
            Dictionary containing atom properties.
        """
        return {
            'atomic_num': atom.GetAtomicNum(),
            'formal_charge': atom.GetFormalCharge(),
            'implicit_valence': atom.GetImplicitValence(),
            'explicit_valence': atom.GetExplicitValence(),
            'aromatic': atom.GetIsAromatic(),
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
    def get_bond_info(bond: Chem.Bond) -> Dict:
        """Extracts detailed information about an RDKit bond.

        Args:
            bond: RDKit Bond object.

        Returns:
            Dictionary containing bond properties.
        """
        return {
            'is_aromatic': bond.GetIsAromatic(),
            'is_conjugated': bond.GetIsConjugated(),
            'in_ring': bond.IsInRing(),
            'stereo': str(bond.GetStereo())
        }

    @staticmethod
    def mol_to_graph_dict(mol: Chem.Mol) -> Dict:
        """Converts an RDKit molecule to a graph dictionary.

        Args:
            mol: RDKit Mol object.

        Returns:
            Dictionary containing graph representation of the molecule.
        """
        if mol is None:
            raise ValueError("Invalid molecule")

        try:
            Chem.SanitizeMol(mol)
            mol = Chem.AddHs(mol)
        except Exception as e:
            raise ValueError(f"Failed to sanitize molecule: {str(e)}")

        nodes = []
        edges = []

        # Add nodes (atoms)
        for atom in mol.GetAtoms():
            nodes.append({
                'id': atom.GetIdx(),
                'element': atom.GetSymbol(),
                **RDKitUtils.get_atom_info(atom)
            })

        # Add edges (bonds)
        for bond in mol.GetBonds():
            edges.append({
                'from_idx': bond.GetBeginAtomIdx(),
                'to_idx': bond.GetEndAtomIdx(),
                'bond_type': str(bond.GetBondType()),
                **RDKitUtils.get_bond_info(bond)
            })

        return {
            'nodes': nodes,
            'edges': edges
        }

    @staticmethod
    def graph_dict_to_mol(graph: Dict) -> Chem.Mol:
        """Converts a graph dictionary back to an RDKit molecule.

        Args:
            graph: Dictionary containing graph representation of molecule.

        Returns:
            RDKit Mol object.

        Raises:
            ValueError: If graph cannot be converted to a valid molecule.
        """
        mol = Chem.RWMol()
        node_to_idx = {}

        # Add atoms
        for node in graph['nodes']:
            if node['element'] != 'H':  # Skip explicit hydrogens
                atom = Chem.Atom(node['element'])
                atom.SetFormalCharge(node['formal_charge'])

                # if node['aromatic']:
                #     atom.SetIsAromatic(True)

                # Set hybridization
                hyb_str = node['hybridization']
                if 'SP3' in hyb_str:
                    atom.SetHybridization(Chem.HybridizationType.SP3)
                elif 'SP2' in hyb_str:
                    atom.SetHybridization(Chem.HybridizationType.SP2)
                elif 'SP' in hyb_str:
                    atom.SetHybridization(Chem.HybridizationType.SP)

                # Set chirality
                if node['chiral']:
                    chiral_tag = node['chiral_tag']
                    if 'CHI_TETRAHEDRAL_CCW' in chiral_tag:
                        atom.SetChiralTag(Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW)
                    elif 'CHI_TETRAHEDRAL_CW' in chiral_tag:
                        atom.SetChiralTag(Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW)

                idx = mol.AddAtom(atom)
                node_to_idx[node['id']] = idx

        # Add bonds
        for edge in graph['edges']:
            # Skip bonds involving hydrogens
            from_node = next(n for n in graph['nodes'] if n['id'] == edge['from_idx'])
            to_node = next(n for n in graph['nodes'] if n['id'] == edge['to_idx'])

            if from_node['element'] != 'H' and to_node['element'] != 'H':
                begin_idx = node_to_idx[edge['from_idx']]
                end_idx = node_to_idx[edge['to_idx']]

                # Set bond type
                if 'SINGLE' in edge['bond_type']:
                    bond_type = Chem.BondType.SINGLE
                elif 'DOUBLE' in edge['bond_type']:
                    bond_type = Chem.BondType.DOUBLE
                elif 'TRIPLE' in edge['bond_type']:
                    bond_type = Chem.BondType.TRIPLE
                elif 'AROMATIC' in edge['bond_type']:
                    bond_type = Chem.BondType.AROMATIC
                else:
                    raise ValueError(f"Unknown bond type: {edge['bond_type']}")

                mol.AddBond(begin_idx, end_idx, bond_type)

                # Set bond properties
                bond = mol.GetBondBetweenAtoms(begin_idx, end_idx)
                """
                if edge['is_aromatic']:
                    bond.SetIsAromatic(True)
                if edge['is_conjugated']:
                    bond.SetIsConjugated(True)
                """
                # Set stereochemistry
                if edge['stereo'] == 'STEREOE':
                    bond.SetStereo(Chem.rdchem.BondStereo.STEREOE)
                elif edge['stereo'] == 'STEREOZ':
                    bond.SetStereo(Chem.rdchem.BondStereo.STEREOZ)

        # Convert to regular molecule and sanitize
        mol = mol.GetMol()
        Chem.SanitizeMol(mol)
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)

        return mol
