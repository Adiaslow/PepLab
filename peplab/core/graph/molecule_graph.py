import networkx as nx
from rdkit import Chem
from typing import Dict, Set, Union


class MoleculeGraph:
    """Graph representation of a molecule with reactive site identification."""

    def __init__(self, mol: Union[Chem.Mol, str]):
        if isinstance(mol, str):
            # Handle SMILES input
            self.mol = Chem.MolFromSmiles(mol)
            if self.mol is None:
                raise ValueError(f"Invalid SMILES: {mol}")
        else:
            self.mol = mol

        self.graph = self._create_graph()
        self.reactive_sites: Dict[str, Set[int]] = {}

    @classmethod
    def from_smiles(cls, smiles: str) -> 'MoleculeGraph':
        """Create MoleculeGraph directly from SMILES."""
        return cls(smiles)

    def _create_graph(self) -> nx.Graph:
        graph = nx.Graph()

        for atom in self.mol.GetAtoms():
            graph.add_node(atom.GetIdx(),
                         atomic_num=atom.GetAtomicNum(),
                         formal_charge=atom.GetFormalCharge(),
                         hybridization=atom.GetHybridization(),
                         is_aromatic=atom.GetIsAromatic())

        for bond in self.mol.GetBonds():
            graph.add_edge(bond.GetBeginAtomIdx(),
                         bond.GetEndAtomIdx(),
                         bond_type=bond.GetBondType(),
                         is_aromatic=bond.GetIsAromatic())

        return graph
