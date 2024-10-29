import logging
from typing import List, Tuple, Optional, Dict
from rdkit import Chem

from ..molecule.atom import GraphNode
from ..molecule.bond import GraphEdge

class MolecularGraph:
    """Graph representation of a molecule."""

    def __init__(self):
        """Initialize empty molecular graph."""
        self.nodes: List[GraphNode] = []
        self.edges: List[GraphEdge] = []
        self.logger = logging.getLogger(self.__class__.__name__)

    @classmethod
    def from_smiles(cls, smiles: str) -> 'MolecularGraph':
        """Create molecular graph from SMILES with improved error handling."""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError(f"Failed to parse SMILES: {smiles}")

            Chem.SanitizeMol(mol)
            mol = Chem.AddHs(mol)

            graph = cls()

            # Add atoms with properties from old code
            for atom in mol.GetAtoms():
                properties = {
                    'atomic_num': atom.GetAtomicNum(),
                    'formal_charge': atom.GetFormalCharge(),
                    'hybridization': str(atom.GetHybridization()),
                    'implicit_valence': atom.GetImplicitValence(),
                    'explicit_valence': atom.GetExplicitValence(),
                    'aromatic': atom.GetIsAromatic(),
                    'num_explicit_hs': atom.GetNumExplicitHs(),
                    'num_implicit_hs': atom.GetNumImplicitHs(),
                    'total_num_hs': atom.GetTotalNumHs(),
                    'degree': atom.GetDegree(),
                    'in_ring': atom.IsInRing(),
                    'chiral': atom.GetChiralTag() != Chem.rdchem.ChiralType.CHI_UNSPECIFIED,
                    'chiral_tag': str(atom.GetChiralTag())
                }

                node = GraphNode(
                    id=atom.GetIdx(),
                    element=atom.GetSymbol(),
                    atomic_num=properties['atomic_num'],
                    formal_charge=properties['formal_charge'],
                    implicit_valence=properties['implicit_valence'],
                    explicit_valence=properties['explicit_valence'],
                    aromatic=properties['aromatic'],
                    hybridization=properties['hybridization'],
                    num_explicit_hs=properties['num_explicit_hs'],
                    num_implicit_hs=properties['num_implicit_hs'],
                    total_num_hs=properties['total_num_hs'],
                    degree=properties['degree'],
                    in_ring=properties['in_ring'],
                    chiral=properties['chiral'],
                    chiral_tag=properties['chiral_tag']
                )
                graph.nodes.append(node)

            # Add bonds with properties from old code
            for bond in mol.GetBonds():
                edge = GraphEdge(
                    from_idx=bond.GetBeginAtomIdx(),
                    to_idx=bond.GetEndAtomIdx(),
                    bond_type=str(bond.GetBondType()),
                    is_aromatic=bond.GetIsAromatic(),
                    is_conjugated=bond.GetIsConjugated(),
                    in_ring=bond.IsInRing(),
                    stereo=str(bond.GetStereo())
                )
                graph.edges.append(edge)

            return graph

        except Exception as e:
            raise ValueError(f"Error creating molecular graph from SMILES: {str(e)}")

    def get_neighbors(self, node_id: int) -> List[Tuple[GraphNode, GraphEdge]]:
        """Get all neighbors of a node with connecting edges."""
        neighbors = []
        for edge in self.edges:
            if edge.from_idx == node_id:
                neighbor = next(n for n in self.nodes if n.id == edge.to_idx)
                neighbors.append((neighbor, edge))
            elif edge.to_idx == node_id:
                neighbor = next(n for n in self.nodes if n.id == edge.from_idx)
                neighbors.append((neighbor, edge))
        return neighbors

    def find_reactive_sites(self, nuc_pattern: str, elec_pattern: str) -> None:
        """Find and mark reactive sites matching old code logic exactly."""
        if nuc_pattern == 'NH2':
            self._find_nh2_pattern()
        elif nuc_pattern == 'NH':
            self._find_nh_pattern()

        if elec_pattern == 'COOH':
            self._find_cooh_pattern()

    def _find_nh2_pattern(self) -> None:
        """Find NH2 groups matching old code exactly."""
        for node in self.nodes:
            if node.element == 'N':
                neighbors = self.get_neighbors(node.id)
                h_count = sum(1 for n, _ in neighbors if n.element == 'H')
                heavy_count = sum(1 for n, _ in neighbors if n.element != 'H')

                if h_count == 2 and heavy_count == 1:
                    node.is_reactive_nuc = True

    def _find_nh_pattern(self) -> None:
        """Find NH groups matching old code exactly."""
        for node in self.nodes:
            if node.element == 'N':
                neighbors = self.get_neighbors(node.id)
                h_count = sum(1 for n, _ in neighbors if n.element == 'H')
                has_carbon = any(n.element == 'C' for n, _ in neighbors)

                if h_count == 1 and has_carbon:
                    node.is_reactive_nuc = True

    def _find_cooh_pattern(self) -> None:
        """Find COOH groups matching old code exactly."""
        for node in self.nodes:
            if node.element == 'C':
                neighbors = self.get_neighbors(node.id)
                double_o = False
                single_oh = False
                other_heavy = False

                for neighbor, edge in neighbors:
                    if neighbor.element == 'O':
                        if edge.bond_type == 'DOUBLE':
                            double_o = True
                        elif edge.bond_type == 'SINGLE':
                            o_neighbors = self.get_neighbors(neighbor.id)
                            h_count = sum(1 for n, _ in o_neighbors if n.element == 'H')
                            if h_count == 1:
                                single_oh = True
                    elif neighbor.element != 'H':
                        if edge.bond_type == 'SINGLE':
                            other_heavy = True

                if double_o and single_oh and other_heavy:
                    node.is_reactive_elec = True

    def to_dict(self) -> Dict:
        """Convert to dictionary format matching old code."""
        return {
            'nodes': [
                {
                    'id': n.id,
                    'element': n.element,
                    'atomic_num': n.atomic_num,
                    'formal_charge': n.formal_charge,
                    'implicit_valence': n.implicit_valence,
                    'explicit_valence': n.explicit_valence,
                    'aromatic': n.aromatic,
                    'hybridization': n.hybridization,
                    'num_explicit_hs': n.num_explicit_hs,
                    'num_implicit_hs': n.num_implicit_hs,
                    'total_num_hs': n.total_num_hs,
                    'degree': n.degree,
                    'in_ring': n.in_ring,
                    'chiral': n.chiral,
                    'chiral_tag': n.chiral_tag,
                    'is_reactive_nuc': n.is_reactive_nuc,
                    'is_reactive_elec': n.is_reactive_elec
                }
                for n in self.nodes
            ],
            'edges': [
                {
                    'from_idx': e.from_idx,
                    'to_idx': e.to_idx,
                    'bond_type': e.bond_type,
                    'is_aromatic': e.is_aromatic,
                    'is_conjugated': e.is_conjugated,
                    'in_ring': e.in_ring,
                    'stereo': e.stereo
                }
                for e in self.edges
            ]
        }
