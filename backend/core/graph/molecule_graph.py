import logging
from typing import List, Tuple, Optional, Dict
from rdkit import Chem

from ..molecule.atom import GraphNode
from ..molecule.bond import GraphEdge

class MolecularGraph:
    def __init__(self):
        self.nodes: List[GraphNode] = []
        self.edges: List[GraphEdge] = []
        self.logger = logging.getLogger(self.__class__.__name__)

    @staticmethod
    def from_smiles(smiles: str) -> 'MolecularGraph':
        """Create molecular graph from SMILES."""
        graph = MolecularGraph()
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError(f"Failed to parse SMILES: {smiles}")

            Chem.SanitizeMol(mol)
            mol = Chem.AddHs(mol)

            # Add atoms
            for atom in mol.GetAtoms():
                node = GraphNode(
                    id=atom.GetIdx(),
                    element=atom.GetSymbol(),
                    atomic_num=atom.GetAtomicNum(),
                    formal_charge=atom.GetFormalCharge(),
                    implicit_valence=atom.GetImplicitValence(),
                    explicit_valence=atom.GetExplicitValence(),
                    aromatic=atom.GetIsAromatic(),
                    hybridization=str(atom.GetHybridization()),
                    num_explicit_hs=atom.GetNumExplicitHs(),
                    num_implicit_hs=atom.GetNumImplicitHs(),
                    total_num_hs=atom.GetTotalNumHs(),
                    degree=atom.GetDegree(),
                    in_ring=atom.IsInRing(),
                    chiral=atom.GetChiralTag() != Chem.rdchem.ChiralType.CHI_UNSPECIFIED,
                    chiral_tag=str(atom.GetChiralTag())
                )
                graph.nodes.append(node)

            # Add bonds
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
        """Find and mark reactive sites."""
        if nuc_pattern == 'NH2':
            self._find_nh2_pattern()
        elif nuc_pattern == 'NH':
            self._find_nh_pattern()
        elif nuc_pattern == 'N3':
            self._find_azide_pattern()

        if elec_pattern == 'COOH':
            self._find_cooh_pattern()
        elif elec_pattern == 'C#C':
            self._find_alkyne_pattern()

    def _find_azide_pattern(self) -> None:
        """Find azide (N3) group and mark first nitrogen as reactive."""
        for node in self.nodes:
            if node.element == 'N':
                # Check for N-N≡N pattern
                neighbors = self.get_neighbors(node.id)
                n2 = None
                for neighbor, edge in neighbors:
                    if neighbor.element == 'N':
                        n2 = neighbor
                        break

                if n2:
                    n2_neighbors = self.get_neighbors(n2.id)
                    n3 = None
                    for neighbor, edge in n2_neighbors:
                        if neighbor.element == 'N' and neighbor.id != node.id:
                            n3 = neighbor
                            break

                    if n3:
                        # Verify N3 pattern: N-N≡N
                        n3_neighbors = self.get_neighbors(n3.id)
                        if (len(neighbors) == 2 and  # First N has 2 connections
                            len(n2_neighbors) == 2 and  # Middle N has 2 connections
                            len(n3_neighbors) == 1):  # Terminal N has 1 connection
                            node.is_reactive_nuc = True

    def _find_alkyne_pattern(self) -> None:
        """Find terminal alkyne (C≡C) group and mark as reactive."""
        for node in self.nodes:
            if node.element == 'C':
                neighbors = self.get_neighbors(node.id)
                for neighbor, edge in neighbors:
                    if (neighbor.element == 'C' and
                        edge.bond_type == 'TRIPLE'):
                        # Check if either carbon is terminal (has only one other connection)
                        neighbor_connections = self.get_neighbors(neighbor.id)
                        if (len(neighbors) <= 2 or  # Current carbon has at most 2 connections (triple bond + 1)
                            len(neighbor_connections) <= 2):  # Partner carbon has at most 2 connections
                            node.is_reactive_elec = True
                            break

    def _find_nh2_pattern(self) -> None:
        for node in self.nodes:
            if node.element == 'N':
                neighbors = self.get_neighbors(node.id)
                h_count = sum(1 for n, _ in neighbors if n.element == 'H')
                heavy_count = sum(1 for n, _ in neighbors if n.element != 'H')

                if h_count == 2 and heavy_count == 1:
                    node.is_reactive_nuc = True

    def _find_nh_pattern(self) -> None:
        for node in self.nodes:
            if node.element == 'N':
                neighbors = self.get_neighbors(node.id)
                h_count = sum(1 for n, _ in neighbors if n.element == 'H')
                has_carbon = any(n.element == 'C' for n, _ in neighbors)

                if h_count == 1 and has_carbon:
                    node.is_reactive_nuc = True

    def _find_cooh_pattern(self) -> None:
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
        """Convert to dictionary format."""
        return {
            'nodes': [n.to_dict() for n in self.nodes],
            'edges': [e.to_dict() for e in self.edges]
        }
