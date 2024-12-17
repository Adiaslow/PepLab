import logging
from typing import List, Tuple, Dict
from uuid import uuid4
from rdkit import Chem

from ..molecule.atom import GraphNode
from ..molecule.bond import GraphEdge
from ..reaction.reactive_pattern import ReactivePattern
from ..reaction.reactive_type import ReactiveType

class MolecularGraph:
    def __init__(self):
        self.nodes: List[GraphNode] = []
        self.edges: List[GraphEdge] = []
        self.logger = logging.getLogger(self.__class__.__name__)
        self._reactive_patterns: List[ReactivePattern] = []

    @staticmethod
    def from_smiles(smiles: str) -> 'MolecularGraph':
        """Create molecular graph from SMILES with UUID generation."""
        graph = MolecularGraph()
        try:
            graph.logger.debug(f"Parsing SMILES: {smiles}")
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError(f"Failed to parse SMILES: {smiles}")

            Chem.SanitizeMol(mol)
            mol = Chem.AddHs(mol)

            # Add atoms with UUIDs
            for atom in mol.GetAtoms():
                idx = atom.GetIdx()
                node = GraphNode(
                    id=str(uuid4()),
                    index=idx,
                    element=atom.GetSymbol(),
                    atomic_num=atom.GetAtomicNum(),
                    formal_charge=atom.GetFormalCharge(),
                    implicit_valence=atom.GetImplicitValence(),
                    explicit_valence=atom.GetExplicitValence(),
                    is_aromatic=atom.GetIsAromatic(),
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

            # Add bonds with UUIDs
            for bond in mol.GetBonds():
                edge = GraphEdge(
                    id=str(uuid4()),
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
            graph.logger.error(f"Error creating molecular graph from SMILES: {str(e)}")
            raise

    def get_neighbors(self, node_index: int) -> List[Tuple[GraphNode, GraphEdge]]:
        """Get all neighbors of a node with connecting edges using node index."""
        return [(n, e) for e in self.edges
                if (e.from_idx == node_index and (n := next(n for n in self.nodes if n.index == e.to_idx)))
                or (e.to_idx == node_index and (n := next(n for n in self.nodes if n.index == e.from_idx)))]

    def find_reactive_sites(self, nuc_pattern: str, elec_pattern: str) -> None:
        """Find and mark reactive sites with comprehensive pattern matching."""
        # Clear existing sites
        self._reactive_patterns.clear()
        for node in self.nodes:
            node.is_reactive_nuc = False
            node.is_reactive_elec = False
            node.is_reactive_click = False

        # Find primary amines (NH2)
        self._find_amine_patterns()

        # Find carboxyl groups (COOH)
        self._find_carboxyl_patterns()

        # Find click chemistry patterns if specified
        if nuc_pattern == 'N3' or elec_pattern == 'N3':
            self._find_azide_patterns()
        if nuc_pattern == 'C#C' or elec_pattern == 'C#C':
            self._find_alkyne_patterns()

        # Log found patterns
        pattern_counts = {}
        for pattern in self._reactive_patterns:
            pattern_counts[pattern.type] = pattern_counts.get(pattern.type, 0) + 1

        if not pattern_counts:
            self.logger.warning("No reactive patterns found")
        else:
            self.logger.debug(f"Found reactive patterns: {pattern_counts}")

    def _find_amine_patterns(self) -> None:
        """Find primary and secondary amine patterns."""
        for node in self.nodes:
            if node.element != 'N':
                continue

            neighbors = self.get_neighbors(node.index)
            h_neighbors = [(n, e) for n, e in neighbors if n.element == 'H']
            heavy_neighbors = [(n, e) for n, e in neighbors if n.element != 'H']

            # Check for NH2 pattern
            if len(h_neighbors) == 2 and len(heavy_neighbors) == 1:
                heavy_atom, heavy_bond = heavy_neighbors[0]
                if heavy_atom.element == 'C' and heavy_bond.bond_type == 'SINGLE':
                    node.is_reactive_nuc = True
                    self._reactive_patterns.append(ReactivePattern(
                        type=ReactiveType.NH2,
                        atoms=[node, *[n for n, _ in h_neighbors], heavy_atom],
                        bonds=[e for _, e in neighbors],
                        pattern_id=str(uuid4())
                    ))

            # Check for NH pattern
            elif len(h_neighbors) == 1 and len(heavy_neighbors) == 2:
                if all(n.element == 'C' and e.bond_type == 'SINGLE'
                      for n, e in heavy_neighbors):
                    node.is_reactive_nuc = True
                    self._reactive_patterns.append(ReactivePattern(
                        type=ReactiveType.NH,
                        atoms=[node, h_neighbors[0][0], *[n for n, _ in heavy_neighbors]],
                        bonds=[e for _, e in neighbors],
                        pattern_id=str(uuid4())
                    ))

    def _find_carboxyl_patterns(self) -> None:
        """Find carboxyl (COOH) patterns."""
        for node in self.nodes:
            if node.element != 'C':
                continue

            neighbors = self.get_neighbors(node.index)
            o_neighbors = [(n, e) for n, e in neighbors if n.element == 'O']
            heavy_neighbors = [(n, e) for n, e in neighbors if n.element not in {'O', 'H'}]

            if len(o_neighbors) != 2 or not heavy_neighbors:
                continue

            # Look for C(=O)OH pattern
            double_o = None
            single_o = None
            for o_atom, o_bond in o_neighbors:
                if o_bond.bond_type == 'DOUBLE':
                    double_o = (o_atom, o_bond)
                elif o_bond.bond_type == 'SINGLE':
                    o_neighbors2 = self.get_neighbors(o_atom.index)
                    h_neighbors = [(n, e) for n, e in o_neighbors2 if n.element == 'H']
                    if len(h_neighbors) == 1:
                        single_o = (o_atom, o_bond)

            if double_o and single_o and len(heavy_neighbors) == 1:
                heavy_atom, heavy_bond = heavy_neighbors[0]
                if heavy_bond.bond_type == 'SINGLE':
                    node.is_reactive_elec = True
                    self._reactive_patterns.append(ReactivePattern(
                        type=ReactiveType.COOH,
                        atoms=[node, double_o[0], single_o[0],
                              next(n for n, _ in self.get_neighbors(single_o[0].index)
                                  if n.element == 'H')],
                        bonds=[double_o[1], single_o[1],
                              next(e for _, e in self.get_neighbors(single_o[0].index)
                                  if next(n for n in self.nodes
                                        if n.index == (e.to_idx if e.from_idx == single_o[0].index
                                                     else e.from_idx)).element == 'H')],
                        pattern_id=str(uuid4())
                    ))

    def _find_azide_patterns(self) -> None:
        """Find azide (N3) patterns."""
        for node in self.nodes:
            if node.element != 'N':
                continue

            neighbors = self.get_neighbors(node.index)
            if len(neighbors) != 2:
                continue

            # Find middle nitrogen
            n2_pair = next(((n, e) for n, e in neighbors
                           if n.element == 'N' and e.bond_type == 'SINGLE'), None)
            if not n2_pair:
                continue

            n2, n1_n2_bond = n2_pair
            n2_neighbors = self.get_neighbors(n2.index)
            if len(n2_neighbors) != 2:
                continue

            # Find terminal nitrogen
            n3_pair = next(((n, e) for n, e in n2_neighbors
                           if n.element == 'N' and n.index != node.index
                           and e.bond_type == 'TRIPLE'), None)
            if not n3_pair:
                continue

            n3, n2_n3_bond = n3_pair
            if len(self.get_neighbors(n3.index)) == 1:
                node.is_reactive_click = True
                self._reactive_patterns.append(ReactivePattern(
                    type=ReactiveType.AZIDE,
                    atoms=[node, n2, n3],
                    bonds=[n1_n2_bond, n2_n3_bond],
                    pattern_id=str(uuid4())
                ))

    def _find_alkyne_patterns(self) -> None:
        """Find terminal alkyne (C≡C) patterns."""
        for node in self.nodes:
            if node.element != 'C':
                continue

            neighbors = self.get_neighbors(node.index)
            triple_bond_pair = next(((n, e) for n, e in neighbors
                                   if n.element == 'C' and e.bond_type == 'TRIPLE'), None)
            if not triple_bond_pair:
                continue

            other_c, triple_bond = triple_bond_pair

            # Get non-triple bond connections
            current_other_bonds = [(n, e) for n, e in neighbors if e != triple_bond]
            other_c_neighbors = self.get_neighbors(other_c.index)
            other_c_other_bonds = [(n, e) for n, e in other_c_neighbors if e != triple_bond]

            # Check for terminal alkyne
            is_current_terminal = len(current_other_bonds) <= 1
            is_other_terminal = len(other_c_other_bonds) <= 1

            if ((is_current_terminal and len(other_c_other_bonds) <= 2) or
                (is_other_terminal and len(current_other_bonds) <= 2)):
                node.is_reactive_click = True
                self._reactive_patterns.append(ReactivePattern(
                    type=ReactiveType.ALKYNE,
                    atoms=[node, other_c],
                    bonds=[triple_bond],
                    pattern_id=str(uuid4())
                ))

    def to_dict(self) -> Dict:
        """Convert to dictionary format."""
        return {
            'nodes': [n.to_dict() for n in self.nodes],
            'edges': [e.to_dict() for e in self.edges]
        }
