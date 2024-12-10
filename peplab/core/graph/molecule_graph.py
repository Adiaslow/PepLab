import logging
from typing import List, Tuple, Optional, Dict
from rdkit import Chem

from ..molecule.atom import GraphNode
from ..molecule.bond import GraphEdge

logging.basicConfig(level=logging.DEBUG)

class MolecularGraph:
    def __init__(self):
        self.nodes: List[GraphNode] = []
        self.edges: List[GraphEdge] = []
        self.logger = logging.getLogger(self.__class__.__name__)

    @staticmethod
    def from_smiles(smiles: str) -> 'MolecularGraph':
        """Create molecular graph from SMILES with debug logging."""
        graph = MolecularGraph()
        try:
            graph.logger.warning(f"\nParsing SMILES: {smiles}")
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError(f"Failed to parse SMILES: {smiles}")

            Chem.SanitizeMol(mol)
            mol = Chem.AddHs(mol)

            # Log molecule composition
            graph.logger.warning(f"Molecule composition:")
            for atom in mol.GetAtoms():
                graph.logger.warning(f"Atom {atom.GetIdx()}: {atom.GetSymbol()}")
                graph.logger.warning(f"  Charge: {atom.GetFormalCharge()}")
                graph.logger.warning(f"  Implicit Hs: {atom.GetNumImplicitHs()}")
                graph.logger.warning(f"  Explicit Hs: {atom.GetNumExplicitHs()}")

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

            # Add bonds and log them
            graph.logger.warning("\nBond information:")
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
                graph.logger.warning(f"Bond {edge.from_idx}-{edge.to_idx}: {edge.bond_type}")

            return graph
        except Exception as e:
            graph.logger.error(f"Error creating molecular graph from SMILES: {str(e)}")
            raise

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
        """Find and mark reactive sites with focused logging."""
        # Clear existing sites
        for node in self.nodes:
            node.is_reactive_nuc = False
            node.is_reactive_elec = False
            node.is_reactive_click = False

        # Search for NH2/COOH
        self._find_peptide_reactive_sites()

        # Check click chemistry patterns
        if nuc_pattern == 'N3' or elec_pattern == 'N3':
            self._find_azide_pattern()
        if nuc_pattern == 'C#C' or elec_pattern == 'C#C':
            self._find_alkyne_pattern()

    def _find_peptide_reactive_sites(self) -> None:
        """Find NH2/NH and COOH groups with minimal logging."""
        found_nh2 = False
        found_cooh = False

        # Find NH2 groups
        for node in self.nodes:
            if node.element == 'N':
                neighbors = self.get_neighbors(node.id)
                h_count = sum(1 for n, _ in neighbors if n.element == 'H')
                non_h_neighbors = [(n, e) for n, e in neighbors if n.element != 'H']

                # Check for NH2 pattern
                if (len(non_h_neighbors) == 1 and h_count == 2 and
                    non_h_neighbors[0][0].element == 'C' and
                    non_h_neighbors[0][1].bond_type == 'SINGLE'):
                    node.is_reactive_nuc = True
                    found_nh2 = True

        # Find COOH groups
        for node in self.nodes:
            if node.element == 'C':
                neighbors = self.get_neighbors(node.id)
                double_o = False
                single_oh = False
                carbon_single = False

                for neighbor, edge in neighbors:
                    if neighbor.element == 'O':
                        if edge.bond_type == 'DOUBLE':
                            double_o = True
                        elif edge.bond_type == 'SINGLE':
                            o_neighbors = self.get_neighbors(neighbor.id)
                            if any(n.element == 'H' for n, _ in o_neighbors):
                                single_oh = True
                    elif neighbor.element == 'C' and edge.bond_type == 'SINGLE':
                        carbon_single = True

                if double_o and single_oh and carbon_single:
                    node.is_reactive_elec = True
                    found_cooh = True

        if not (found_nh2 and found_cooh):
            self.logger.warning(f"Missing reactive sites - NH2: {found_nh2}, COOH: {found_cooh}")

    def _find_azide_pattern(self) -> None:
        """Find azide (N3) group."""
        # Implementation stays the same as before, but marks additional reactive site
        for node in self.nodes:
            if node.element == 'N':
                neighbors = self.get_neighbors(node.id)

                # First N should have exactly two connections
                if len(neighbors) != 2:
                    continue

                # Find N2 connected by single bond
                n2_pair = next(((n, e) for n, e in neighbors
                              if n.element == 'N' and e.bond_type == 'SINGLE'), None)
                if not n2_pair:
                    continue

                n2, n1_n2_bond = n2_pair
                n2_neighbors = self.get_neighbors(n2.id)

                # N2 should have exactly two connections
                if len(n2_neighbors) != 2:
                    continue

                # Find N3 connected to N2 by triple bond
                n3_pair = next(((n, e) for n, e in n2_neighbors
                              if n.element == 'N' and n.id != node.id
                              and e.bond_type == 'TRIPLE'), None)
                if not n3_pair:
                    continue

                n3, n2_n3_bond = n3_pair

                # Check N3 has exactly one connection
                n3_neighbors = self.get_neighbors(n3.id)
                if len(n3_neighbors) == 1:
                    # Add azide as additional reactive site
                    node.is_reactive_click = True  # Need to add this field to GraphNode

    def _find_alkyne_pattern(self) -> None:
        """Find terminal alkyne (C≡C) group."""
        # Implementation stays the same as before, but marks additional reactive site
        for node in self.nodes:
            if node.element == 'C':
                neighbors = self.get_neighbors(node.id)

                # Find triple-bonded carbon
                triple_bond_pair = next(((n, e) for n, e in neighbors
                                       if n.element == 'C' and e.bond_type == 'TRIPLE'), None)
                if not triple_bond_pair:
                    continue

                other_c, triple_bond = triple_bond_pair

                # Get non-triple bond connections for both carbons
                current_other_bonds = [(n, e) for n, e in neighbors if e != triple_bond]
                other_c_neighbors = self.get_neighbors(other_c.id)
                other_c_other_bonds = [(n, e) for n, e in other_c_neighbors if e != triple_bond]

                # Check if either carbon is terminal
                is_current_terminal = len(current_other_bonds) <= 1
                is_other_terminal = len(other_c_other_bonds) <= 1

                if ((is_current_terminal and len(other_c_other_bonds) <= 2) or
                    (is_other_terminal and len(current_other_bonds) <= 2)):
                    # Add alkyne as additional reactive site
                    node.is_reactive_click = True  # Need to add this field to GraphNode

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

    def print_node_info(self) -> None:
        """Print detailed information about all nodes for debugging."""
        self.logger.debug("\nNode Information:")
        for node in self.nodes:
            neighbors = self.get_neighbors(node.id)
            neighbor_info = [(n.element, e.bond_type) for n, e in neighbors]
            self.logger.debug(f"Node {node.id}:")
            self.logger.debug(f"  Element: {node.element}")
            self.logger.debug(f"  Neighbors: {neighbor_info}")
            self.logger.debug(f"  Is reactive nuc: {node.is_reactive_nuc}")
            self.logger.debug(f"  Is reactive elec: {node.is_reactive_elec}")
            self.logger.debug(f"  Is reactive click: {node.is_reactive_click}")
