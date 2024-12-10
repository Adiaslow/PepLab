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
        """Find and mark reactive sites with comprehensive debug output."""
        self.logger.warning(f"\nAnalyzing molecule for reactive sites:")
        self.logger.warning(f"Looking for nucleophile: {nuc_pattern}")
        self.logger.warning(f"Looking for electrophile: {elec_pattern}")

        # Print molecule composition
        elements = [node.element for node in self.nodes]
        self.logger.warning(f"Molecule contains {len(self.nodes)} atoms: {elements}")

        # Print all bonds
        for edge in self.edges:
            from_atom = next(n for n in self.nodes if n.id == edge.from_idx)
            to_atom = next(n for n in self.nodes if n.id == edge.to_idx)
            self.logger.warning(f"Bond: {from_atom.element}{edge.from_idx}-{to_atom.element}{edge.to_idx} ({edge.bond_type})")

        # Check for NH2 groups
        self.logger.warning("\nSearching for NH2 groups:")
        for node in self.nodes:
            if node.element == 'N':
                neighbors = self.get_neighbors(node.id)
                h_count = sum(1 for n, _ in neighbors if n.element == 'H')
                non_h_neighbors = [(n, e) for n, e in neighbors if n.element != 'H']

                self.logger.warning(f"Nitrogen {node.id}:")
                self.logger.warning(f"  H count: {h_count}")
                self.logger.warning(f"  Non-H neighbors: {[(n.element, e.bond_type) for n, e in non_h_neighbors]}")

                if h_count == 2 and len(non_h_neighbors) == 1:
                    n_c = non_h_neighbors[0]
                    if n_c[0].element == 'C' and n_c[1].bond_type == 'SINGLE':
                        node.is_reactive_nuc = True
                        self.logger.warning(f"  Found NH2 group at N{node.id}")

        # Check for COOH groups
        self.logger.warning("\nSearching for COOH groups:")
        for node in self.nodes:
            if node.element == 'C':
                neighbors = self.get_neighbors(node.id)
                self.logger.warning(f"Carbon {node.id}:")
                self.logger.warning(f"  Neighbors: {[(n.element, e.bond_type) for n, e in neighbors]}")

                # Count different types of connections
                double_o = False
                single_oh = False
                carbon_single = False

                for neighbor, edge in neighbors:
                    if neighbor.element == 'O':
                        if edge.bond_type == 'DOUBLE':
                            double_o = True
                            self.logger.warning(f"  Found C=O bond")
                        elif edge.bond_type == 'SINGLE':
                            o_neighbors = self.get_neighbors(neighbor.id)
                            if any(n.element == 'H' for n, _ in o_neighbors):
                                single_oh = True
                                self.logger.warning(f"  Found C-OH group")
                    elif neighbor.element == 'C' and edge.bond_type == 'SINGLE':
                        carbon_single = True
                        self.logger.warning(f"  Found C-C bond")

                if double_o and single_oh and carbon_single:
                    node.is_reactive_elec = True
                    self.logger.warning(f"  Found complete COOH group at C{node.id}")

        # Final summary
        nuc_sites = [n for n in self.nodes if n.is_reactive_nuc]
        elec_sites = [n for n in self.nodes if n.is_reactive_elec]
        self.logger.warning("\nSummary:")
        self.logger.warning(f"Found {len(nuc_sites)} nucleophilic sites: {[(n.id, n.element) for n in nuc_sites]}")
        self.logger.warning(f"Found {len(elec_sites)} electrophilic sites: {[(n.id, n.element) for n in elec_sites]}")

    def _find_peptide_reactive_sites(self) -> None:
        """Find NH2/NH and COOH groups with detailed logging."""
        # Find NH2 groups
        for node in self.nodes:
            if node.element == 'N':
                neighbors = self.get_neighbors(node.id)
                non_h_neighbors = [(n, e) for n, e in neighbors if n.element != 'H']
                h_neighbors = [(n, e) for n, e in neighbors if n.element == 'H']

                self.logger.debug(f"Nitrogen atom {node.id}: {len(non_h_neighbors)} non-H neighbors, {len(h_neighbors)} H neighbors")

                # Check for NH2 pattern
                if (len(non_h_neighbors) == 1 and  # One non-H neighbor
                    len(h_neighbors) == 2 and      # Two H neighbors
                    non_h_neighbors[0][0].element == 'C' and  # Connected to carbon
                    non_h_neighbors[0][1].bond_type == 'SINGLE'):  # By single bond
                    node.is_reactive_nuc = True
                    self.logger.debug(f"Found NH2 group at nitrogen {node.id}")
                    continue

                # Check for NH pattern (secondary amine)
                if (len(non_h_neighbors) == 2 and  # Two non-H neighbors
                    len(h_neighbors) == 1 and      # One H neighbor
                    all(n.element == 'C' and e.bond_type == 'SINGLE'
                        for n, e in non_h_neighbors)):  # Connected to carbons by single bonds
                    node.is_reactive_nuc = True
                    self.logger.debug(f"Found NH group at nitrogen {node.id}")

        # Find COOH groups
        for node in self.nodes:
            if node.element == 'C':
                neighbors = self.get_neighbors(node.id)

                # Log all neighbors for debugging
                neighbor_info = [(n.element, e.bond_type) for n, e in neighbors]
                self.logger.debug(f"Carbon atom {node.id} neighbors: {neighbor_info}")

                # Count different types of connections
                double_o = False
                single_oh = False
                carbon_single = False

                for neighbor, edge in neighbors:
                    if neighbor.element == 'O':
                        if edge.bond_type == 'DOUBLE':
                            double_o = True
                            self.logger.debug(f"Found C=O at carbon {node.id}")
                        elif edge.bond_type == 'SINGLE':
                            # Check if this O has a hydrogen
                            o_neighbors = self.get_neighbors(neighbor.id)
                            if any(n.element == 'H' for n, _ in o_neighbors):
                                single_oh = True
                                self.logger.debug(f"Found C-OH at carbon {node.id}")
                    elif (neighbor.element == 'C' and
                            edge.bond_type == 'SINGLE'):
                        carbon_single = True
                        self.logger.debug(f"Found C-C at carbon {node.id}")

                if double_o and single_oh and carbon_single:
                    node.is_reactive_elec = True
                    self.logger.debug(f"Found complete COOH group at carbon {node.id}")

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
