# @title Peptide Builder

"""Peptide builder with threading support and progress tracking."""

import logging
from typing import List, Dict, Optional
from dataclasses import dataclass
import itertools
from itertools import product
import copy
import concurrent.futures
from tqdm import tqdm
import threading
from queue import Queue

from ..graph.molecule_graph import MolecularGraph
from ..molecule.atom import GraphNode
from ..molecule.bond import GraphEdge
from ..molecule.peptide import PeptideInfo
from .library import LibraryInfo

class PeptideBuilder:
    """Handles construction of peptides with threading support and click chemistry."""

    def __init__(self, max_workers: Optional[int] = None):
        """Initialize builder with logging and threading settings."""
        self.logger = logging.getLogger(self.__class__.__name__)
        self.max_workers = max_workers
        self._progress_queue = Queue()
        self._total_combinations = 0

    def _process_combination(self, combo: tuple) -> Optional[PeptideInfo]:
        """Process a single combination of residues."""
        try:
            residue_graphs = [res['graph'] for res in combo]
            peptide = self.build_linear_peptide(residue_graphs)

            # Check if click chemistry is possible
            if self._cyclize and self._has_click_pairs(combo):
                peptide = self.click_cyclize_peptide(peptide)
                peptide_type = 'click-cyclic'
            elif self._cyclize:
                peptide = self.cyclize_peptide(peptide)
                peptide_type = 'cyclic'
            else:
                peptide_type = 'linear'

            self._progress_queue.put(1)  # Signal progress

            return PeptideInfo(
                graph=peptide.to_dict(),
                sequence=[res['id'] for res in combo],
                peptide_type=peptide_type,
                smiles=None
            )
        except Exception as e:
            self.logger.warning(f"Error generating peptide: {str(e)}")
            self._progress_queue.put(1)  # Signal progress even for failures
            return None

    def _progress_tracker(self, total: int):
        """Track progress using tqdm."""
        with tqdm(total=total, desc="Generating peptides") as pbar:
            completed = 0
            while completed < total:
                self._progress_queue.get()
                completed += 1
                pbar.update(1)

    def enumerate_library(self, library_info: LibraryInfo, cyclize: bool = False, click_cyclize: bool = False) -> List[PeptideInfo]:
        """Enumerate peptide library with threading support."""
        peptides = []
        self._cyclize = cyclize
        self._click_cyclize = click_cyclize  # Store click_cyclize parameter


        try:
            positions = sorted(library_info.positions.keys())
            residue_options = []

            # Process each position
            for pos in positions:
                pos_residues = []
                for res_key, res_info in library_info.positions[pos].residues.items():
                    try:
                        graph = MolecularGraph()
                        graph = graph.from_smiles(res_info.smiles)
                        graph.find_reactive_sites(
                            res_info.nucleophile,
                            res_info.electrophile
                        )

                        pos_residues.append({
                            'graph': graph,
                            'id': res_info.name,
                            'nuc': res_info.nucleophile,
                            'elec': res_info.electrophile
                        })
                    except Exception as e:
                        self.logger.error(
                            f"Error processing residue {res_info.name}: {str(e)}"
                        )

                if not pos_residues:
                    raise ValueError(f"No valid residues for position {pos}")

                residue_options.append(pos_residues)

            # Calculate total combinations
            self._total_combinations = 1
            for options in residue_options:
                self._total_combinations *= len(options)

            # Start progress tracking thread
            progress_thread = threading.Thread(
                target=self._progress_tracker,
                args=(self._total_combinations,)
            )
            progress_thread.start()

            # Try multithreading first
            try:
                with concurrent.futures.ThreadPoolExecutor(max_workers=self.max_workers) as executor:
                    futures = []
                    for combo in product(*residue_options):
                        futures.append(executor.submit(self._process_combination, combo))

                    for future in concurrent.futures.as_completed(futures):
                        result = future.result()
                        if result is not None:
                            peptides.append(result)

            except (ImportError, RuntimeError) as e:
                self.logger.warning(f"Threading unavailable: {str(e)}. Falling back to single-threaded mode.")
                # Clear progress queue and restart progress tracking for single-threaded mode
                while not self._progress_queue.empty():
                    self._progress_queue.get()

                # Process combinations sequentially
                for combo in product(*residue_options):
                    result = self._process_combination(combo)
                    if result is not None:
                        peptides.append(result)

            # Wait for progress tracking to complete
            progress_thread.join()

            self.logger.info(
                f"Generated {len(peptides)} out of {self._total_combinations} possible peptides"
            )
            return peptides

        except Exception as e:
            self.logger.error(f"Error enumerating library: {str(e)}")
            return []

    def build_linear_peptide(self, residue_graphs: List[MolecularGraph]) -> MolecularGraph:
        """Build linear peptide from residue graphs."""
        if not residue_graphs:
            raise ValueError("No residue graphs provided")

        peptide = copy.deepcopy(residue_graphs[0])

        for i in range(1, len(residue_graphs)):
            peptide = self._form_peptide_bond(peptide, residue_graphs[i])

        return peptide

    def _form_peptide_bond(self, res1: MolecularGraph, res2: MolecularGraph) -> MolecularGraph:
        """Form peptide bond between residues with detailed debug output."""
        self.logger.warning("\nAttempting to form peptide bond:")

        # Log first residue details
        self.logger.warning("First residue details:")
        for node in res1.nodes:
            if node.element == 'N':
                neighbors = res1.get_neighbors(node.id)
                self.logger.warning(f"N{node.id}:")
                self.logger.warning(f"  Is reactive: {node.is_reactive_nuc}")
                self.logger.warning(f"  Neighbors: {[(n.element, e.bond_type) for n, e in neighbors]}")

        # Log second residue details
        self.logger.warning("\nSecond residue details:")
        for node in res2.nodes:
            if node.element == 'C':
                neighbors = res2.get_neighbors(node.id)
                self.logger.warning(f"C{node.id}:")
                self.logger.warning(f"  Is reactive: {node.is_reactive_elec}")
                self.logger.warning(f"  Neighbors: {[(n.element, e.bond_type) for n, e in neighbors]}")

        # Find reactive sites
        nuc_site = next(
            (n for n in res1.nodes if n.is_reactive_nuc and not n.is_reactive_click),
            None
        )
        elec_site = next(
            (n for n in res2.nodes if n.is_reactive_elec and not n.is_reactive_click),
            None
        )

        # Log found sites
        self.logger.warning("\nReactive sites found:")
        self.logger.warning(f"Nucleophilic site: {nuc_site.id if nuc_site else None}")
        self.logger.warning(f"Electrophilic site: {elec_site.id if elec_site else None}")

        if not (nuc_site and elec_site):
            raise ValueError("Could not find required reactive sites for peptide bond formation")

        res1_mod = copy.deepcopy(res1)
        res2_mod = copy.deepcopy(res2)

        res1_mod = self._remove_h_from_nh(res1_mod, nuc_site.id)
        res2_mod = self._remove_oh_from_cooh(res2_mod, elec_site.id)

        max_id = max(n.id for n in res1_mod.nodes) + 1
        for node in res2_mod.nodes:
            node.id += max_id
        for edge in res2_mod.edges:
            edge.from_idx += max_id
            edge.to_idx += max_id

        combined = MolecularGraph()
        combined.nodes = res1_mod.nodes + res2_mod.nodes
        combined.edges = res1_mod.edges + res2_mod.edges

        combined.edges.append(GraphEdge(
            from_idx=nuc_site.id,
            to_idx=elec_site.id + max_id,
            bond_type='SINGLE',
            is_aromatic=False,
            is_conjugated=False,
            in_ring=False,
            stereo='NONE'
        ))

        for node in combined.nodes:
            if node.id == nuc_site.id:
                node.is_reactive_nuc = False
            if node.id == elec_site.id + max_id:
                node.is_reactive_elec = False

        return self._reindex_graph(combined)

    def _remove_h_from_nh(self, graph: MolecularGraph, node_id: int) -> MolecularGraph:
        """Remove hydrogen from NH/NH2 group."""
        mod_graph = copy.deepcopy(graph)

        for edge in mod_graph.edges[:]:  # Use slice to allow modification during iteration
            if edge.from_idx == node_id:
                h_node = next((n for n in mod_graph.nodes if n.id == edge.to_idx and n.element == 'H'), None)
                if h_node:
                    mod_graph.nodes.remove(h_node)
                    mod_graph.edges.remove(edge)
                    break
            elif edge.to_idx == node_id:
                h_node = next((n for n in mod_graph.nodes if n.id == edge.from_idx and n.element == 'H'), None)
                if h_node:
                    mod_graph.nodes.remove(h_node)
                    mod_graph.edges.remove(edge)
                    break

        return mod_graph

    def _remove_oh_from_cooh(self, graph: MolecularGraph, node_id: int) -> MolecularGraph:
        """Remove OH group from COOH."""
        mod_graph = copy.deepcopy(graph)

        for edge in mod_graph.edges[:]:
            if edge.bond_type == 'SINGLE':
                if edge.from_idx == node_id:
                    o_node = next((n for n in mod_graph.nodes if n.id == edge.to_idx and n.element == 'O'), None)
                elif edge.to_idx == node_id:
                    o_node = next((n for n in mod_graph.nodes if n.id == edge.from_idx and n.element == 'O'), None)
                else:
                    continue

                if o_node:
                    # Find H attached to O
                    for h_edge in mod_graph.edges[:]:
                        if h_edge.from_idx == o_node.id:
                            h_node = next((n for n in mod_graph.nodes if n.id == h_edge.to_idx and n.element == 'H'), None)
                        elif h_edge.to_idx == o_node.id:
                            h_node = next((n for n in mod_graph.nodes if n.id == h_edge.from_idx and n.element == 'H'), None)
                        else:
                            continue

                        if h_node:
                            mod_graph.nodes.remove(o_node)
                            mod_graph.nodes.remove(h_node)
                            mod_graph.edges.remove(edge)
                            mod_graph.edges.remove(h_edge)
                            return mod_graph

        return mod_graph

    def _has_click_pairs(self, combo: tuple) -> bool:
        """Check if the combination has valid azide and alkyne pairs for click chemistry."""
        # Find residues with potential azide/alkyne groups
        azide_residues = [res for res in combo
                            if res.get('nuc') == 'N3' or res.get('elec') == 'N3']
        alkyne_residues = [res for res in combo
                            if res.get('nuc') == 'C#C' or res.get('elec') == 'C#C']

        if not (azide_residues and alkyne_residues):
            return False

        # For each potential pair, validate the reactive sites
        for azide_res in azide_residues:
            graph = azide_res['graph']
            azide_n = next((n for n in graph.nodes if n.is_reactive_click), None)
            if not azide_n:
                continue

            for alkyne_res in alkyne_residues:
                if azide_res == alkyne_res:  # Skip same residue
                    continue

                graph = alkyne_res['graph']
                alkyne_c = next((n for n in graph.nodes if n.is_reactive_click), None)
                if not alkyne_c:
                    continue

                # If we find a valid pair, return True
                if self.validate_click_chemistry_pair(azide_n, alkyne_c, graph):
                    return True

        return False

    def click_cyclize_peptide(self, linear_peptide: MolecularGraph) -> MolecularGraph:
        """Form triazole ring through click chemistry cyclization."""
        cyclic = copy.deepcopy(linear_peptide)

        # Find azide and alkyne sites (using click reactive sites)
        azide_site = next((n for n in cyclic.nodes if n.is_reactive_click
                            and n.element == 'N'), None)
        alkyne_site = next((n for n in cyclic.nodes if n.is_reactive_click
                            and n.element == 'C'), None)

        if not azide_site or not alkyne_site:
            raise ValueError("Cannot find required click chemistry reactive sites")

        # Find the complete azide chain
        azide_chain = []
        current = azide_site
        while len(azide_chain) < 3:
            azide_chain.append(current)
            for edge in cyclic.edges:
                if edge.from_idx == current.id:
                    next_node = next((n for n in cyclic.nodes if n.id == edge.to_idx
                                    and n.element == 'N' and n not in azide_chain), None)
                    if next_node:
                        current = next_node
                        break
                elif edge.to_idx == current.id:
                    next_node = next((n for n in cyclic.nodes if n.id == edge.from_idx
                                    and n.element == 'N' and n not in azide_chain), None)
                    if next_node:
                        current = next_node
                        break
            else:
                break

        # Find the alkyne carbons
        alkyne_chain = []
        current = alkyne_site
        for edge in cyclic.edges:
            if edge.bond_type == 'TRIPLE':
                if edge.from_idx == current.id:
                    next_node = next((n for n in cyclic.nodes if n.id == edge.to_idx
                                    and n.element == 'C'), None)
                    if next_node:
                        alkyne_chain = [current, next_node]
                        break
                elif edge.to_idx == current.id:
                    next_node = next((n for n in cyclic.nodes if n.id == edge.from_idx
                                    and n.element == 'C'), None)
                    if next_node:
                        alkyne_chain = [current, next_node]
                        break

        if len(azide_chain) != 3 or len(alkyne_chain) != 2:
            raise ValueError("Could not find complete azide or alkyne groups")

        # Remove all existing bonds between atoms in the chains
        edges_to_remove = []
        for edge in cyclic.edges:
            from_in_azide = any(n.id == edge.from_idx for n in azide_chain)
            to_in_azide = any(n.id == edge.to_idx for n in azide_chain)
            from_in_alkyne = any(n.id == edge.from_idx for n in alkyne_chain)
            to_in_alkyne = any(n.id == edge.to_idx for n in alkyne_chain)

            if (from_in_azide and to_in_azide) or (from_in_alkyne and to_in_alkyne):
                edges_to_remove.append(edge)

        for edge in edges_to_remove:
            if edge in cyclic.edges:
                cyclic.edges.remove(edge)

        # Remove hydrogens from all atoms involved
        for node in azide_chain + alkyne_chain:
            edges_to_remove = []
            nodes_to_remove = []
            for edge in cyclic.edges:
                if edge.from_idx == node.id or edge.to_idx == node.id:
                    other_id = edge.to_idx if edge.from_idx == node.id else edge.from_idx
                    other_node = next((n for n in cyclic.nodes if n.id == other_id), None)
                    if other_node and other_node.element == 'H':
                        edges_to_remove.append(edge)
                        nodes_to_remove.append(other_node)

            for edge in edges_to_remove:
                if edge in cyclic.edges:
                    cyclic.edges.remove(edge)
            for node_to_remove in nodes_to_remove:
                if node_to_remove in cyclic.nodes:
                    cyclic.nodes.remove(node_to_remove)

        # Create the triazole ring with correct alternating single/double bonds
        # Pattern: N1-N2=N3-C4=C5-N1
        new_edges = [
            GraphEdge(  # N1-N2
                from_idx=azide_chain[0].id,
                to_idx=azide_chain[1].id,
                bond_type='SINGLE',
                is_aromatic=False,
                is_conjugated=True,
                in_ring=True,
                stereo='NONE'
            ),
            GraphEdge(  # N2=N3
                from_idx=azide_chain[1].id,
                to_idx=azide_chain[2].id,
                bond_type='DOUBLE',
                is_aromatic=False,
                is_conjugated=True,
                in_ring=True,
                stereo='NONE'
            ),
            GraphEdge(  # N3-C4
                from_idx=azide_chain[2].id,
                to_idx=alkyne_chain[0].id,
                bond_type='SINGLE',
                is_aromatic=False,
                is_conjugated=True,
                in_ring=True,
                stereo='NONE'
            ),
            GraphEdge(  # C4=C5
                from_idx=alkyne_chain[0].id,
                to_idx=alkyne_chain[1].id,
                bond_type='DOUBLE',
                is_aromatic=False,
                is_conjugated=True,
                in_ring=True,
                stereo='NONE'
            ),
            GraphEdge(  # C5-N1
                from_idx=alkyne_chain[1].id,
                to_idx=azide_chain[0].id,
                bond_type='SINGLE',
                is_aromatic=False,
                is_conjugated=True,
                in_ring=True,
                stereo='NONE'
            )
        ]

        # Add new bonds after verifying they don't already exist
        for new_edge in new_edges:
            if not any(e.from_idx == new_edge.from_idx and e.to_idx == new_edge.to_idx for e in cyclic.edges):
                cyclic.edges.append(new_edge)

        # Update atomic properties
        for node in cyclic.nodes:
            if node in azide_chain or node in alkyne_chain:
                node.is_reactive_nuc = False
                node.is_reactive_elec = False
                node.implicit_h_count = 0
                node.explicit_h_count = 0
                node.formal_charge = 0
                node.is_aromatic = False
                node.is_conjugated = True  # Set conjugated to true for all ring atoms
                if node.element == 'N':
                    node.valence = 3
                elif node.element == 'C':
                    node.valence = 4

        # Reindex and clean up the graph
        return self._reindex_graph(cyclic)

    def _get_azide_chain(self, graph: MolecularGraph, start_id: int) -> List[GraphNode]:
        """Get the complete azide chain starting from the first nitrogen."""
        chain = [next(n for n in graph.nodes if n.id == start_id)]
        current_id = start_id
        visited = {current_id}

        while True:
            next_n = None
            for edge in graph.edges:
                if edge.from_idx == current_id or edge.to_idx == current_id:
                    other_id = edge.to_idx if edge.from_idx == current_id else edge.from_idx
                    node = next((n for n in graph.nodes if n.id == other_id), None)
                    if node and node.element == 'N' and node.id not in visited:
                        next_n = node
                        chain.append(node)
                        visited.add(node.id)
                        current_id = node.id
                        break
            if not next_n:
                break

        return chain

    def _get_alkyne_chain(self, graph: MolecularGraph, start_id: int) -> List[GraphNode]:
        """Get the complete alkyne chain starting from the first carbon."""
        chain = [next(n for n in graph.nodes if n.id == start_id)]

        # Find the carbon connected by triple bond
        for edge in graph.edges:
            if edge.bond_type == 'TRIPLE' and (edge.from_idx == start_id or edge.to_idx == start_id):
                other_id = edge.to_idx if edge.from_idx == start_id else edge.from_idx
                other_node = next((n for n in graph.nodes if n.id == other_id), None)
                if other_node and other_node.element == 'C':
                    chain.append(other_node)
                    break

        return chain

    def _is_azide_nitrogen(self, node: GraphNode, graph: MolecularGraph) -> bool:
        """
        Check if node is the first nitrogen of an azide group.
        Verifies N-N≡N pattern with correct bond types and connectivity.
        """
        if node.element != 'N':
            return False

        # Get connected atoms and their bonds
        connected = []
        for edge in graph.edges:
            if edge.from_idx == node.id:
                atom = next(n for n in graph.nodes if n.id == edge.to_idx)
                connected.append((atom, edge))
            elif edge.to_idx == node.id:
                atom = next(n for n in graph.nodes if n.id == edge.from_idx)
                connected.append((atom, edge))

        # First N should have exactly two connections
        if len(connected) != 2:
            return False

        # Find the second nitrogen (N2)
        n2_pair = next(((atom, edge) for atom, edge in connected
                        if atom.element == 'N' and edge.bond_type == 'SINGLE'), None)
        if not n2_pair:
            return False

        n2, n1_n2_bond = n2_pair

        # Check N2 connections
        n2_connected = []
        for edge in graph.edges:
            if edge.from_idx == n2.id:
                atom = next(n for n in graph.nodes if n.id == edge.to_idx)
                n2_connected.append((atom, edge))
            elif edge.to_idx == n2.id:
                atom = next(n for n in graph.nodes if n.id == edge.from_idx)
                n2_connected.append((atom, edge))

        # N2 should have exactly two connections (N1 and N3)
        if len(n2_connected) != 2:
            return False

        # Find N3 (should be connected to N2 with a triple bond)
        n3_pair = next(((atom, edge) for atom, edge in n2_connected
                        if atom.element == 'N' and atom.id != node.id
                        and edge.bond_type == 'TRIPLE'), None)
        if not n3_pair:
            return False

        n3, n2_n3_bond = n3_pair

        # Check N3 has only one connection (to N2)
        n3_connections = sum(1 for edge in graph.edges
                            if edge.from_idx == n3.id or edge.to_idx == n3.id)
        if n3_connections != 1:
            return False

        return True

    def _is_alkyne_carbon(self, node: GraphNode, graph: MolecularGraph) -> bool:
        """
        Check if node is part of a terminal alkyne group.
        Verifies C≡C pattern with correct bond types and connectivity.
        """
        if node.element != 'C':
            return False

        # Get all connections and their bond types
        connected = []
        for edge in graph.edges:
            if edge.from_idx == node.id:
                atom = next(n for n in graph.nodes if n.id == edge.to_idx)
                connected.append((atom, edge))
            elif edge.to_idx == node.id:
                atom = next(n for n in graph.nodes if n.id == edge.from_idx)
                connected.append((atom, edge))

        # Find triple bond and connected carbon
        triple_bond_pair = next(((atom, edge) for atom, edge in connected
                                if atom.element == 'C' and edge.bond_type == 'TRIPLE'), None)
        if not triple_bond_pair:
            return False

        other_c, triple_bond = triple_bond_pair

        # Check connections of both carbons
        node_other_connections = [(atom, edge) for atom, edge in connected if edge != triple_bond]
        other_c_connections = []
        for edge in graph.edges:
            if edge != triple_bond:
                if edge.from_idx == other_c.id:
                    atom = next(n for n in graph.nodes if n.id == edge.to_idx)
                    other_c_connections.append((atom, edge))
                elif edge.to_idx == other_c.id:
                    atom = next(n for n in graph.nodes if n.id == edge.from_idx)
                    other_c_connections.append((atom, edge))

        # At least one carbon must be terminal (only one other connection)
        # and the other should have at most two other connections
        is_node_terminal = len(node_other_connections) <= 1
        is_other_terminal = len(other_c_connections) <= 1

        if not (is_node_terminal or is_other_terminal):
            return False

        # The non-terminal carbon should have at most 2 other connections
        if is_node_terminal and len(other_c_connections) > 2:
            return False
        if is_other_terminal and len(node_other_connections) > 2:
            return False
        return True

    def validate_click_chemistry_pair(self, azide_n: GraphNode, alkyne_c: GraphNode,
                                        graph: MolecularGraph) -> bool:
            """Validate that the given pair can participate in click chemistry."""
            if not (azide_n.is_reactive_click and alkyne_c.is_reactive_click):
                return False

            return True

    def _get_connected_atoms(self, node_id: int, graph: MolecularGraph) -> List[GraphNode]:
        """Get all atoms connected to the given node."""
        connected = []
        for edge in graph.edges:
            if edge.from_idx == node_id:
                connected.append(next(n for n in graph.nodes if n.id == edge.to_idx))
            elif edge.to_idx == node_id:
                connected.append(next(n for n in graph.nodes if n.id == edge.from_idx))
        return connected

    def _get_azide_nitrogens(self, graph: MolecularGraph, start_nitrogen_id: int) -> List[GraphNode]:
        """Get all nitrogen atoms in an azide group starting from the first nitrogen."""
        nitrogens = [next(n for n in graph.nodes if n.id == start_nitrogen_id)]
        current = nitrogens[0]

        while True:
            connected = self._get_connected_atoms(current.id, graph)
            next_n = next((n for n in connected if n.element == 'N' and n not in nitrogens), None)
            if not next_n:
                break
            nitrogens.append(next_n)
            current = next_n

        return nitrogens

    def cyclize_peptide(self, linear_peptide: MolecularGraph) -> MolecularGraph:
        """Cyclize linear peptide using peptide bond formation."""
        # Find peptide bond reactive sites (not click sites)
        nuc_site = next(
            (n for n in linear_peptide.nodes
                if n.is_reactive_nuc and not n.is_reactive_click),
            None
        )
        elec_site = next(
            (n for n in linear_peptide.nodes
                if n.is_reactive_elec and not n.is_reactive_click),
            None
        )

        if not (nuc_site and elec_site):
            raise ValueError("Cannot find terminal reactive sites for cyclization")

        cyclic = copy.deepcopy(linear_peptide)
        cyclic = self._remove_h_from_nh(cyclic, nuc_site.id)
        cyclic = self._remove_oh_from_cooh(cyclic, elec_site.id)

        cyclic.edges.append(GraphEdge(
            from_idx=nuc_site.id,
            to_idx=elec_site.id,
            bond_type='SINGLE',
            is_aromatic=False,
            is_conjugated=False,
            in_ring=True,
            stereo='NONE'
        ))

        for node in cyclic.nodes:
            if node.id == nuc_site.id:
                node.is_reactive_nuc = False
            if node.id == elec_site.id:
                node.is_reactive_elec = False

        return self._reindex_graph(cyclic)

    def _reindex_graph(self, graph: MolecularGraph) -> MolecularGraph:
        """Reindex graph nodes and edges sequentially."""
        new_graph = MolecularGraph()
        old_to_new = {}

        # Reindex nodes
        for i, node in enumerate(sorted(graph.nodes, key=lambda x: x.id)):
            new_node = copy.deepcopy(node)
            new_node.id = i
            old_to_new[node.id] = i
            new_graph.nodes.append(new_node)

        # Update edges
        for edge in graph.edges:
            new_edge = copy.deepcopy(edge)
            new_edge.from_idx = old_to_new[edge.from_idx]
            new_edge.to_idx = old_to_new[edge.to_idx]
            new_graph.edges.append(new_edge)

        return new_graph
