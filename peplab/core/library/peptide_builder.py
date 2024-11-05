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
        """Form peptide bond between residues."""
        nuc_site = next(
            (n for n in res1.nodes if n.is_reactive_nuc),
            None
        )
        elec_site = next(
            (n for n in res2.nodes if n.is_reactive_elec),
            None
        )

        if not (nuc_site and elec_site):
            raise ValueError("Could not find required reactive sites")

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
        """Check if the combination has matching azide and alkyne pairs."""
        has_azide = any(res.get('nuc') == 'N3' or res.get('elec') == 'N3' for res in combo)
        has_alkyne = any(res.get('nuc') == 'C#C' or res.get('elec') == 'C#C' for res in combo)
        return has_azide and has_alkyne

    def click_cyclize_peptide(self, linear_peptide: MolecularGraph) -> MolecularGraph:
        """Form triazole ring by directly replacing reactive sites and cleaning up properly."""
        # Find reactive sites
        azide_site = next((n for n in linear_peptide.nodes if n.is_reactive_nuc), None)
        alkyne_site = next((n for n in linear_peptide.nodes if n.is_reactive_elec), None)

        if not azide_site or not alkyne_site:
            raise ValueError("Cannot find required reactive sites")

        cyclic = copy.deepcopy(linear_peptide)

        # Step 1: Find all nodes to remove from reactive groups
        nodes_to_remove = []  # Changed from set to list
        edges_to_remove = []  # Using a list for edges

        # Collect all nodes and edges connected to reactive sites
        def find_connected_nodes(start_id: int, visited_ids: list):
            for edge in cyclic.edges:
                if edge.from_idx == start_id:
                    # Check if edge is already in edges_to_remove by comparing indices
                    if not any(e.from_idx == edge.from_idx and e.to_idx == edge.to_idx for e in edges_to_remove):
                        edges_to_remove.append(edge)
                    other_id = edge.to_idx
                    if other_id not in visited_ids:
                        visited_ids.append(other_id)
                        if not any(n.id == other_id for n in nodes_to_remove):
                            nodes_to_remove.append(next(n for n in cyclic.nodes if n.id == other_id))
                        find_connected_nodes(other_id, visited_ids)
                elif edge.to_idx == start_id:
                    # Check if edge is already in edges_to_remove by comparing indices
                    if not any(e.from_idx == edge.from_idx and e.to_idx == edge.to_idx for e in edges_to_remove):
                        edges_to_remove.append(edge)
                    other_id = edge.from_idx
                    if other_id not in visited_ids:
                        visited_ids.append(other_id)
                        if not any(n.id == other_id for n in nodes_to_remove):
                            nodes_to_remove.append(next(n for n in cyclic.nodes if n.id == other_id))
                        find_connected_nodes(other_id, visited_ids)

        # Find all reactive group atoms
        find_connected_nodes(azide_site.id, [azide_site.id])
        find_connected_nodes(alkyne_site.id, [alkyne_site.id])

        # Keep the reactive sites themselves but remove everything else
        nodes_to_remove = [n for n in nodes_to_remove if n.id != azide_site.id and n.id != alkyne_site.id]

        # Remove all collected nodes and edges
        cyclic.nodes = [n for n in cyclic.nodes if not any(remove_node.id == n.id for remove_node in nodes_to_remove)]
        cyclic.edges = [e for e in cyclic.edges if not any(
            remove_edge.from_idx == e.from_idx and remove_edge.to_idx == e.to_idx
            for remove_edge in edges_to_remove
        )]

        # Step 2: Create new atoms for triazole
        max_id = max(n.id for n in cyclic.nodes)

        # Create new ring carbon (position 4)
        c4 = GraphNode(
            id=max_id + 1,
            element='C',
            atomic_num=6,
            formal_charge=0,
            implicit_valence=3,
            explicit_valence=3,
            aromatic=True,
            hybridization='SP2',
            num_explicit_hs=0,
            num_implicit_hs=1,
            total_num_hs=1,
            degree=2,
            in_ring=True
        )

        # Create middle nitrogen (position 2)
        n2 = GraphNode(
            id=max_id + 2,
            element='N',
            atomic_num=7,
            formal_charge=0,
            implicit_valence=2,
            explicit_valence=2,
            aromatic=True,
            hybridization='SP2',
            num_explicit_hs=0,
            num_implicit_hs=0,
            total_num_hs=0,
            degree=2,
            in_ring=True
        )

        # Update alkyne site to ring carbon (position 5)
        for node in cyclic.nodes:
            if node.id == alkyne_site.id:
                node.is_reactive_elec = False
                node.element = 'C'
                node.atomic_num = 6
                node.formal_charge = 0
                node.implicit_valence = 3
                node.explicit_valence = 3
                node.aromatic = True
                node.hybridization = 'SP2'
                node.num_explicit_hs = 0
                node.num_implicit_hs = 0
                node.total_num_hs = 0
                node.degree = 2
                node.in_ring = True

            elif node.id == azide_site.id:
                # Convert to NH nitrogen (position 1)
                node.is_reactive_nuc = False
                node.element = 'N'
                node.atomic_num = 7
                node.formal_charge = 0
                node.implicit_valence = 2
                node.explicit_valence = 2
                node.aromatic = True
                node.hybridization = 'SP2'
                node.num_explicit_hs = 0
                node.num_implicit_hs = 0
                node.total_num_hs = 0
                node.degree = 2
                node.in_ring = True

        # Add new atoms
        cyclic.nodes.extend([c4, n2])

        # Create triazole ring bonds
        triazole_bonds = [
            GraphEdge(
                from_idx=azide_site.id,
                to_idx=alkyne_site.id,
                bond_type='SINGLE',
                is_aromatic=True,
                is_conjugated=True,
                in_ring=True
            ),
            GraphEdge(
                from_idx=alkyne_site.id,
                to_idx=n2.id,
                bond_type='SINGLE',
                is_aromatic=True,
                is_conjugated=True,
                in_ring=True
            ),
            GraphEdge(
                from_idx=n2.id,
                to_idx=c4.id,
                bond_type='SINGLE',
                is_aromatic=True,
                is_conjugated=True,
                in_ring=True
            ),
            GraphEdge(
                from_idx=c4.id,
                to_idx=azide_site.id,
                bond_type='SINGLE',
                is_aromatic=True,
                is_conjugated=True,
                in_ring=True
            )
        ]

        cyclic.edges.extend(triazole_bonds)

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
        """Check if node is the first nitrogen of an azide group."""
        if node.element != 'N':
            return False

        # Check for N-N≡N pattern
        connected = self._get_connected_atoms(node.id, graph)
        if len(connected) != 2:  # Should have one bond to carbon chain and one to next nitrogen
            return False

        n2 = next((n for n in connected if n.element == 'N'), None)
        if not n2:
            return False

        n2_connected = self._get_connected_atoms(n2.id, graph)
        n3 = next((n for n in n2_connected if n.element == 'N' and n.id != node.id), None)

        return n3 is not None

    def _is_alkyne_carbon(self, node: GraphNode, graph: MolecularGraph) -> bool:
        """Check if node is part of an alkyne group."""
        if node.element != 'C':
            return False

        # Find triple bond
        triple_bond = next(
            (e for e in graph.edges
                if (e.from_idx == node.id or e.to_idx == node.id)
                and e.bond_type == 'TRIPLE'),
            None
        )

        if not triple_bond:
            return False

        # Find the other carbon of the triple bond
        other_carbon_id = triple_bond.to_idx if triple_bond.from_idx == node.id else triple_bond.from_idx
        other_carbon = next((n for n in graph.nodes if n.id == other_carbon_id), None)

        if not other_carbon or other_carbon.element != 'C':
            return False

        # Count connections to both carbons
        node_connections = len(self._get_connected_atoms(node.id, graph))
        other_connections = len(self._get_connected_atoms(other_carbon_id, graph))

        # At least one of the carbons should have only one connection besides the triple bond
        return node_connections <= 2 or other_connections <= 2

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
        """Cyclize linear peptide."""
        nuc_site = next(
            (n for n in linear_peptide.nodes if n.is_reactive_nuc),
            None
        )
        elec_site = next(
            (n for n in linear_peptide.nodes if n.is_reactive_elec),
            None
        )

        if not (nuc_site and elec_site):
            raise ValueError("Cannot find terminal reactive sites")

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
