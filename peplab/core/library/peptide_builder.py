# @title Peptide Builder

import logging
from itertools import product
import copy
import concurrent.futures
from tqdm import tqdm
import threading
from typing import List, Optional
from queue import Queue
from uuid import uuid4

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
            print(":) 1")
            residue_graphs = [res['graph'] for res in combo]
            peptide = self.build_linear_peptide(residue_graphs)
            print(":) 2")
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
            print(":) 0")
            return PeptideInfo(
                graph=peptide.to_dict(),
                id=str(uuid4()),
                sequence=[res['name'] for res in combo],
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
        self._cyclize = cyclize or click_cyclize
        self._click_cyclize = click_cyclize

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
                            'name': res_info.name,
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
            print(":) 3")
            peptide = self._form_peptide_bond(peptide, residue_graphs[i])

        return peptide

    def _form_peptide_bond(self, res1: MolecularGraph, res2: MolecularGraph) -> MolecularGraph:
        """Form peptide bond between residues, handling both primary and secondary amines."""
        nuc_site = next(
            (n for n in res1.nodes if n.is_reactive_nuc and not n.is_reactive_click),
            None
        )
        elec_site = next(
            (n for n in res2.nodes if n.is_reactive_elec and not n.is_reactive_click),
            None
        )

        if not (nuc_site and elec_site):
            self.logger.warning(
                f"Failed to find reactive sites - Amine: {nuc_site is not None}, "
                f"COOH: {elec_site is not None}"
            )
            raise ValueError("Could not find required reactive sites for peptide bond formation")

        res1_mod = copy.deepcopy(res1)
        res2_mod = copy.deepcopy(res2)

        # Remove appropriate number of hydrogens based on amine type
        res1_mod = self._remove_h_from_nh(res1_mod, nuc_site.index )
        res2_mod = self._remove_oh_from_cooh(res2_mod, elec_site.index )

        # Combine graphs with updated atom indexing
        max_index = max(n.index  for n in res1_mod.nodes) + 1
        for node in res2_mod.nodes:
            node.index  += max_index
        for edge in res2_mod.edges:
            edge.from_idx += max_index
            edge.to_idx += max_index

        combined = MolecularGraph()
        combined.nodes = res1_mod.nodes + res2_mod.nodes
        combined.edges = res1_mod.edges + res2_mod.edges

        # Form the new peptide bond
        combined.edges.append(GraphEdge(
            id=str(uuid4()),
            from_idx=nuc_site.index,
            to_idx=elec_site.index  + max_index,
            bond_type='SINGLE',
            is_aromatic=False,
            is_conjugated=False,
            in_ring=False,
            stereo='NONE'
        ))

        # Update reactive sites
        for node in combined.nodes:
            if node.index  == nuc_site.index :
                node.is_reactive_nuc = False
            if node.index  == elec_site.index  + max_index:
                node.is_reactive_elec = False
        print(":) 4")
        return self._reindex_graph(combined)

    def _get_amine_type(self, node: GraphNode, graph: MolecularGraph) -> str:
        """
        Determine if a nitrogen is a primary or secondary amine.
        Returns 'NH2' for primary or 'NH' for secondary amine.
        """
        if not node.is_reactive_nuc:
            return ''

        h_count = 0
        heavy_count = 0

        # Count connected atoms
        for neighbor, edge in graph.get_neighbors(node.index ):
            if neighbor.element == 'H':
                h_count += 1
            else:
                heavy_count += 1

        if h_count == 2 and heavy_count == 1:
            return 'NH2'
        elif h_count == 1 and heavy_count == 2:
            return 'NH'
        return ''

    def _remove_h_from_nh(self, graph: MolecularGraph, node_id: int) -> MolecularGraph:
        """
        Remove appropriate number of hydrogens based on amine type.
        For NH2: removes one H
        For NH: removes the single H
        """
        mod_graph = copy.deepcopy(graph)
        amine_type = self._get_amine_type(
            next(n for n in mod_graph.nodes if n.index  == node_id),
            mod_graph
        )

        h_to_remove = 1  # Remove one H for either type
        h_removed = 0

        edges_to_remove = []
        nodes_to_remove = []

        for edge in mod_graph.edges:
            if h_removed >= h_to_remove:
                break

            if edge.from_idx == node_id:
                h_node = next((n for n in mod_graph.nodes
                              if n.index  == edge.to_idx and n.element == 'H'), None)
                if h_node:
                    edges_to_remove.append(edge)
                    nodes_to_remove.append(h_node)
                    h_removed += 1
            elif edge.to_idx == node_id:
                h_node = next((n for n in mod_graph.nodes
                              if n.index  == edge.from_idx and n.element == 'H'), None)
                if h_node:
                    edges_to_remove.append(edge)
                    nodes_to_remove.append(h_node)
                    h_removed += 1

        for edge in edges_to_remove:
            if edge in mod_graph.edges:
                mod_graph.edges.remove(edge)
        for node in nodes_to_remove:
            if node in mod_graph.nodes:
                mod_graph.nodes.remove(node)

        return mod_graph

    def _remove_oh_from_cooh(self, graph: MolecularGraph, node_id: int) -> MolecularGraph:
        """Remove OH group from COOH."""
        mod_graph = copy.deepcopy(graph)

        for edge in mod_graph.edges[:]:
            if edge.bond_type == 'SINGLE':
                if edge.from_idx == node_id:
                    o_node = next((n for n in mod_graph.nodes if n.index  == edge.to_idx and n.element == 'O'), None)
                elif edge.to_idx == node_id:
                    o_node = next((n for n in mod_graph.nodes if n.index  == edge.from_idx and n.element == 'O'), None)
                else:
                    continue

                if o_node:
                    # Find H attached to O
                    for h_edge in mod_graph.edges[:]:
                        if h_edge.from_idx == o_node.index :
                            h_node = next((n for n in mod_graph.nodes if n.index  == h_edge.to_idx and n.element == 'H'), None)
                        elif h_edge.to_idx == o_node.index :
                            h_node = next((n for n in mod_graph.nodes if n.index  == h_edge.from_idx and n.element == 'H'), None)
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
        # Debug logging
        self.logger.debug("Checking for click pairs...")

        # Find residues with potential azide/alkyne groups
        azide_residues = []
        alkyne_residues = []

        for res in combo:
            # Debug logging for each residue
            self.logger.debug(f"Checking residue: {res.get('id')}")
            self.logger.debug(f"Nuc: {res.get('nuc')}, Elec: {res.get('elec')}")

            if res.get('nuc') == 'N3' or res.get('elec') == 'N3':
                azide_residues.append(res)
                self.logger.debug(f"Found azide residue: {res.get('id')}")
            if res.get('nuc') == 'C#C' or res.get('elec') == 'C#C':
                alkyne_residues.append(res)
                self.logger.debug(f"Found alkyne residue: {res.get('id')}")

        if not (azide_residues and alkyne_residues):
            self.logger.debug("No azide/alkyne pairs found")
            return False

        # For each potential pair, validate the reactive sites
        for azide_res in azide_residues:
            graph = azide_res['graph']
            # Try to find any nitrogen that's part of an azide group
            for node in graph.nodes:
                if node.element == 'N' and self._is_azide_nitrogen(node, graph):
                    self.logger.debug(f"Found valid azide site in {azide_res.get('index')}")
                    azide_n = node

                    for alkyne_res in alkyne_residues:
                        if azide_res == alkyne_res:  # Skip same residue
                            continue

                        graph = alkyne_res['graph']
                        alkyne_c = next((n for n in graph.nodes if n.is_reactive_click and
                                       self._is_alkyne_carbon(n, graph)), None)

                        if alkyne_c:
                            self.logger.debug(f"Found valid alkyne site in {alkyne_res.get('id')}")
                            return True
                        else:
                            self.logger.debug(f"No valid alkyne site found in {alkyne_res.get('id')}")

            self.logger.debug(f"No valid azide site found in {azide_res.get('id')}")

        self.logger.debug("No valid click pairs found after checking all combinations")
        return False

    def click_cyclize_peptide(self, linear_peptide: MolecularGraph) -> MolecularGraph:
        """
        Form triazole ring through copper-catalyzed azide-alkyne cycloaddition (CuAAC).
        Creates a 1,4-disubstituted 1,2,3-triazole ring.
        Handles both N-N≡N and N=[N+]=[N-] azide representations.
        """
        # Track initial state
        initial_atom_count = len(linear_peptide.nodes)
        self.logger.debug(f"Initial atom count: {initial_atom_count}")

        cyclic = copy.deepcopy(linear_peptide)

        def count_atoms_by_element(graph):
            """Helper to count atoms by element"""
            counts = {}
            for node in graph.nodes:
                counts[node.element] = counts.get(node.element, 0) + 1
            return counts

        initial_counts = count_atoms_by_element(linear_peptide)
        self.logger.debug(f"Initial atom counts by element: {initial_counts}")

        # Find azide nitrogen chain
        azide_start = None
        azide_chain = []

        # Look for any nitrogen that's part of an azide group
        for node in cyclic.nodes:
            if node.element == 'N' and self._is_azide_nitrogen(node, cyclic):
                azide_start = node
                break

        if not azide_start:
            raise ValueError("Could not find azide group starting nitrogen")

        # Get connected nitrogens to find full azide chain
        visited = {azide_start.index }
        current = azide_start
        azide_chain = [current]

        while len(azide_chain) < 3:
            next_n = None
            for edge in cyclic.edges:
                if edge.from_idx == current.index :
                    node = next((n for n in cyclic.nodes if n.index  == edge.to_idx
                               and n.element == 'N' and n.index  not in visited), None)
                    if node:
                        next_n = node
                        break
                elif edge.to_idx == current.index :
                    node = next((n for n in cyclic.nodes if n.index  == edge.from_idx
                               and n.element == 'N' and n.index  not in visited), None)
                    if node:
                        next_n = node
                        break

            if not next_n:
                break

            visited.add(next_n.index )
            azide_chain.append(next_n)
            current = next_n

        # Find alkyne carbons
        alkyne_chain = []
        for node in cyclic.nodes:
            if node.element == 'C' and self._is_alkyne_carbon(node, cyclic):
                alkyne_chain = [node]
                # Find the connected triple-bonded carbon
                for edge in cyclic.edges:
                    if edge.bond_type == 'TRIPLE':
                        if edge.from_idx == node.index:
                            other_c = next((n for n in cyclic.nodes if n.index == edge.to_idx
                                          and n.element == 'C'), None)
                            if other_c:
                                alkyne_chain.append(other_c)
                                break
                        elif edge.to_idx == node.index:
                            other_c = next((n for n in cyclic.nodes if n.index == edge.from_idx
                                          and n.element == 'C'), None)
                            if other_c:
                                alkyne_chain.append(other_c)
                                break
                if len(alkyne_chain) == 2:
                    break

        if len(azide_chain) != 3 or len(alkyne_chain) != 2:
            raise ValueError("Could not find complete azide or alkyne groups")

        self.logger.debug(f"Found azide chain with {len(azide_chain)} N atoms")
        self.logger.debug(f"Found alkyne chain with {len(alkyne_chain)} C atoms")

        # Track atoms to be removed
        atoms_to_remove = set()
        edges_to_remove = []

        # Remove existing bonds between azide nitrogens while keeping track of removed atoms
        for edge in cyclic.edges:
            from_in_azide = any(n.index  == edge.from_idx for n in azide_chain)
            to_in_azide = any(n.index  == edge.to_idx for n in azide_chain)
            from_in_alkyne = any(n.index  == edge.from_idx for n in alkyne_chain)
            to_in_alkyne = any(n.index  == edge.to_idx for n in alkyne_chain)

            if (from_in_azide and to_in_azide) or (from_in_alkyne and to_in_alkyne):
                edges_to_remove.append(edge)

                # If this is a bond to a hydrogen, mark the hydrogen for removal
                for node_id in [edge.from_idx, edge.to_idx]:
                    node = next((n for n in cyclic.nodes if n.index  == node_id), None)
                    if node and node.element == 'H':
                        atoms_to_remove.add(node)

        # Remove tracked edges
        for edge in edges_to_remove:
            if edge in cyclic.edges:
                cyclic.edges.remove(edge)

        # Remove hydrogens from reaction sites while tracking removals
        atoms_to_clean = azide_chain + alkyne_chain
        for atom in atoms_to_clean:
            edges_to_remove = []
            for edge in cyclic.edges:
                if edge.from_idx == atom.index  or edge.to_idx == atom.index :
                    other_id = edge.to_idx if edge.from_idx == atom.index  else edge.from_idx
                    other_node = next((n for n in cyclic.nodes if n.index  == other_id), None)
                    if other_node and other_node.element == 'H':
                        edges_to_remove.append(edge)
                        atoms_to_remove.add(other_node)

            for edge in edges_to_remove:
                if edge in cyclic.edges:
                    cyclic.edges.remove(edge)

        # Remove all tracked atoms
        for node in atoms_to_remove:
            if node in cyclic.nodes:
                cyclic.nodes.remove(node)

        # Log intermediate state
        intermediate_counts = count_atoms_by_element(cyclic)
        self.logger.debug(f"Atom counts after cleanup: {intermediate_counts}")

        # Create triazole ring with correct alternating single/double bonds
        new_edges = [
            GraphEdge(  # N1-N2
                id=str(uuid4()),
                from_idx=azide_chain[0].index,
                to_idx=azide_chain[1].index,
                bond_type='SINGLE',
                is_aromatic=True,
                is_conjugated=True,
                in_ring=True,
                stereo='NONE'
            ),
            GraphEdge(  # N2=N3
                id=str(uuid4()),
                from_idx=azide_chain[1].index,
                to_idx=azide_chain[2].index,
                bond_type='DOUBLE',
                is_aromatic=True,
                is_conjugated=True,
                in_ring=True,
                stereo='NONE'
            ),
            GraphEdge(  # N3-C4
                id=str(uuid4()),
                from_idx=azide_chain[2].index,
                to_idx=alkyne_chain[0].index,
                bond_type='SINGLE',
                is_aromatic=True,
                is_conjugated=True,
                in_ring=True,
                stereo='NONE'
            ),
            GraphEdge(  # C4=C5
                id=str(uuid4()),
                from_idx=alkyne_chain[0].index,
                to_idx=alkyne_chain[1].index,
                bond_type='DOUBLE',
                is_aromatic=True,
                is_conjugated=True,
                in_ring=True,
                stereo='NONE'
            ),
            GraphEdge(  # C5-N1
                id=str(uuid4()),
                from_idx=alkyne_chain[1].index,
                to_idx=azide_chain[0].index,
                bond_type='SINGLE',
                is_aromatic=True,
                is_conjugated=True,
                in_ring=True,
                stereo='NONE'
            )
        ]

        # Add new bonds
        for new_edge in new_edges:
            cyclic.edges.append(new_edge)

        # Update properties for atoms in the triazole ring
        for node in cyclic.nodes:
            if node in azide_chain or node in alkyne_chain:
                node.is_reactive_nuc = False
                node.is_reactive_elec = False
                node.is_reactive_click = False
                node.num_implicit_hs = 0
                node.num_explicit_hs = 0
                node.formal_charge = 0
                node.is_aromatic = True
                if node.element == 'N':
                    node.explicit_valence = 3
                elif node.element == 'C':
                    node.explicit_valence = 4

        # Verify final atom count
        final_graph = self._reindex_graph(cyclic)
        final_counts = count_atoms_by_element(final_graph)
        self.logger.debug(f"Final atom counts by element: {final_counts}")

        # Validate atom conservation
        if len(final_graph.nodes) != initial_atom_count - 2:  # We should lose 2 hydrogens in the click reaction
            self.logger.error(f"Atom count mismatch! Initial: {initial_atom_count}, Final: {len(final_graph.nodes)}")
            self.logger.error("Initial counts: " + str(initial_counts))
            self.logger.error("Final counts: " + str(final_counts))
            raise ValueError("Atom count mismatch in click reaction")

        return final_graph

    def _get_azide_chain(self, graph: MolecularGraph, start_id: int) -> List[GraphNode]:
        """Get the complete azide chain starting from the first nitrogen."""
        chain = [next(n for n in graph.nodes if n.index  == start_id)]
        current_id = start_id
        visited = {current_id}

        while True:
            next_n = None
            for edge in graph.edges:
                if edge.from_idx == current_id or edge.to_idx == current_id:
                    other_id = edge.to_idx if edge.from_idx == current_id else edge.from_idx
                    node = next((n for n in graph.nodes if n.index  == other_id), None)
                    if node and node.element == 'N' and node.index  not in visited:
                        next_n = node
                        chain.append(node)
                        visited.add(node.index )
                        current_id = node.index
                        break
            if not next_n:
                break

        return chain

    def _get_alkyne_chain(self, graph: MolecularGraph, start_id: int) -> List[GraphNode]:
        """Get the complete alkyne chain starting from the first carbon."""
        chain = [next(n for n in graph.nodes if n.index  == start_id)]

        # Find the carbon connected by triple bond
        for edge in graph.edges:
            if edge.bond_type == 'TRIPLE' and (edge.from_idx == start_id or edge.to_idx == start_id):
                other_id = edge.to_idx if edge.from_idx == start_id else edge.from_idx
                other_node = next((n for n in graph.nodes if n.index  == other_id), None)
                if other_node and other_node.element == 'C':
                    chain.append(other_node)
                    break

        return chain

    def _is_azide_nitrogen(self, node: GraphNode, graph: MolecularGraph) -> bool:
        """
        Check if node is part of an azide group.
        Handles both N-N≡N and N=[N+]=[N-] representations.
        """
        if node.element != 'N':
            return False

        # Get connected atoms and their bonds
        connected = []
        for edge in graph.edges:
            if edge.from_idx == node.index :
                atom = next(n for n in graph.nodes if n.index  == edge.to_idx)
                connected.append((atom, edge))
            elif edge.to_idx == node.index :
                atom = next(n for n in graph.nodes if n.index  == edge.from_idx)
                connected.append((atom, edge))

        # Check for ionic azide representation (N=[N+]=[N-])
        def check_ionic_azide():
            # First N should have one or two connections (one to next N, possibly one to carbon chain)
            if not (1 <= len(connected) <= 2):
                return False

            # Find the second nitrogen (N2) connected by double bond
            n2_pair = next(((atom, edge) for atom, edge in connected
                           if atom.element == 'N' and edge.bond_type == 'DOUBLE'), None)
            if not n2_pair:
                return False

            n2, n1_n2_bond = n2_pair

            # Check N2 connections
            n2_connected = []
            for edge in graph.edges:
                if edge.from_idx == n2.index :
                    atom = next(n for n in graph.nodes if n.index  == edge.to_idx)
                    n2_connected.append((atom, edge))
                elif edge.to_idx == n2.index :
                    atom = next(n for n in graph.nodes if n.index  == edge.from_idx)
                    n2_connected.append((atom, edge))

            # N2 should have exactly two connections
            if len(n2_connected) != 2:
                return False

            # Find N3 (should be connected to N2 with a double bond)
            n3_pair = next(((atom, edge) for atom, edge in n2_connected
                           if atom.element == 'N' and atom.index  != node.index
                           and edge.bond_type == 'DOUBLE'), None)
            if not n3_pair:
                return False

            n3, n2_n3_bond = n3_pair

            # N3 should only have one connection (to N2)
            n3_connections = sum(1 for edge in graph.edges
                               if edge.from_idx == n3.index  or edge.to_idx == n3.index )
            return n3_connections == 1

        # Check for covalent azide representation (N-N≡N)
        def check_covalent_azide():
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
                if edge.from_idx == n2.index :
                    atom = next(n for n in graph.nodes if n.index  == edge.to_idx)
                    n2_connected.append((atom, edge))
                elif edge.to_idx == n2.index :
                    atom = next(n for n in graph.nodes if n.index  == edge.from_idx)
                    n2_connected.append((atom, edge))

            # N2 should have exactly two connections
            if len(n2_connected) != 2:
                return False

            # Find N3 (should be connected to N2 with a triple bond)
            n3_pair = next(((atom, edge) for atom, edge in n2_connected
                           if atom.element == 'N' and atom.index  != node.index
                           and edge.bond_type == 'TRIPLE'), None)
            if not n3_pair:
                return False

            n3, n2_n3_bond = n3_pair

            # N3 should only have one connection (to N2)
            n3_connections = sum(1 for edge in graph.edges
                               if edge.from_idx == n3.index  or edge.to_idx == n3.index )
            return n3_connections == 1

        # Try both representations
        return check_ionic_azide() or check_covalent_azide()

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
            if edge.from_idx == node.index :
                atom = next(n for n in graph.nodes if n.index  == edge.to_idx)
                connected.append((atom, edge))
            elif edge.to_idx == node.index :
                atom = next(n for n in graph.nodes if n.index  == edge.from_idx)
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
                if edge.from_idx == other_c.index :
                    atom = next(n for n in graph.nodes if n.index  == edge.to_idx)
                    other_c_connections.append((atom, edge))
                elif edge.to_idx == other_c.index :
                    atom = next(n for n in graph.nodes if n.index  == edge.from_idx)
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
                connected.append(next(n for n in graph.nodes if n.index  == edge.to_idx))
            elif edge.to_idx == node_id:
                connected.append(next(n for n in graph.nodes if n.index  == edge.from_idx))
        return connected

    def _get_azide_nitrogens(self, graph: MolecularGraph, start_nitrogen_id: int) -> List[GraphNode]:
        """Get all nitrogen atoms in an azide group starting from the first nitrogen."""
        nitrogens = [next(n for n in graph.nodes if n.index  == start_nitrogen_id)]
        current = nitrogens[0]

        while True:
            connected = self._get_connected_atoms(current.index , graph)
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
        cyclic = self._remove_h_from_nh(cyclic, nuc_site.index )
        cyclic = self._remove_oh_from_cooh(cyclic, elec_site.index )
        cyclic.edges.append(GraphEdge(
            id=str(uuid4),
            from_idx=nuc_site.index ,
            to_idx=elec_site.index ,
            bond_type='SINGLE',
            is_aromatic=False,
            is_conjugated=False,
            in_ring=True,
            stereo='NONE'
        ))

        for node in cyclic.nodes:
            if node.index  == nuc_site.index :
                node.is_reactive_nuc = False
            if node.index  == elec_site.index :
                node.is_reactive_elec = False

        return self._reindex_graph(cyclic)

    def _reindex_graph(self, graph: MolecularGraph) -> MolecularGraph:
        """Reindex graph nodes and edges sequentially."""
        new_graph = MolecularGraph()
        old_to_new = {}

        # Reindex nodes
        for i, node in enumerate(sorted(graph.nodes, key=lambda x: x.index )):
            new_node = copy.deepcopy(node)
            old_to_new[node.index] = i
            new_graph.nodes.append(new_node)
        print(":) 5")
        # Update edges
        for edge in graph.edges:
            new_edge = copy.deepcopy(edge)
            new_edge.from_idx = old_to_new[edge.from_idx]
            new_edge.to_idx = old_to_new[edge.to_idx]
            new_graph.edges.append(new_edge)
        print(":) 6")
        return new_graph
