# @title Peptide Builder

"""Peptide builder with improved error handling."""

import logging
from typing import List, Dict, Optional
from dataclasses import dataclass
import itertools
import copy

from ..graph.molecule_graph import MolecularGraph
from ..molecule.atom import GraphNode
from ..molecule.bond import GraphEdge
from ..molecule.peptide import PeptideInfo
from .library import LibraryInfo

class PeptideBuilder:
    """Handles construction of peptides."""

    def __init__(self):
        self.logger = logging.getLogger(self.__class__.__name__)

    def build_linear_peptide(self, residue_graphs: List[MolecularGraph]) -> MolecularGraph:
        """Build linear peptide exactly as in old code."""
        if not residue_graphs:
            raise ValueError("No residue graphs provided")

        peptide = copy.deepcopy(residue_graphs[0])

        for i in range(1, len(residue_graphs)):
            peptide = self._form_peptide_bond(peptide, residue_graphs[i])

        return peptide

    def _form_peptide_bond(self, graph1: MolecularGraph, graph2: MolecularGraph) -> MolecularGraph:
        """Form peptide bond matching old code exactly."""
        # Find reactive sites
        nuc_site = next((n for n in graph1.nodes if n.is_reactive_nuc), None)
        elec_site = next((n for n in graph2.nodes if n.is_reactive_elec), None)

        if not (nuc_site and elec_site):
            raise ValueError("Could not find required reactive sites")

        # Make copies to modify
        g1 = copy.deepcopy(graph1)
        g2 = copy.deepcopy(graph2)

        # Remove reactive groups
        g1 = self._remove_h_from_nh(g1, nuc_site.id)
        g2 = self._remove_oh_from_cooh(g2, elec_site.id)

        # Adjust indices for second graph
        max_id = max(n.id for n in g1.nodes) + 1
        for node in g2.nodes:
            node.id += max_id
        for edge in g2.edges:
            edge.from_idx += max_id
            edge.to_idx += max_id

        # Combine graphs
        combined = MolecularGraph()
        combined.nodes = g1.nodes + g2.nodes
        combined.edges = g1.edges + g2.edges

        # Add peptide bond
        peptide_bond = GraphEdge(
            from_idx=nuc_site.id,
            to_idx=elec_site.id + max_id,
            bond_type='SINGLE',
            is_aromatic=False,
            is_conjugated=False,
            in_ring=False,
            stereo='NONE'
        )
        combined.edges.append(peptide_bond)

        # Update reactive flags
        nuc_node = next(n for n in combined.nodes if n.id == nuc_site.id)
        elec_node = next(n for n in combined.nodes if n.id == elec_site.id + max_id)
        nuc_node.is_reactive_nuc = False
        elec_node.is_reactive_elec = False

        return self._reindex_graph(combined)

    def _remove_h_from_nh(self, graph: MolecularGraph, n_id: int) -> MolecularGraph:
        """Remove H from NH/NH2 exactly as in old code."""
        mod_graph = copy.deepcopy(graph)

        for edge in mod_graph.edges:
            if edge.from_idx == n_id:
                h_node = next((n for n in mod_graph.nodes if n.id == edge.to_idx), None)
                if h_node and h_node.element == 'H':
                    mod_graph.nodes.remove(h_node)
                    mod_graph.edges.remove(edge)
                    break
            elif edge.to_idx == n_id:
                h_node = next((n for n in mod_graph.nodes if n.id == edge.from_idx), None)
                if h_node and h_node.element == 'H':
                    mod_graph.nodes.remove(h_node)
                    mod_graph.edges.remove(edge)
                    break

        return mod_graph

    def _remove_oh_from_cooh(self, graph: MolecularGraph, c_id: int) -> MolecularGraph:
        """Remove OH from COOH exactly as in old code."""
        mod_graph = copy.deepcopy(graph)

        # Find O with single bond to C
        for edge in mod_graph.edges:
            if edge.bond_type == 'SINGLE':
                if edge.from_idx == c_id:
                    o_node = next((n for n in mod_graph.nodes if n.id == edge.to_idx), None)
                elif edge.to_idx == c_id:
                    o_node = next((n for n in mod_graph.nodes if n.id == edge.from_idx), None)
                else:
                    continue

                if o_node and o_node.element == 'O':
                    # Find H attached to this O
                    for h_edge in mod_graph.edges:
                        if h_edge.from_idx == o_node.id:
                            h_node = next((n for n in mod_graph.nodes if n.id == h_edge.to_idx), None)
                        elif h_edge.to_idx == o_node.id:
                            h_node = next((n for n in mod_graph.nodes if n.id == h_edge.from_idx), None)
                        else:
                            continue

                        if h_node and h_node.element == 'H':
                            # Remove both O and H
                            mod_graph.nodes.remove(o_node)
                            mod_graph.nodes.remove(h_node)
                            mod_graph.edges.remove(edge)
                            mod_graph.edges.remove(h_edge)
                            return mod_graph

        return mod_graph

    def _reindex_graph(self, graph: MolecularGraph) -> MolecularGraph:
        """Reindex graph nodes and edges exactly as in old code."""
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

    def enumerate_library(self, library_info: LibraryInfo, cyclize: bool = False) -> List[PeptideInfo]:
        """Enumerate all possible peptides with exact old code functionality."""
        peptides = []
        total_attempts = 0
        successful = 0

        try:
            # Get positions in sorted order
            positions = sorted(library_info.positions.keys())
            residue_options = []

            # Process each position
            for pos in positions:
                pos_residues = []
                for res_key, res_info in library_info.positions[pos].residues.items():
                    try:
                        # Create and mark residue graph
                        graph = MolecularGraph.from_smiles(res_info.smiles)
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

            # Generate all possible combinations
            for combo in itertools.product(*residue_options):
                total_attempts += 1
                try:
                    # Get residue graphs for this combination
                    residue_graphs = [res['graph'] for res in combo]

                    # Build linear peptide
                    peptide = self.build_linear_peptide(residue_graphs)

                    if cyclize:
                        try:
                            peptide = self.cyclize_peptide(peptide)
                            peptide_type = 'cyclic'
                        except Exception as e:
                            self.logger.warning(
                                f"Error cyclizing peptide combination {total_attempts}: {str(e)}"
                            )
                            continue
                    else:
                        peptide_type = 'linear'

                    # Create PeptideInfo object
                    peptides.append(PeptideInfo(
                        graph=peptide.to_dict(),
                        sequence=[res['id'] for res in combo],
                        peptide_type=peptide_type,
                        smiles=None  # Will be generated later
                    ))
                    successful += 1

                except Exception as e:
                    self.logger.warning(
                        f"Error generating peptide combination {total_attempts}: {str(e)}"
                    )

            self.logger.info(
                f"Generated {successful} out of {total_attempts} possible peptides"
            )
            return peptides

        except Exception as e:
            self.logger.error(f"Error enumerating library: {str(e)}")
            return []

    def cyclize_peptide(self, linear_peptide: MolecularGraph) -> MolecularGraph:
        """Cyclize linear peptide exactly as in old code."""
        # Find terminal reactive sites
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

        # Make copy to modify
        cyclic = copy.deepcopy(linear_peptide)

        # Remove terminal groups
        cyclic = self._remove_h_from_nh(cyclic, nuc_site.id)
        cyclic = self._remove_oh_from_cooh(cyclic, elec_site.id)

        # Add cyclizing bond
        cyclic.edges.append(GraphEdge(
            from_idx=nuc_site.id,
            to_idx=elec_site.id,
            bond_type='SINGLE',
            is_aromatic=False,
            is_conjugated=False,
            in_ring=True,
            stereo='NONE'
        ))

        # Update reactive flags
        for node in cyclic.nodes:
            if node.id == nuc_site.id:
                node.is_reactive_nuc = False
            if node.id == elec_site.id:
                node.is_reactive_elec = False

        return self._reindex_graph(cyclic)
