# peplab/chemistry/strategies/bonding/amide_bonding.py

"""
Module for amide bonding strategy.

This module provides the AmideBonding class which implements the BondingStrategy protocol
for forming amide bonds between molecular graphs.
"""

import copy
from ....core.graph import Edge, MolecularGraph, Node
from .. import BondingStrategy

class AmideBonding(BondingStrategy):
    """Concrete strategy for forming amide bonds."""

    def form_bond(self, res1: MolecularGraph, res2: MolecularGraph) -> MolecularGraph:
        """
        Form an amide bond between two molecular graphs.

        Args:
            res1 (MolecularGraph): The first molecular graph.
            res2 (MolecularGraph): The second molecular graph.

        Returns:
            MolecularGraph: The molecular graph with the newly formed amide bond.
        """
        nuc_site = next((n for n in res1.nodes if n.is_reactive_nuc and not n.is_reactive_click), None)
        elec_site = next((n for n in res2.nodes if n.is_reactive_elec and not n.is_reactive_click), None)

        if not (nuc_site and elec_site):
            raise ValueError("Could not find required reactive sites for peptide bond formation")

        res1_mod = copy.deepcopy(res1)
        res2_mod = copy.deepcopy(res2)

        res1_mod = self._remove_h_from_nh(res1_mod, nuc_site.index)
        res2_mod = self._remove_oh_from_cooh(res2_mod, elec_site.index)

        max_index = max(n.index for n in res1_mod.nodes) + 1
        for node in res2_mod.nodes:
            node.index += max_index
        for edge in res2_mod.edges:
            edge.from_idx += max_index
            edge.to_idx += max_index

        combined = MolecularGraph()
        combined.nodes = res1_mod.nodes + res2_mod.nodes
        combined.edges = res1_mod.edges + res2_mod.edges

        combined.edges.append(Edge(
            from_idx=nuc_site.index,
            to_idx=elec_site.index + max_index,
            properties={'bond_type': 'SINGLE'}
        ))

        for node in combined.nodes:
            if node.index == nuc_site.index:
                node.is_reactive_nuc = False
            if node.index == elec_site.index + max_index:
                node.is_reactive_elec = False

        return self._reindex_graph(combined)

    def _remove_h_from_nh(self, graph: MolecularGraph, node_id: int) -> MolecularGraph:
        """
        Remove appropriate number of hydrogens based on amine type.

        Args:
            graph (MolecularGraph): The molecular graph.
            node_id (int): The index of the nitrogen atom to modify.

        Returns:
            MolecularGraph: The modified molecular graph.
        """
        mod_graph = copy.deepcopy(graph)
        node = next(n for n in mod_graph.nodes if n.index == node_id)

        h_to_remove = 1  # For simplicity, assume removing one H
        h_removed = 0

        edges_to_remove = []
        nodes_to_remove = []

        for edge in mod_graph.edges:
            if h_removed >= h_to_remove:
                break
            if edge.from_idx == node_id:
                h_node = next((n for n in mod_graph.nodes if n.index == edge.to_idx and n.element == 'H'), None)
                if h_node:
                    edges_to_remove.append(edge)
                    nodes_to_remove.append(h_node)
                    h_removed += 1
            elif edge.to_idx == node_id:
                h_node = next((n for n in mod_graph.nodes if n.index == edge.from_idx and n.element == 'H'), None)
                if h_node:
                    edges_to_remove.append(edge)
                    nodes_to_remove.append(h_node)
                    h_removed += 1

        for edge in edges_to_remove:
            mod_graph.edges.remove(edge)
        for node in nodes_to_remove:
            mod_graph.nodes.remove(node)

        return mod_graph

    def _remove_oh_from_cooh(self, graph: MolecularGraph, node_id: int) -> MolecularGraph:
        """
        Remove OH group from COOH.

        Args:
            graph (MolecularGraph): The molecular graph.
            node_id (int): The index of the carbon atom to modify.

        Returns:
            MolecularGraph: The modified molecular graph.
        """
        mod_graph = copy.deepcopy(graph)
        node = next(n for n in mod_graph.nodes if n.index == node_id)

        o_node = next(
            (n for n in mod_graph.nodes if n.element == 'O' and any(e.from_idx == n.index or e.to_idx == n.index
                                                                    for e in mod_graph.edges if e.from_idx == node_id or e.to_idx == node_id)),
            None)

        if o_node:
            h_node = next((n for n in mod_graph.nodes if n.element == 'H' and any(e.from_idx == n.index or e.to_idx == n.index
                                                                                   for e in mod_graph.edges if e.from_idx == o_node.index or e.to_idx == o_node.index)),
                          None)
            if h_node:
                mod_graph.nodes.remove(o_node)
                mod_graph.nodes.remove(h_node)
                mod_graph.edges = [e for e in mod_graph.edges if e.from_idx not in {o_node.index, h_node.index} and e.to_idx not in {o_node.index, h_node.index}]

        return mod_graph

    def _reindex_graph(self, graph: MolecularGraph) -> MolecularGraph:
        """
        Reindex graph nodes and edges sequentially.

        Args:
            graph (MolecularGraph): The molecular graph to reindex.

        Returns:
            MolecularGraph: The reindexed molecular graph.
        """
        new_graph = MolecularGraph()
        old_to_new = {}

        for i, node in enumerate(sorted(graph.nodes, key=lambda x: x.index)):
            new_node = copy.deepcopy(node)
            old_to_new[node.index] = i
            new_node.index = i
            new_graph.nodes.append(new_node)

        for edge in graph.edges:
            new_edge = copy.deepcopy(edge)
            new_edge.from_idx = old_to_new[edge.from_idx]
            new_edge.to_idx = old_to_new[edge.to_idx]
            new_graph.edges.append(new_edge)

        return new_graph
