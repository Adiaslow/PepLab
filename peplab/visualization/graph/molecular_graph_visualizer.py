# peplab/visualization/graph/molecular_graph_visualizer.py

import networkx as nx
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import rdDepictor
from peplab.chemistry.reactivity.reactive_group_identifier import ReactiveGroupIdentifier
from peplab.core.graph.composites import MolecularGraph
from peplab.core.graph.leaves import Node, Edge

class MolecularGraphVisualizer:
    def __init__(self, smiles: str):
        self.smiles = smiles
        self.mol = Chem.MolFromSmiles(smiles)
        if not self.mol:
            raise ValueError(f"Invalid SMILES string: {smiles}")
        rdDepictor.Compute2DCoords(self.mol)
        self.graph = self._convert_to_molecular_graph()
        self.nx_graph = self._create_networkx_graph()

    def _convert_to_molecular_graph(self) -> MolecularGraph:
        graph = MolecularGraph()

        # Add nodes (atoms)
        for atom in self.mol.GetAtoms():
            node = Node(
                index=atom.GetIdx(),
                element=atom.GetSymbol(),
                properties={
                    'num_implicit_hs': atom.GetNumImplicitHs(),
                    'aromatic': atom.GetIsAromatic(),
                    'charge': atom.GetFormalCharge()
                }
            )
            graph.add_node(node)

        # Add edges (bonds)
        for bond in self.mol.GetBonds():
            edge = Edge(
                from_idx=bond.GetBeginAtomIdx(),
                to_idx=bond.GetEndAtomIdx(),
                properties={'bond_type': bond.GetBondType().name}
            )
            graph.add_edge(edge)

        return graph

    def _create_networkx_graph(self) -> nx.Graph:
        G = nx.Graph()
        for node in self.graph.nodes:
            G.add_node(node.index, element=node.element, **node.properties)
        for edge in self.graph.edges:
            G.add_edge(edge.from_idx, edge.to_idx, bond_type=edge.properties['bond_type'])
        return G

    def show_molecule(self, file_path: str):
        """
        Display the molecule and save it to a file.
        """
        pos = nx.spring_layout(self.nx_graph)
        labels = {node: data['element'] for node, data in self.nx_graph.nodes(data=True)}
        nx.draw(self.nx_graph, pos, with_labels=True, labels=labels, node_color='skyblue', node_size=500, font_size=10)
        plt.savefig(file_path)
        plt.show()
        plt.close()

    def highlight_reactive_groups(self, file_path: str):
        """
        Highlight reactive groups in the molecule using the existing code for finding reactive groups.
        Save the highlighted molecule to a file. Unknown reactive groups are skipped.
        """
        identifier = ReactiveGroupIdentifier(self.graph)
        reactive_groups = identifier.identify_reactive_groups()
        highlight_colors = {}
        color_palette = ['red', 'green', 'blue', 'yellow', 'magenta', 'cyan']
        legend_labels = []

        # Filter out None values and then process valid groups
        valid_groups = [group for group in reactive_groups if group is not None]

        for i, group in enumerate(valid_groups):
            color = color_palette[i % len(color_palette)]
            legend_labels.append((color, group.name))
            for node in group.nodes:
                highlight_colors[node.index] = color

        pos = nx.spring_layout(self.nx_graph, seed=137, k=0.5, iterations=100)
        labels = {node: data['element'] for node, data in self.nx_graph.nodes(data=True)}
        node_colors = [highlight_colors.get(node, 'skyblue') for node in self.nx_graph.nodes()]
        nx.draw(self.nx_graph, pos, with_labels=True, labels=labels,
                node_color=node_colors, node_size=500, font_size=10)

        # Add legend
        for color, label in legend_labels:
            plt.plot([], [], color=color, label=label)
        plt.legend(loc='upper right')

        plt.savefig(file_path)
        plt.show()
        plt.close()
