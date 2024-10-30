from typing import Tuple, Optional
import matplotlib.pyplot as plt
import networkx as nx

from peplab.visualization.cpk_colors import CPKColors

class GraphVisualizer:
    """Visualizes molecular graphs using NetworkX."""

    @staticmethod
    def create_graph_plot(
        graph_dict: dict,
        size: Tuple[int, int] = (10, 10),
        title: Optional[str] = None    ) -> plt.Figure:
        """Creates a plot of the molecular graph.

        Args:
            graph_dict: Dictionary containing graph data.
            size: Figure size in inches.
            title: Optional title for the plot.

        Returns:
            Matplotlib figure object.
        """
        # Create NetworkX graph
        G = nx.Graph()

        # Add nodes
        for node in graph_dict['nodes']:
            if node['element'] != 'H':  # Skip explicit hydrogens
                label = f"{node['element']}{node['id']}"
                if node['total_num_hs'] > 0:
                    label += f"H{node['total_num_hs']}"
                if node['formal_charge'] != 0:
                    label += f"{'+' if node['formal_charge'] > 0 else ''}{node['formal_charge']}"

                G.add_node(
                    node['id'],
                    label=label,
                    color=CPKColors.get_color(node['element']),
                    is_reactive_nuc=node.get('is_reactive_nuc', False),
                    is_reactive_elec=node.get('is_reactive_elec', False)
                )

        # Add edges
        for edge in graph_dict['edges']:
            if (edge['from_idx'] in G.nodes and
                edge['to_idx'] in G.nodes):
                style = 'solid'
                width = 1
                if 'DOUBLE' in edge['bond_type']:
                    width = 2
                elif 'TRIPLE' in edge['bond_type']:
                    width = 3
                elif 'AROMATIC' in edge['bond_type']:
                    style = 'dashed'

                G.add_edge(
                    edge['from_idx'],
                    edge['to_idx'],
                    style=style,
                    width=width
                )

        # Create figure
        fig = plt.figure(figsize=size)

        # Create layout
        pos = nx.spring_layout(G, k=0.1, iterations=1000, threshold=1e-10)

        # Draw different node types
        GraphVisualizer._draw_nodes(G, pos)

        # Draw edges
        GraphVisualizer._draw_edges(G, pos)

        # Draw labels
        nx.draw_networkx_labels(
            G, pos,
            labels=nx.get_node_attributes(G, 'label'),
            font_size=8
        )

        # Add legend
        GraphVisualizer._add_legend()

        if title:
            plt.title(title)
        plt.axis('off')
        plt.tight_layout()

        return fig

    @staticmethod
    def _draw_nodes(G: nx.Graph, pos: dict) -> None:
        """Draws nodes with appropriate styles."""
        # Regular nodes
        regular_nodes = [
            n for n, attr in G.nodes(data=True)
            if not (attr['is_reactive_nuc'] or attr['is_reactive_elec'])
        ]
        if regular_nodes:
            nx.draw_networkx_nodes(
                G, pos,
                nodelist=regular_nodes,
                node_size=500,
                node_color=[G.nodes[n]['color'] for n in regular_nodes],
                edgecolors='black',
                linewidths=1.0
            )

        # Nucleophilic sites
        nuc_nodes = [
            n for n, attr in G.nodes(data=True)
            if attr['is_reactive_nuc']
        ]
        if nuc_nodes:
            nx.draw_networkx_nodes(
                G, pos,
                nodelist=nuc_nodes,
                node_size=500,
                node_color=[G.nodes[n]['color'] for n in nuc_nodes],
                edgecolors='blue',
                linewidths=3.0
            )

        # Electrophilic sites
        elec_nodes = [
            n for n, attr in G.nodes(data=True)
            if attr['is_reactive_elec']
        ]
        if elec_nodes:
            nodes = nx.draw_networkx_nodes(
                G, pos,
                nodelist=elec_nodes,
                node_size=500,
                node_color=[G.nodes[n]['color'] for n in elec_nodes],
                edgecolors='none'
            )
            # Add dashed border
            nodes.set_edgecolor('red')
            nodes.set_linewidth(3.0)
            nodes.set_linestyle('--')

    @staticmethod
    def _draw_edges(G: nx.Graph, pos: dict) -> None:
        """Draws edges with appropriate styles."""
        edge_styles = nx.get_edge_attributes(G, 'style')
        edge_widths = nx.get_edge_attributes(G, 'width')

        # Draw solid edges
        solid_edges = [
            (u, v) for (u, v), style in edge_styles.items()
            if style == 'solid'
        ]
        if solid_edges:
            nx.draw_networkx_edges(
                G, pos,
                edgelist=solid_edges,
                width=[edge_widths[(u, v)] for (u, v) in solid_edges],
                edge_color='black'
            )

        # Draw dashed edges (aromatic)
        dashed_edges = [
            (u, v) for (u, v), style in edge_styles.items()
            if style == 'dashed'
        ]
        if dashed_edges:
            nx.draw_networkx_edges(
                G, pos,
                edgelist=dashed_edges,
                style='dashed',
                width=1,
                edge_color='black'
            )

    @staticmethod
    def _add_legend() -> None:
        """Adds a legend to the plot."""
        legend_elements = [
            plt.Line2D(
                [0], [0],
                marker='o',
                color='w',
                markerfacecolor='gray',
                markeredgecolor='blue',
                markeredgewidth=3,
                markersize=15,
                label='Nucleophilic Site'
            ),
            plt.Line2D(
                [0], [0],
                marker='o',
                color='w',
                markerfacecolor='gray',
                markeredgecolor='red',
                markeredgewidth=3,
                markersize=15,
                label='Electrophilic Site',
                dashes=[2, 2]
            )
        ]
        plt.legend(handles=legend_elements, loc='upper right')
