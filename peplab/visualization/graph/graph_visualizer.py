from typing import Tuple, Optional
import matplotlib.pyplot as plt
import networkx as nx

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
        print("Creating graph plot...")
        try:
            from peplab.visualization.cpk_colors import CPKColors
            COLORS = CPKColors.COLORS
        except ImportError:
            COLORS = {
                'H': '#FFFFFF',   # White
                'He': '#D9FFFF',  # Light cyan
                'Li': '#CC80FF',  # Light purple
                'Be': '#C2FF00',  # Lime green
                'B': '#FFB5B5',   # Pink
                'C': '#909090',   # Grey
                'N': '#3050F8',   # Blue
                'O': '#FF0D0D',   # Red
                'F': '#90E050',   # Light green
                'Ne': '#B3E3F5',  # Cyan
                'Na': '#AB5CF2',  # Purple
                'Mg': '#8AFF00',  # Bright green
                'Al': '#BFA6A6',  # Light brown
                'Si': '#F0C8A0',  # Peach
                'P': '#FF8000',   # Orange
                'S': '#FFFF30',   # Yellow
                'Cl': '#1FF01F',  # Green
                'Ar': '#80D1E3',  # Light blue
                'K': '#8F40D4',   # Purple
                'Ca': '#3DFF00',  # Bright green
                'Sc': '#E6E6E6',  # Light grey
                'Ti': '#BFC2C7',  # Grey
                'V': '#A6A6AB',   # Grey
                'Cr': '#8A99C7',  # Grey-blue
                'Mn': '#9C7AC7',  # Purple-grey
                'Fe': '#E06633',  # Orange-brown
                'Co': '#F090A0',  # Pink
                'Ni': '#50D050',  # Green
                'Cu': '#C88033',  # Brown
                'Zn': '#7D80B0',  # Grey-blue
                'Ga': '#C28F8F',  # Pink-brown
                'Ge': '#668F8F',  # Grey-green
                'As': '#BD80E3',  # Purple
                'Se': '#FFA100',  # Orange
                'Br': '#A62929',  # Dark red
                'Kr': '#5CB8D1',  # Light blue
                'Rb': '#702EB0',  # Dark purple
                'Sr': '#00FF00',  # Green
                'Y': '#94FFFF',   # Light cyan
                'Zr': '#94E0E0',  # Light blue
                'Nb': '#73C2C9',  # Light blue
                'Mo': '#54B5B5',  # Blue-green
                'Tc': '#3B9E9E',  # Blue-green
                'Ru': '#248F8F',  # Blue-green
                'Rh': '#0A7D8C',  # Dark blue-green
                'Pd': '#006985',  # Dark blue
                'Ag': '#C0C0C0',  # Silver
                'Cd': '#FFD98F',  # Light yellow
                'In': '#A67573',  # Brown
                'Sn': '#668080',  # Grey-blue
                'Sb': '#9E63B5',  # Purple
                'Te': '#D47A00',  # Brown
                'I': '#940094',   # Purple
                'Xe': '#429EB0',  # Blue
                'Cs': '#57178F',  # Dark purple
                'Ba': '#00C900',  # Green
                'La': '#70D4FF',  # Light blue
                'Ce': '#FFFFC7',  # Very light yellow
                'Pr': '#D9FFC7',  # Very light green
                'Nd': '#C7FFC7',  # Light green
                'Pm': '#A3FFC7',  # Light green
                'Sm': '#8FFFC7',  # Light green
                'Eu': '#61FFC7',  # Light green
                'Gd': '#45FFC7',  # Light green
                'Tb': '#30FFC7',  # Light green
                'Dy': '#1FFFC7',  # Light green
                'Ho': '#00FF9C',  # Green
                'Er': '#00E675',  # Green
                'Tm': '#00D452',  # Green
                'Yb': '#00BF38',  # Green
                'Lu': '#00AB24',  # Green
                'Hf': '#4DC2FF',  # Light blue
                'Ta': '#4DA6FF',  # Light blue
                'W': '#2194D6',   # Blue
                'Re': '#267DAB',  # Blue
                'Os': '#266696',  # Blue
                'Ir': '#175487',  # Dark blue
                'Pt': '#D0D0E0',  # Grey
                'Au': '#FFD123',  # Gold
                'Hg': '#B8B8D0',  # Grey
                'Tl': '#A6544D',  # Brown
                'Pb': '#575961',  # Dark grey
                'Bi': '#9E4FB5',  # Purple
                'Po': '#AB5C00',  # Brown
                'At': '#754F45',  # Dark brown
                'Rn': '#428296',  # Blue
                'Fr': '#420066',  # Dark purple
                'Ra': '#007D00',  # Dark green
                'Ac': '#70ABFA',  # Light blue
                'Th': '#00BAFF',  # Light blue
                'Pa': '#00A1FF',  # Light blue
                'U': '#008FFF',   # Blue
                'Np': '#0080FF',  # Blue
                'Pu': '#006BFF',  # Blue
                'Am': '#545CF2',  # Blue
                'Cm': '#785CE3',  # Blue
                'Bk': '#8A4FE3',  # Purple
                'Cf': '#A136D4',  # Purple
                'Es': '#B31FD4',  # Purple
                'Fm': '#B31FBA',  # Purple
                'Md': '#B30DA6',  # Purple
                'No': '#BD0D87',  # Purple
                'Lr': '#C70066',  # Dark pink
                'Rf': '#CC0059',  # Dark pink
                'Db': '#D1004F',  # Dark pink
                'Sg': '#D90045',  # Dark pink
                'Bh': '#E00038',  # Dark pink
                'Hs': '#E6002E',  # Dark pink
                'Mt': '#EB0026',  # Dark pink
            }

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
                    color=COLORS[node['element']],
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
