# peplab/core/graph/factories/edge_factory.py

"""
Factory module for creating edges.
"""

from ..leaves import Edge

class EdgeFactory:
    """
    Factory class for creating edges.
    """

    @staticmethod
    def create_edge(from_idx, to_idx, properties) -> Edge:
        """
        Creates an Edge instance.

        Args:
            from_idx (int): The starting node index of the edge.
            to_idx (int): The ending node index of the edge.
            properties (dict): Additional properties of the edge.

        Returns:
            Edge: An instance of Edge.
        """
        return Edge(from_idx=from_idx, to_idx=to_idx, properties=properties)
