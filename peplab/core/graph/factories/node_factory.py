# peplab/core/graph/factories/node_factory.py

"""
Factory module for creating nodes.
"""

from ..leaves import Node

class NodeFactory:
    """
    Factory class for creating nodes.
    """

    @staticmethod
    def create_node(index, element, properties) -> Node:
        """
        Creates a Node instance.

        Args:
            index (int): The index of the node.
            element (str): The chemical element of the node.
            properties (dict): Additional properties of the node.

        Returns:
            Node: An instance of Node.
        """
        return Node(index=index, element=element, properties=properties)
