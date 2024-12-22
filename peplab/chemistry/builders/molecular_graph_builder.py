# peplab/chemistry/builders/molecular_graph_builder.py

"""
Module for building molecular graphs.

This module provides the MolecularGraphBuilder class which serves as a builder for creating MolecularGraph instances
from SMILES strings.
"""

from rdkit import Chem
from ...core.graph.leaves import Node, Edge
from ...core.graph.composites import MolecularGraph
from ...analysis.properties.molecular_graph_property_calculator import MolecularGraphPropertyCalculator

class MolecularGraphBuilder:
    """
    Builder class for creating MolecularGraph instances from SMILES strings.
    """

    def __init__(self):
        self.graph = MolecularGraph()

    def from_smiles(self, smiles: str) -> MolecularGraph:
        """
        Constructs a molecular graph from a SMILES string.

        Args:
            smiles (str): The SMILES string representing the molecule.

        Returns:
            MolecularGraph: An instance of MolecularGraph.

        Raises:
            ValueError: If the SMILES string cannot be parsed or is empty.
        """
        if not smiles:
            raise ValueError("SMILES string is empty")

        try:
            self.graph.logger.debug(f"Parsing SMILES: {smiles}")
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError(f"Failed to parse SMILES: {smiles}")

            Chem.SanitizeMol(mol)
            mol = Chem.AddHs(mol)

            # Add atoms
            for atom in mol.GetAtoms():
                idx = atom.GetIdx()
                node = Node(
                    index=idx,
                    element=atom.GetSymbol(),
                    properties=MolecularGraphPropertyCalculator.get_atom_properties(atom)
                )
                self.graph.add_node(node)

            # Add bonds
            for bond in mol.GetBonds():
                edge = Edge(
                    from_idx=bond.GetBeginAtomIdx(),
                    to_idx=bond.GetEndAtomIdx(),
                    properties=MolecularGraphPropertyCalculator.get_bond_properties(bond)
                )
                self.graph.add_edge(edge)

            return self.graph

        except Exception as e:
            self.graph.logger.error(f"Error creating molecular graph from SMILES: {str(e)}")
            raise
