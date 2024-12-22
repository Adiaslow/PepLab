# tests/test_reactive_group_identifier.py

import pytest
from peplab.chemistry.reactivity.reactive_group_identifier import ReactiveGroupIdentifier
from peplab.chemistry.builders.molecular_graph_builder import MolecularGraphBuilder
from peplab.utils.logging import Logger

def test_identify_reactive_groups():
    builder = MolecularGraphBuilder()
    logger = Logger().get_logger()

    # Define SMILES strings for different molecules
    benzene_smiles = "c1ccccc1"
    ethyne_smiles = "C#C"
    methane_smiles = "C"
    primary_amine_smiles = "CN"
    hydroxyl_smiles = "CO"

    # Build molecular graphs from SMILES strings
    benzene_graph = builder.from_smiles(benzene_smiles)
    ethyne_graph = builder.from_smiles(ethyne_smiles)
    methane_graph = builder.from_smiles(methane_smiles)
    primary_amine_graph = builder.from_smiles(primary_amine_smiles)
    hydroxyl_graph = builder.from_smiles(hydroxyl_smiles)

    # Create a list of graphs for easier iteration
    graphs = [
        benzene_graph,
        ethyne_graph,
        methane_graph,
        primary_amine_graph,
        hydroxyl_graph
    ]

    # Expected results
    expected_counts = {
        'aromatic_carbon': 6,
        'alkyne': 1,
        'methyl': 1,
        'primary_amine': 1,
        'hydroxyl': 1
    }

    for graph in graphs:
        identifier = ReactiveGroupIdentifier(graph)
        reactive_groups = identifier.identify_reactive_groups()

        aromatic_carbons = [group for group in reactive_groups if group.name == 'aromatic_carbon']
        alkynes = [group for group in reactive_groups if group.name == 'alkyne']
        methyls = [group for group in reactive_groups if group.name == 'methyl']
        primary_amines = [group for group in reactive_groups if group.name == 'primary_amine']
        hydroxyls = [group for group in reactive_groups if group.name == 'hydroxyl']

        logger.debug(f"Aromatic carbons: {len(aromatic_carbons)}")
        logger.debug(f"Alkynes: {len(alkynes)}")
        logger.debug(f"Methyls: {len(methyls)}")
        logger.debug(f"Primary amines: {len(primary_amines)}")
        logger.debug(f"Hydroxyls: {len(hydroxyls)}")

        assert len(aromatic_carbons) == expected_counts['aromatic_carbon'], "All aromatic carbons should be identified correctly."
        assert len(alkynes) == expected_counts['alkyne'], "The alkyne should be identified correctly."
        assert len(methyls) == expected_counts['methyl'], "The methyl group should be identified correctly."
        assert len(primary_amines) == expected_counts['primary_amine'], "The primary amine should be identified correctly."
        assert len(hydroxyls) == expected_counts['hydroxyl'], "The hydroxyl group should be identified correctly."

if __name__ == "__main__":
    pytest.main()
