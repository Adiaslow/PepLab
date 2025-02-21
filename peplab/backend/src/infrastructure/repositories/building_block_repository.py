# peplab/backend/src/infrastructure/repositories/building_block_repository.py

# External imports
from typing import List, Dict

# Internal imports
from peplab.backend.src.domain.models.building_block import BuildingBlock
from peplab.backend.src.core.exceptions.domain_exceptions import (
    InvalidBuildingBlockError,
)
from peplab.backend.src.core.exceptions.infrastructure_exceptions import (
    EmptyBuildingBlockRepositoryError,
)


class BuildingBlockRepository:
    """Repository for building blocks.

    A singleton repository for building blocks.

    Attributes:
        building_blocks: The building blocks in the repository.

    Methods:
        add_building_block: Add a building block to the repository.
        get_building_block: Get a building block from the repository.
        get_all_building_blocks: Get all building blocks from the repository.
    """

    def __init__(self, building_blocks: List[BuildingBlock] = []) -> None:
        """Initialize the building block repository.

        Args:
            building_blocks: The building blocks to add to the repository.
        """
        self.building_blocks: List[BuildingBlock] = building_blocks
        self.building_block_table: Dict[str, BuildingBlock] = {}

    def index_building_blocks(self) -> None:
        """Index the building blocks in the repository."""

    def add_building_block(self, building_block: BuildingBlock) -> None:
        """Add a building block to the repository.

        Args:
            building_block: The building block to add to the repository.
        """
        self.building_blocks.append(building_block)
        self.building_block_table[building_block.id] = building_block

    def get_building_block_by_id(self, id: str) -> BuildingBlock:
        """Get a building block from the repository by id."""
        for building_block in self.building_blocks:
            if building_block.id == id:
                return building_block
        raise InvalidBuildingBlockError(f"Building block with id {id} not found.")

    def get_building_block_by_name(self, name: str) -> BuildingBlock:
        """Get a building block from the repository by name.

        Args:
            name: The name of the building block to get.

        Returns:
            The building block with the given name.

        Raises:
            InvalidBuildingBlockError: If the building block with the given name is not found.
        """
        for building_block in self.building_blocks:
            if building_block.name == name:
                return building_block
        raise InvalidBuildingBlockError(f"Building block with name {name} not found.")

    def get_all_building_blocks(self) -> List[BuildingBlock]:
        """Get all building blocks from the repository.

        Returns:
            The list of building blocks in the repository.

        Raises:
            EmptyBuildingBlockRepositoryError: If the repository is empty.
        """
        if not self.building_blocks:
            raise EmptyBuildingBlockRepositoryError(
                "No building blocks in the repository."
            )
        return self.building_blocks

    def remove_building_block_by_id(self, id: str) -> None:
        """Remove a building block from the repository by id.

        Args:
            id: The id of the building block to remove.

        Raises:
            InvalidBuildingBlockError: If the building block with the given id is not found.
        """
        try:
            self.building_blocks = [
                building_block
                for building_block in self.building_blocks
                if building_block.id != id
            ]
        except Exception as e:
            raise InvalidBuildingBlockError(
                f"Building block with id {id} not found."
            ) from e

    def remove_building_block_by_name(self, name: str) -> None:
        """Remove a building block from the repository by name.

        Args:
            name: The name of the building block to remove.

        Raises:
            InvalidBuildingBlockError: If the building block with the given name is not found.
        """
        try:
            self.building_blocks = [
                building_block
                for building_block in self.building_blocks
                if building_block.name != name
            ]
        except Exception as e:
            raise InvalidBuildingBlockError(
                f"Building block with name {name} not found."
            ) from e
