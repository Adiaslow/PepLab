# peplab/backend/src/domain/models/library.py
"""
The library model.

Classes:
    Library: The library model.

Todo:
    Fix Library lookup tables.
    Fix Library.__getitem__ method.
    Fix Library__index_peptides method.

"""
# Standard Library Imports
from typing import Dict, List, Union, Tuple
from uuid import UUID

# Third Party Imports
from pydantic import BaseModel, Field

# Internal Imports
from peplab.backend.src.domain.models.building_block import BuildingBlock
from peplab.backend.src.domain.models.peptide import Peptide
from peplab.backend.src.core.exceptions.domain_exceptions import (
    InvalidBuildingBlocksError,
    InvalidLibraryIndexError,
    InvalidPeptideNameOrIDError,
    InvalidPeptideUUIDError,
)


class Library(BaseModel):
    """A library of peptides.

    Attributes:
        peptides: The peptides in the library.
        peptides_by_id_lookup: A lookup table for the peptides in the library.
        peptides_by_name_lookup: A lookup table for the peptides in the library by name.
        peptides_by_building_blocks_lookup: A lookup table for the peptides in the library by building blocks.

    Methods:
        __len__: Return the number of peptides in the library.
        __getitem__: Return the peptide at the given index.

    Todo:
        Fix the __getitem__ method.
    """

    peptides: List[Peptide] = Field(..., description="The peptides in the library.")
    peptides_by_id_lookup: Dict[int, Peptide] = Field(
        default_factory=dict,
        description="A lookup table for the peptides in the library.",
    )
    peptides_by_name_lookup: Dict[str, Peptide] = Field(
        default_factory=dict,
        description="A lookup table for the peptides in the library by name.",
    )
    peptides_by_building_blocks_lookup: Dict[Tuple[BuildingBlock, ...], Peptide] = (
        Field(
            default_factory=dict,
            description="A lookup table for the peptides in the library by building blocks.",
        )
    )

    def __len__(self) -> int:
        """Return the number of peptides in the library.

        Returns:
            The number of peptides in the library.
        """
        return len(self.peptides)

    # FIXME
    def __getitem__(self, key: Union[int, str, UUID, List[BuildingBlock]]) -> Peptide:
        """Return the peptide at the given index.

        Args:
            key: The key to get the peptide from.

        Returns:
            The peptide at the given index.
        """
        if isinstance(key, int):
            if key < 0 or key >= len(self.peptides):
                raise InvalidLibraryIndexError("Index out of bounds")
            return self.peptides[key]
        if isinstance(key, str):
            if (
                key not in self.peptides_by_name_lookup
                and key not in self.peptides_by_id_lookup
            ):
                raise InvalidPeptideNameOrIDError(
                    f"Peptide with name or ID {key} not found in library."
                )
            try:
                return self.peptides_by_name_lookup[key]
            except:
                try:
                    return self.peptides_by_id_lookup[int(key)]
                except:
                    raise InvalidPeptideNameOrIDError(
                        f"Peptide with name or ID {key} not found in library."
                    )
        if isinstance(key, UUID):
            try:
                return self.peptides_by_id_lookup[int(key)]
            except:
                raise InvalidPeptideUUIDError(
                    f"Peptide with UUID {key} not found in library."
                )
        if isinstance(key, List):
            if len(key) == 0:
                raise InvalidBuildingBlocksError("No building blocks provided.")
            try:
                return self.peptides_by_building_blocks_lookup[tuple(key)]
            except KeyError:
                raise InvalidBuildingBlocksError(
                    f"Peptide with building blocks {key} not found in library."
                )

    def __contains__(self, key: Union[int, str, UUID, List[BuildingBlock]]) -> bool:
        """Check if a peptide is in the library.

        Args:
            key: The key to check if the peptide is in the library.
        """
        if isinstance(key, int):
            return key in self.peptides_by_id_lookup
        if isinstance(key, str):
            return key in self.peptides_by_name_lookup
        if isinstance(key, UUID):
            return key in self.peptides_by_id_lookup
        if isinstance(key, List):
            return key in self.peptides_by_building_blocks_lookup
        return False

    def __add__(self, peptide: Peptide) -> "Library":
        """Add a peptide to the library.

        Args:
            peptide: The peptide to add to the library.
        """
        self.peptides.append(peptide)
        return self

    def __iadd__(self, peptide: Peptide) -> "Library":
        """Add a peptide to the library.

        Args:
            peptide: The peptide to add to the library.
        """
        return self + peptide

    def __sub__(self, peptide: Peptide) -> "Library":
        """Remove a peptide from the library.

        Args:
            peptide: The peptide to remove from the library.
        """
        self.peptides.remove(peptide)
        return self

    def __isub__(self, peptide: Peptide) -> "Library":
        """Remove a peptide from the library.

        Args:
            peptide: The peptide to remove from the library.
        """
        return self - peptide

    # FIXME
    def __index_peptides(self) -> None:
        """Index the peptides in the library lookup tables.

        Returns:
            None.
        """
        self.peptides_by_id_lookup = {}
        self.peptides_by_name_lookup = {}
        self.peptides_by_building_blocks_lookup = {}
        for peptide in self.peptides:
            self.peptides_by_id_lookup[peptide.id] = peptide  # type: ignore
            self.peptides_by_name_lookup[peptide.name] = peptide  # type: ignore
            self.peptides_by_building_blocks_lookup[tuple(peptide.building_blocks)] = (  # type: ignore
                peptide
            )

    def count_building_blocks(self) -> int:
        """Count the number of building blocks in the library.

        Returns:
            The number of building blocks in the library.
        """
        return sum(
            len(peptide.building_blocks)
            for peptide in self.peptides
            if peptide.building_blocks is not None
        )

    def count_building_block(self, building_block: BuildingBlock) -> int:
        """Count the number of times a building block appears in the library.

        Args:
            building_block: The building block to count.
        """
        return sum(
            building_block in peptide.building_blocks
            for peptide in self.peptides
            if peptide.building_blocks is not None
        )

    def count_building_blocks_at_position(self, position: int) -> int:
        """Count the number of building blocks at a given position in the library.

        Args:
            position: The position to count the building blocks at.
        """
        return sum(
            len(peptide.building_blocks)
            for peptide in self.peptides
            if peptide.building_blocks is not None
            and peptide.building_blocks[position] is not None
        )

    def count_building_block_at_position(
        self, building_block: BuildingBlock, position: int
    ) -> int:
        """Count the number of times a building block appears at a given position in the library.

        Args:
            building_block: The building block to count.
            position: The position to count the building block at.
        """
        return sum(
            building_block == peptide.building_blocks[position]
            for peptide in self.peptides
            if peptide.building_blocks is not None
        )

    def count_building_blocks_at_positions(self, positions: List[int]) -> int:
        """Count the number of building blocks at a given positions in the library.

        Args:
            positions: The positions to count the building blocks at.
        """
        return sum(
            self.count_building_blocks_at_position(position) for position in positions
        )

    def count_building_block_at_positions(
        self, building_block: BuildingBlock, positions: List[int]
    ) -> int:
        """Count the number of times a building block appears at a given positions in the library.

        Args:
            building_block: The building block to count.
            positions: The positions to count the building block at.
        """
        return sum(
            self.count_building_block_at_position(building_block, position)
            for position in positions
        )
