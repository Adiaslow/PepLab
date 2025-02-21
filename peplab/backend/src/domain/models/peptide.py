# peplab/backend/src/domain/models/peptide.py
"""
This module contains the peptide model.

Classes:
    Peptide: The peptide model.
"""

# Standard Library Imports
from typing import Any, Dict, List, Optional, Union
from uuid import UUID

# Third Party Imports
from pydantic import BaseModel, Field

# Internal imports
from peplab.backend.src.domain.models.building_block import BuildingBlock
from peplab.backend.src.core.exceptions.domain_exceptions import (
    InvalidDirectionError,
    InvalidPeptideError,
    NoBuildingBlocksError,
    InvalidBuildingBlockIndexError,
    NoPropertyOrMetadataError,
)


class Peptide(BaseModel):
    """The peptide model.

    Attributes:
        id: The id of the peptide.
        building_blocks: The building blocks of the peptide.
        properties: The properties of the peptide.
        metadata: The metadata of the peptide.

    Methods:
        __post_init__: Post initialization hook.
        sequence: Return the sequence of the peptide.
        __str__: Return the string representation of the peptide.
        __repr__: Return the representation of the peptide.
        __hash__: Return the hash of the peptide.
        __eq__: Return the equality of the peptide.
        __ne__: Return the inequality of the peptide.
    """

    id: Optional[Union[int, str, UUID]] = Field(
        default=None, description="The id of the peptide."
    )
    name: Optional[str] = Field(default=None, description="The name of the peptide.")
    building_blocks: Optional[List[BuildingBlock]] = Field(
        default=None, description="The building blocks of the peptide."
    )
    embeddings: Optional[Dict[str, Any]] = Field(
        default=None,
        description="The embeddings of the peptide i.e. molecular graph, SMILES, etc.",
    )
    encodings: Optional[Dict[str, Any]] = Field(
        default=None,
        description="The encodings of the peptide i.e. one-hot encoding, fingerprint, etc.",
    )
    properties: Optional[Dict[str, Any]] = Field(
        default_factory=dict, description="The properties of the peptide."
    )
    metadata: Optional[Dict[str, Any]] = Field(
        default_factory=dict, description="The metadata of the peptide."
    )

    def __post_init__(self) -> None:
        """Post initialization hook."""
        if self.id is None:
            self.id = str(uuid.uuid4())

    @property
    def sequence(
        self, direction: Literal["n_to_c", "c_to_n"] = "n_to_c", delimiter: str = "_"
    ) -> Optional[str]:
        """Return the sequence of the peptide.

        Args:
            direction: The direction of the sequence.
            delimiter: The delimiter between the building blocks.

        Returns:
            The sequence of the peptide.

        Raises:
            ValueError: If the direction is invalid.
        """
        if direction == "n_to_c":
            return (
                delimiter.join(
                    [building_block.name for building_block in self.building_blocks]
                )
                if self.building_blocks
                else None
            )
        elif direction == "c_to_n":
            return (
                delimiter.join(
                    [
                        building_block.name
                        for building_block in self.building_blocks[::-1]
                    ]
                )
                if self.building_blocks
                else None
            )
        else:
            raise InvalidDirectionError(
                "Invalid direction. Must be 'n_to_c' or 'c_to_n'."
            )

    def __str__(self) -> str:
        """Return the string representation of the peptide.

        Returns:
            The string representation of the peptide.
        """
        return f"Peptide(sequence={self.sequence}, id={self.id})"

    def __repr__(self) -> str:
        """Return the representation of the peptide.

        Returns:
            The representation of the peptide.
        """
        return f"Peptide(sequence={self.sequence}, id={self.id})"

    def __hash__(self) -> int:
        """Return the hash of the peptide.

        Returns:
            The hash of the peptide.
        """
        return hash(self.id)

    def __eq__(self, other: Any) -> bool:
        """Return the equality of the peptide.

        Returns:
            The equality of the peptide.
        """
        return self.sequence == other.sequence

    def __ne__(self, other: Any) -> bool:
        """Return the inequality of the peptide.

        Returns:
            The inequality of the peptide.
        """
        return not self.__eq__(other)

    def __len__(self) -> int:
        """Return the length of the peptide.

        Returns:
            The length of the peptide.
        """
        if self.building_blocks is None:
            raise NoBuildingBlocksError("No building blocks set for the peptide.")
        return len(self.building_blocks)

    def __getitem__(self, value: Union[int, str]) -> Union[BuildingBlock, Any]:
        """Return the building block at the given index or get a property or metadata value.

        Args:
            value: The index of the building block or a property or metadata key.

        Returns:
            The building block at the given index or the value of the property or metadata.

        Raises:
            NoBuildingBlocksError: If no building blocks are set for the peptide.
            InvalidBuildingBlockIndexError: If the index is invalid.
            NoPropertyOrMetadataError: If the property or metadata key is invalid.
        """
        if isinstance(value, int):
            if self.building_blocks is None:
                raise NoBuildingBlocksError("No building blocks set for the peptide.")
            if value < 0 or value >= len(self.building_blocks):
                raise InvalidBuildingBlockIndexError(
                    f"Invalid index. Must be between 0 and {len(self.building_blocks)}."
                )
            return self.building_blocks[value]
        elif isinstance(value, str) and (
            self.properties is not None or self.metadata is not None
        ):
            if self.properties is not None and value in self.properties.keys():
                return self.properties[value]
            if self.metadata is not None and value in self.metadata.keys():
                return self.metadata[value]

            raise NoPropertyOrMetadataError(
                "Invalid property or metadata key. Must be one of "
                + f"{self.properties.keys() if self.properties is not None else []} "
                + f"or {self.metadata.keys() if self.metadata is not None else []}."
            )

    def __setitem__(self, index: int, value: BuildingBlock) -> None:
        """Set the building block at the given index.

        Args:
            index: The index of the building block to set.
            value: The building block to set.

        Raises:
            InvalidBuildingBlockIndexError: If the index is invalid.
        """
        if index is 0 and self.building_blocks is None:
            self.building_blocks = [value]
        elif self.building_blocks is not None and (
            index > len(self.building_blocks) + 1 or index < 0
        ):
            raise InvalidBuildingBlockIndexError(
                f"Invalid index. Must be between 0 and {len(self.building_blocks)}."
            )
        elif self.building_blocks is not None:
            self.building_blocks[index] = value
        else:
            raise InvalidBuildingBlockIndexError(
                "Invalid index. Must be between 0 and "
                + f"{len(self.building_blocks)+1 if self.building_blocks is not None else 1}."
            )

    def __delitem__(self, index: int) -> None:
        """Delete the building block at the given index.

        Args:
            index: The index of the building block to delete.

        Raises:
            NoBuildingBlocksError: If no building blocks are set for the peptide.
            InvalidBuildingBlockIndexError: If the index is invalid.
        """
        if self.building_blocks is None:
            raise NoBuildingBlocksError("No building blocks set for the peptide.")
        if index < 0 or index >= len(self.building_blocks):
            raise InvalidBuildingBlockIndexError(
                f"Invalid index. Must be between 0 and {len(self.building_blocks)}."
            )
        del self.building_blocks[index]

    def __contains__(self, item: BuildingBlock) -> bool:
        """Return True if the building block is in the peptide.

        Returns:
            True if the building block is in the peptide.

        Raises:
            NoBuildingBlocksError: If no building blocks are set for the peptide.
        """
        if self.building_blocks is None:
            raise NoBuildingBlocksError("No building blocks set for the peptide.")
        return item in self.building_blocks

    def __add__(self, other: Union[BuildingBlock, "Peptide"]) -> "Peptide":
        """Add a building block or peptide to the peptide.

        Returns:
            The peptide with the added building block.
        """
        if self.building_blocks is None and isinstance(other, BuildingBlock):
            self.building_blocks = [other]
            return self
        if self.building_blocks is not None and isinstance(other, BuildingBlock):
            self.building_blocks.append(other)
            return self
        if (
            self.building_blocks is None
            and isinstance(other, Peptide)
            and other.building_blocks is not None
        ):
            self.building_blocks = other.building_blocks
            return self
        if (
            self.building_blocks is not None
            and isinstance(other, Peptide)
            and other.building_blocks is not None
        ):
            self.building_blocks += other.building_blocks
            return self
        raise InvalidPeptideError("Does not form a valid peptide.")

    def __iadd__(self, other: Union[BuildingBlock, "Peptide"]) -> "Peptide":
        """Append a building block or peptide to the peptide.

        Keeps the original peptide object.

        Returns:
            The peptide with the appended building block.
        """
        if self.building_blocks is None and isinstance(other, BuildingBlock):
            self.building_blocks = [other]
            return self
        if self.building_blocks is not None and isinstance(other, BuildingBlock):
            self.building_blocks.append(other)
            return self
        if (
            self.building_blocks is None
            and isinstance(other, Peptide)
            and other.building_blocks is not None
        ):
            self.building_blocks = other.building_blocks
            return self
        if (
            self.building_blocks is not None
            and isinstance(other, Peptide)
            and other.building_blocks is not None
        ):
            self.building_blocks += other.building_blocks
            return self
        raise InvalidPeptideError("Does not form a valid peptide.")
