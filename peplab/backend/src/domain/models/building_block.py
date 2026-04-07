# peplab/backend/src/domain/models/building_block.py
"""
This module contains the building block model.

Classes:
    BuildingBlock: The building block model.
"""

# Standard Library Imports
from typing import Any, Dict, Optional, Union
from uuid import UUID, uuid4

# Third Party Imports
from pydantic import BaseModel, Field

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from peplab.backend.src.domain.models.peptide import Peptide


class BuildingBlock(BaseModel):
    """The building block model."""

    id: Optional[Union[int, str, UUID]] = Field(
        default=None, description="The id of the building block."
    )
    name: str = Field(..., description="The name of the building block.")
    embeddings: Optional[Dict[str, Any]] = Field(
        default=None,
        description="The embeddings of the building block i.e. molecular graph, SMILES, etc.",
    )
    encodings: Optional[Dict[str, Any]] = Field(
        default=None,
        description="The encodings of the building block i.e. one-hot encoding, fingerprint, etc.",
    )
    properties: Optional[Dict[str, Any]] = Field(
        default_factory=dict, description="The properties of the building block."
    )
    metadata: Optional[Dict[str, Any]] = Field(
        default_factory=dict, description="The metadata of the building block."
    )

    def __post_init__(self) -> None:
        """Post initialization hook."""
        if self.id is None:
            self.id = str(uuid4())

    def __hash__(self) -> int:
        """Hash the building block."""
        return hash(self.id)

    def __eq__(self, other: Any) -> bool:
        """Check if the building block is equal to another object."""
        if not isinstance(other, BuildingBlock):
            return False
        return self.name == other.name

    def __ne__(self, other: Any) -> bool:
        """Check if the building block is not equal to another object."""
        if not isinstance(other, BuildingBlock):
            return True
        return self.name != other.name

    def __str__(self) -> str:
        """Return the string representation of the building block."""
        return self.name

    def __repr__(self) -> str:
        """Return the representation of the building block."""
        return self.name

    def __add__(self, other: "BuildingBlock") -> "Peptide":
        """Add two building blocks together to get a Peptide."""
        from peplab.backend.src.domain.models.peptide import Peptide
        return Peptide(building_blocks=[self, other])

    def __iadd__(self, other: "BuildingBlock") -> "Peptide":
        """Add two building blocks together to get a Peptide."""
        from peplab.backend.src.domain.models.peptide import Peptide
        return Peptide(building_blocks=[self, other])
