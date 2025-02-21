# peplab/backend/src/core/exceptions/__init__.py
"""
The core exceptions module.
"""

# Standard Library Imports
from typing import List

# Internal Imports
from peplab.backend.src.core.exceptions.domain_exceptions import (
    InvalidBuildingBlockError,
    InvalidBuildingBlockIndexError,
    InvalidPeptideError,
    InvalidDirectionError,
    NoBuildingBlocksError,
    NoPropertyOrMetadataError,
    InvalidLibraryError,
    InvalidLibraryIndexError,
    InvalidPeptideNameOrIDError,
    InvalidBuildingBlocksError,
    InvalidPeptideUUIDError,
)

from peplab.backend.src.core.exceptions.infrastructure_exceptions import (
    InfrastructureException,
    EmptyBuildingBlockRepositoryError,
)

__all__: List[str] = [
    "InvalidBuildingBlockError",
    "InvalidBuildingBlockIndexError",
    "InvalidPeptideError",
    "InvalidDirectionError",
    "NoBuildingBlocksError",
    "NoPropertyOrMetadataError",
    "InvalidLibraryError",
    "InvalidLibraryIndexError",
    "InvalidPeptideNameOrIDError",
    "InvalidBuildingBlocksError",
    "InvalidPeptideUUIDError",
    "InfrastructureException",
    "EmptyBuildingBlockRepositoryError",
]
