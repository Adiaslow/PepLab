# peplab/frontend/src/core/types/initialization_types.py
"""
This module contains the types for the initialization state of the application.
"""

# Standard Library Imports
from enum import Enum
from typing import List


class InitializationPhase(Enum):
    """The phase of the initialization state."""

    START = "start"
    CHECKING = "checking"
    COMPLETE = "complete"
    FAILED = "failed"

    @classmethod
    def get_all_phases(cls) -> List[str]:
        """Get all phases."""
        return [phase.value for phase in cls]
