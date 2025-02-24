# peplab/frontend/src/application/states/analysis_state.py
"""
This module handles the analysis state of the application.

The analysis state is the state of the application when the user is analyzing the library.

Classes:
    AnalysisState: The analysis state of the application.
"""

from typing import Optional
from peplab.frontend.src.core.types.analysis_types import AnalysisType
from peplab.frontend.src.infrastructure.interfaces.state import State


class AnalysisState(State[AnalysisType]):
    """Manages analysis state and transitions."""

    def __init__(self) -> None:
        """Initialize the analysis state."""
        super().__init__()

    def is_valid_substate(self, value: AnalysisType) -> bool:
        """Validate analysis substate."""
        return isinstance(value, AnalysisType)

    def handle(self) -> None:
        """Handle analysis state behavior."""
        pass

    def get_route(self) -> str:
        """Get the current route based on state and substate."""
        base_route = "/analysis"
        if self.substate:
            return f"{base_route}/{self.substate.value}"
        return base_route

    def __str__(self) -> str:
        """Return string representation of the state."""
        return f"Analysis State: {self.substate.value}"
