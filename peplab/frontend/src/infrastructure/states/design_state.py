# peplab/frontend/src/application/states/design_state.py
"""
This module contains the DesignState class which manages design-related state.
"""

from typing import Optional
from peplab.frontend.src.core.types.design_types import DesignType
from peplab.frontend.src.infrastructure.interfaces.state import State


class DesignState(State[DesignType]):
    """
    Represents the design state of the application.
    Handles design-specific behavior and transitions.
    """

    def __init__(self) -> None:
        """Initialize the design state."""
        super().__init__()

    def is_valid_substate(self, value: DesignType) -> bool:
        """Validate design substate."""
        return isinstance(value, DesignType)

    def handle(self) -> None:
        """Handle design state behavior."""
        # Implementation of design state handling
        pass

    def get_route(self) -> str:
        """Get the current route based on state and substate."""
        base_route = "/design"
        if self.substate:
            return f"{base_route}/{self.substate.value}"
        return base_route

    def __str__(self) -> str:
        """Return string representation of the state."""
        return f"Design State: {self.substate.value}"
