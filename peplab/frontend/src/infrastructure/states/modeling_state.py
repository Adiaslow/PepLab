# frontend/src/infrastructure/states/modeling_state.py
"""
This module defines the ModelingState class, which is responsible for managing
the state of the modeling process.
"""

from typing import Optional

from peplab.frontend.src.infrastructure.interfaces.state import State
from peplab.frontend.src.core.types.modeling_types import ModelingType


class ModelingState(State):
    """
    State class for modeling process.
    Manages the current modeling type and associated data.
    """

    def __init__(self) -> None:
        """Initialize the modeling state."""
        super().__init__()
        self._substate: Optional[ModelingType] = None

    @property
    def substate(self) -> Optional[ModelingType]:
        """Get the current modeling substate.

        Returns:
            Optional[ModelingType]: The current modeling type
        """
        return self._substate

    @substate.setter
    def substate(self, value: ModelingType) -> None:
        """Set the current modeling substate.

        Args:
            value: The modeling type to set
        """
        if not isinstance(value, ModelingType):
            raise ValueError("Substate must be a ModelingType")
        self._substate = value

    def handle(self) -> None:
        """Handle the current modeling state."""
        if self._substate is None:
            return
        # Handle specific modeling type logic here

    def get_route(self) -> str:
        """Get the route for this state.

        Returns:
            str: The route for the modeling state
        """
        return "modeling.modeling"
