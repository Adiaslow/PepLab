"""
This module defines the OptimizationState class, which is responsible for managing
the state of the optimization process.
"""

from typing import Optional

from peplab.frontend.src.infrastructure.interfaces.state import State
from peplab.frontend.src.core.types.optimization_types import OptimizationType


class OptimizationState(State):
    """
    State class for optimization process.
    Manages the current optimization type and associated data.
    """

    def __init__(self) -> None:
        """Initialize the optimization state."""
        super().__init__()
        self._substate: Optional[OptimizationType] = None

    @property
    def substate(self) -> Optional[OptimizationType]:
        """Get the current optimization substate.

        Returns:
            Optional[OptimizationType]: The current optimization type
        """
        return self._substate

    @substate.setter
    def substate(self, value: OptimizationType) -> None:
        """Set the current optimization substate.

        Args:
            value: The optimization type to set
        """
        if not isinstance(value, OptimizationType):
            raise ValueError("Substate must be an OptimizationType")
        self._substate = value

    def handle(self) -> None:
        """Handle the current optimization state."""
        if self._substate is None:
            return
        # Handle specific optimization type logic here

    def get_route(self) -> str:
        """Get the route for this state.

        Returns:
            str: The route for the optimization state
        """
        return "optimization.optimization"
