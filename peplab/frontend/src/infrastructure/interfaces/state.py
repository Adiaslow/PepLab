# peplab/frontend/src/application/interfaces/state.py
"""
This module defines the base State interface with hierarchical state support.

Classes:
    State: The State interface.
"""

# Standard Library Imports
from abc import ABC, abstractmethod
from typing import Any, List, Optional, TypeVar, Generic

# Internal Imports
from peplab.frontend.src.core.exceptions.state_exceptions import (
    InvalidSubstateError,
    NoSubstateError,
)

T = TypeVar("T")  # For type of substate (e.g., DesignType, AnalysisType)


class State(ABC, Generic[T]):
    """Base State interface with hierarchical state support."""

    def __init__(self) -> None:
        """Initialize state with substate support."""
        self._substate: Optional[T] = None
        self._parent_state: Optional["State"] = None
        self._child_states: List["State"] = []

    @property
    def substate(self) -> Optional[T]:
        """Get current substate safely."""
        return self._substate

    @substate.setter
    def substate(self, value: T) -> None:
        """Set substate with validation."""
        if not self.is_valid_substate(value):
            raise ValueError(f"Invalid substate {value} for {self.__class__.__name__}")
        self._substate = value

    def is_valid_substate(self, value: T) -> bool:
        """Validate if a substate is valid for this state."""
        return True  # Override in concrete classes

    @property
    def parent_state(self) -> Optional["State"]:
        """Get parent state."""
        return self._parent_state

    @parent_state.setter
    def parent_state(self, state: "State") -> None:
        """Set parent state."""
        self._parent_state = state

    def add_child_state(self, state: "State") -> None:
        """Add a child state."""
        state.parent_state = self
        self._child_states.append(state)

    def get_state_path(self) -> List[str]:
        """Get the full path of states from root to this state."""
        path = []
        current = self
        while current:
            state_name = current.__class__.__name__.lower().replace("state", "")
            path.insert(0, state_name)
            current = current.parent_state
        return path

    @abstractmethod
    def handle(self) -> None:
        """Handle the state's behavior."""
        pass

    def get_substate(self) -> Any:
        """Get the substate.

        This method is responsible for getting the substate.
        """
        if self.substate is None:
            raise NoSubstateError("No substate found")
        return self.substate

    def set_substate(self, substate: Any) -> None:
        """Set the substate.

        This method is responsible for setting the substate.

        Args:
            substate: The substate to set.

        Raises:
            InvalidSubstateError: If the substate is invalid.
        """
        self.__validate_substate(substate)
        self._substate = substate

    def __validate_substate(self, substate: Any) -> bool:
        """Validate the substate.

        This method is responsible for validating the substate.

        Args:
            substate: The substate to validate.

        Returns:
            bool: True if the substate is valid, False otherwise.
        """
        if substate not in self._child_states:
            raise InvalidSubstateError(f"Invalid substate: {substate}")
        return True

    @abstractmethod
    def get_route(self) -> str:
        """Get the current route for this state.

        Returns:
            str: The route for this state
        """
        pass
