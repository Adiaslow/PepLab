# peplab/frontend/src/application/interfaces/state.py
"""
This module is responsible for defining the State interface.

Classes:
    State: The State interface.
"""

# Standard Library Imports
from abc import ABC, abstractmethod
from typing import Any, List, Union

# Internal Imports
from core.exceptions.state_exceptions import InvalidSubstateError, NoSubstateError


class State(ABC):
    """The State interface.

    Methods:
        handle: Handle the state.
    """

    def __init__(self, substate: Any | None = None) -> None:
        """Initialize the state.

        This method is responsible for initializing the state.
        """
        self.substates: List[Any] = []
        self.substate: Union[Any, None] = substate
        self.substate_history: List[Any] = []

    @abstractmethod
    def handle(self) -> None:
        """Handle the state.

        This method is responsible for handling the state.
        """
        ...

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
        self.substate_history.append(self.substate)
        self.substate = substate

    def __validate_substate(self, substate: Any) -> bool:
        """Validate the substate.

        This method is responsible for validating the substate.

        Args:
            substate: The substate to validate.

        Returns:
            bool: True if the substate is valid, False otherwise.
        """
        if substate not in self.substates:
            raise InvalidSubstateError(f"Invalid substate: {substate}")
        return True
