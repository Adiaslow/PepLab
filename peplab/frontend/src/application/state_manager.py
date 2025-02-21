# peplab/frontend/src/application/state_manager.py
"""
This module is responsible for managing the state of the application.
It is responsible for setting the state of the application and handling the state.

Classes:
    StateManager: Manages the state of the application.
"""

# Internal imports
from .interfaces.state import State


class StateManager:
    """This class is responsible for managing the state of the application.
    It is responsible for setting the state of the application and handling the state.

    Attributes:
        _state: The current state of the application.
        _state_history: The history of the states.

    Methods:
        __init__: Initialize the StateManager.
        set_state: Set the state of the application.
        get_state: Get the current state of the application.
        state: Get the current state of the application.
        state_history: Get the history of the states.
        handle: Handle the current state.
        _validate_state: Validate the state.
    """

    def __init__(self, state: State | None = None) -> None:
        """Initialize the StateManager.

        Args:
            state: The initial state of the application.

        Raises:
            ValueError: If the state is not a valid state.
        """
        if state is not None and not isinstance(state, State):
            raise ValueError("State must be a valid state.")
        self.__state: State | None = state
        self.__state_history: list[State] = []

    def set_state(self, state: State) -> None:
        """Set the state of the application.

        This method is responsible for setting the state of the application.
        It will validate the state and set the state of the application.

        Args:
            state: The state to set the application to.

        Raises:
            ValueError: If the state is not a valid state.
        """
        try:
            self.__validate_state(state)
            self.__state = state
        except:
            ...

    def get_state(self) -> State | None:
        """Get the current state of the application.

        Returns:
            The current state of the application.

        Raises:
            ValueError: If the state is not set.
        """
        if self.__state is None:
            raise ValueError("State is not set.")
        return self.__state

    @property
    def state(self) -> State | None:
        """Get the current state of the application.

        Returns:
            The current state of the application.
        """
        return self.__state

    @property
    def state_history(self) -> list[State]:
        """Get the history of the states.

        Returns:
            The history of the states.
        """
        return self.__state_history

    def handle(self) -> None:
        """Handle the current state.

        This method is responsible for handling the current state.

        Raises:
            ValueError: If the state is not set.
        """
        if self.__state is not None:
            self.__state.handle()
        else:
            raise ValueError("State is not set.")

    def __validate_state(self, state: State) -> bool:
        """Validate the state.

        Args:
            state: The state to validate.

        Returns:
            True if the state is valid, False otherwise.

        Raises:
            ValueError: If the state is not a valid state.
        """
        if not isinstance(state, State):
            raise ValueError("State must be a valid state.")
        return True
