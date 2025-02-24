# peplab/frontend/src/infrastructure/managers/state_manager.py
"""
This module contains the StateManager class which manages application state transitions.
"""

from typing import Optional
from peplab.frontend.src.infrastructure.interfaces.state import State
from peplab.frontend.src.infrastructure.states.initialization_state import (
    InitializationState,
)


class StateManager:
    """
    Manages application state transitions using the State pattern.
    Implements Singleton pattern to ensure only one state manager exists.
    """

    _instance: Optional["StateManager"] = None
    _current_state: Optional[State] = None

    def __new__(cls) -> "StateManager":
        """Create or return the singleton instance."""
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._instance._current_state = InitializationState()
        return cls._instance

    def __init__(self) -> None:
        """Initialize with default state."""
        if not getattr(self, "initialized", False):
            self._current_state: Optional[State] = None
            # Set initial state to InitializationState
            self.set_state(InitializationState())
            self.initialized = True

    @property
    def current_state(self) -> Optional[State]:
        """Get the current state.

        Returns:
            Optional[State]: The current state
        """
        return self._current_state

    def set_state(self, state: State) -> None:
        """Set the current state.

        Args:
            state: The state to set
        """
        self._current_state = state

    def handle(self) -> None:
        """Handle the current state."""
        if self._current_state:
            self._current_state.handle()

    @classmethod
    def reset(cls) -> None:
        """Reset the singleton instance."""
        cls._instance = None
