# peplab/frontend/src/application/states/home_state.py
"""
This module is responsible for handling the home state or splash screen of the application.

Classes:
    HomeState: The home state of the application.
"""

# Internal imports
from ..interfaces.state import State
from typing import Optional, Literal


HomeSubstate = Literal["main", "about", "contact"]


class HomeState(State[HomeSubstate]):
    """
    Represents the home/splash screen state of the application.
    """

    def __init__(self) -> None:
        """Initialize the home state."""
        super().__init__()
        self.substate: HomeSubstate = "main"

    def is_valid_substate(self, value: HomeSubstate) -> bool:
        """Validate home substate."""
        return value in ["main", "about", "contact"]

    def handle(self) -> None:
        """Handle any home state specific logic."""
        pass

    def get_route(self) -> str:
        """Get the current route based on state and substate."""
        if self.substate and self.substate != "main":
            return f"/{self.substate}"
        return "/"

    def __str__(self) -> str:
        """Return string representation of the state."""
        return f"Home State: {self.substate}"
