# infrastructure/states/dashboard_state.py
"""
This module handles the dashboard state of the application.

Classes:
    DashboardState: The dashboard state of the application.
"""

from typing import Optional, Literal
from peplab.frontend.src.infrastructure.interfaces.state import State


DashboardSubstate = Literal["overview", "settings", "history"]


class DashboardState(State[DashboardSubstate]):
    """Manages dashboard state and transitions."""

    def __init__(self) -> None:
        """Initialize the dashboard state."""
        super().__init__()
        self.substate: DashboardSubstate = "overview"

    def is_valid_substate(self, value: DashboardSubstate) -> bool:
        """Validate dashboard substate."""
        return value in ["overview", "settings", "history"]

    def handle(self) -> None:
        """Handle dashboard state behavior."""
        pass

    def get_route(self) -> str:
        """Get the current route based on state and substate."""
        base_route = "/dashboard"
        if self.substate and self.substate != "overview":
            return f"{base_route}/{self.substate}"
        return base_route

    def __str__(self) -> str:
        """Return string representation of the state."""
        return f"Dashboard State: {self.substate}"
