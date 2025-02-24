# peplab/frontend/src/infrastructure/router.py
"""
This module maintains route state information for the frontend application.
Actual routing is handled by the web framework in presentation/app/routes/.
"""

# Internal imports
from peplab.frontend.src.core.types.design_types import DesignType
from typing import List


class Router:
    """Router for handling navigation in the application.

    Attributes:
        _current_route: The current route of the application.
        _route_history: The history of routes visited.

    Methods:
        __init__: Initialize router with default route.
        current_route: Get the current route.
        history: Get the history of routes visited.
        navigate_to: Navigate to a new route.
        handle_home: Handle home route.
        handle_design: Handle design route.
    """

    def __init__(self) -> None:
        """Initialize router with default route."""
        self._current_route: str = "/"
        self._route_history: List[str] = ["/"]

    @property
    def current_route(self) -> str:
        """Get current route.

        Returns:
            str: The current route.
        """
        return self._current_route

    @property
    def history(self) -> List[str]:
        """Get route history.

        Returns:
            List[str]: The history of routes visited.
        """
        return self._route_history.copy()

    def navigate_to(self, route: str) -> None:
        """Navigate to a new route.

        Args:
            route: The route to navigate to.
        """
        self._current_route = route
        self._route_history.append(route)

    def handle_home(self) -> None:
        """Update route state to home.

        Returns:
            None
        """
        self._current_route = "/"

    def handle_design(self, design_type: DesignType) -> None:
        """Update route state for design pages

        Args:
            design_type (DesignType): The design type to navigate to
        """
        self._current_route = f"/{design_type.value.lower()}"
