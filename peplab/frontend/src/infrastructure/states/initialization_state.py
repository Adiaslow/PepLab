# peplab/frontend/src/application/states/initialization_state.py
"""
This module is responsible for handling the initialization state of the application.

The initialization state is the first state of the application.
It is responsible for initializing the infrastructure of the application.

Classes:
    InitializationState: The initialization state of the application.
"""

# Internal imports
from peplab.frontend.src.infrastructure.interfaces.state import State
from peplab.frontend.src.core.types.initialization_types import InitializationPhase
from peplab.frontend.src.infrastructure.services.backend_service import BackendService


class InitializationState(State[InitializationPhase]):
    """Manages initialization state and transitions.

    Attributes:
        substate: The current phase of the initialization state.
        _backend_service: The backend service for the application.

    Methods:
        is_valid_substate: Validate initialization phase.
        handle: Handle initialization state behavior.
        get_route: Get the current route based on state and substate.
        __str__: Return string representation of the state.
    """

    def __init__(self) -> None:
        """Initialize the initialization state."""
        super().__init__()
        self.substate: InitializationPhase = InitializationPhase.START
        self._backend_service = BackendService()

    def is_valid_substate(self, value: InitializationPhase) -> bool:
        """Validate initialization phase."""
        return isinstance(value, InitializationPhase)

    def handle(self) -> None:
        """Handle initialization state behavior."""
        self._backend_service.verify_connection()

    def get_route(self) -> str:
        """Get the current route based on state and substate."""
        return "/"

    def __str__(self) -> str:
        """Return string representation of the state."""
        return (
            f"Initialization State: {self.substate.value if self.substate else 'none'}"
        )
