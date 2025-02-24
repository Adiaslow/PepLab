# peplab/frontend/src/infrastructure/orchestrator.py
"""
This module serves as the main orchestrator for the application, coordinating
between the state manager, router, and backend services.

Classes:
    ApplicationContext: Represents the current state of the application.
    ApplicationOrchestrator: Orchestrates the application state and routing.
"""

from typing import Optional
from dataclasses import dataclass
from peplab.frontend.src.infrastructure.managers.state_manager import StateManager
from peplab.frontend.src.infrastructure.states.initialization_state import (
    InitializationState,
    InitializationPhase,
)
from peplab.frontend.src.infrastructure.states.home_state import HomeState
from peplab.frontend.src.infrastructure.states.design_state import DesignState
from peplab.frontend.src.infrastructure.router import Router
from peplab.frontend.src.infrastructure.services.backend_service import BackendService
from peplab.frontend.src.core.types.design_types import DesignType
from peplab.frontend.src.infrastructure.states.analysis_state import AnalysisState
from peplab.frontend.src.core.types.analysis_types import AnalysisType
from peplab.frontend.src.infrastructure.interfaces.state import State
from peplab.frontend.src.infrastructure.states.dashboard_state import DashboardState


@dataclass
class ApplicationContext:
    """Represents the current state of the application.

    Attributes:
        current_route: The current route of the application.
        current_state: The current state of the application.
        backend_status: The status of the backend connection.
    """

    current_route: str = "/"
    current_state: str = "initialization"  # Start at initialization state
    backend_status: str = "inactive"  # Initial backend status


class ApplicationOrchestrator:
    """Orchestrates the application state and routing.

    This class is responsible for coordinating the application state and routing.
    It is responsible for initializing the application, handling state changes,
    and getting the current application context.

    Attributes:
        _instance: The singleton instance of the ApplicationOrchestrator.
        _initialized: Whether the application has been initialized.
        _context: The current application context.
        router: The router for the application.
        state_manager: The state manager for the application.
        backend_service: The backend service for the application.

    Methods:
        __new__: The singleton pattern implementation.
        __init__: The constructor for the ApplicationOrchestrator.
        initialize_application: Initializes the application.
        handle_state_change: Handles state changes in the application.
        get_application_context: Gets the current application context.
        reset: Resets the singleton instance.
    """

    _instance: Optional["ApplicationOrchestrator"] = None

    def __new__(cls) -> "ApplicationOrchestrator":
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._instance._initialized = False
        return cls._instance

    def __init__(self) -> None:
        """Initialize orchestrator with default configuration."""
        if not self._initialized:
            self._context = ApplicationContext()
            self.router = Router()
            self.state_manager = (
                StateManager()
            )  # This now sets InitializationState by default
            self.backend_service = BackendService()
            self._initialized = False

    def initialize_application(self) -> None:
        """Initialize the application state with initialization checks."""
        if not self._initialized:
            try:
                # Start with initialization state
                if not self.state_manager.current_state:
                    init_state = InitializationState()
                    init_state.substate = InitializationPhase.START
                    self.state_manager.set_state(init_state)

                # Perform initialization checks
                init_state = self.state_manager.current_state
                if isinstance(init_state, InitializationState):
                    init_state.substate = InitializationPhase.CHECKING
                    init_state.handle()

                    if self.backend_service.is_connected:
                        # Mark initialization as complete
                        init_state.substate = InitializationPhase.COMPLETE

                        # Transition to home/splash state
                        self.state_manager.set_state(HomeState())
                        self.router.navigate_to("/")
                        self._initialized = True
                    else:
                        init_state.substate = InitializationPhase.FAILED
                        raise RuntimeError(
                            "Backend connection failed during initialization"
                        )
            except Exception as e:
                print(f"Initialization failed: {str(e)}")
                self._initialized = False

    def handle_state_change(self, new_state: str) -> None:
        """Handle state changes in the application."""
        if not self.backend_service.is_connected:
            return

        try:
            new_state_lower = new_state.lower()

            # Handle main navigation flows
            if new_state_lower == "home":
                home_state = HomeState()
                self._transition_to_hub_state(home_state)
                return

            if new_state_lower == "dashboard":
                dashboard_state = DashboardState()
                self._transition_to_hub_state(dashboard_state)
                return

            # Handle workflow states through dashboard
            if new_state_lower == "design":
                # Ensure we're in dashboard first
                if not isinstance(self.state_manager.current_state, DashboardState):
                    self._transition_to_hub_state(DashboardState())
                design_state = DesignState()
                self._transition_to_hub_state(design_state)
                return

            if new_state_lower == "analysis":
                # Ensure we're in dashboard first
                if not isinstance(self.state_manager.current_state, DashboardState):
                    self._transition_to_hub_state(DashboardState())
                analysis_state = AnalysisState()
                self._transition_to_hub_state(analysis_state)
                return

            # Handle specific design types
            if new_state_lower in [state.value for state in DesignType]:
                self._ensure_hub_state(DesignState)
                if isinstance(self.state_manager.current_state, DesignState):
                    self.state_manager.current_state.substate = DesignType(
                        new_state_lower
                    )
                    self.router.navigate_to(
                        self.state_manager.current_state.get_route()
                    )

            # Handle specific analysis types
            elif new_state_lower in [state.value for state in AnalysisType]:
                self._ensure_hub_state(AnalysisState)
                if isinstance(self.state_manager.current_state, AnalysisState):
                    self.state_manager.current_state.substate = AnalysisType(
                        new_state_lower
                    )
                    self.router.navigate_to(
                        self.state_manager.current_state.get_route()
                    )

        except ValueError:
            # Invalid state - no change
            pass

    def _ensure_hub_state(self, state_class: type) -> None:
        """Ensure we're in the correct hub state."""
        if not isinstance(self.state_manager.current_state, state_class):
            # First ensure we're in dashboard
            if not isinstance(self.state_manager.current_state, DashboardState):
                dashboard_state = DashboardState()
                self._transition_to_hub_state(dashboard_state)
            # Then transition to requested hub
            new_state = state_class()
            self._transition_to_hub_state(new_state)

    def _transition_to_hub_state(self, state: State) -> None:
        """Transition to a hub state.

        Args:
            state: The state to transition to
        """
        self.state_manager.set_state(state)
        self.router.navigate_to(state.get_route())

    def get_application_context(self) -> ApplicationContext:
        """Get the current application context.

        Returns:
            ApplicationContext: The current application context.
        """
        # Ensure we have a state
        if not self.state_manager.current_state:
            self.state_manager.set_state(InitializationState())

        # Get the class name, remove 'State' suffix and convert to lowercase
        current_state: str = self.state_manager.current_state.__class__.__name__.lower()
        current_state: str = current_state.replace("state", "")

        # Determine backend status
        backend_status: str = (
            "active" if self.backend_service.is_connected else "inactive"
        )

        context: ApplicationContext = ApplicationContext(
            current_route=self.router.current_route,
            current_state=current_state,
            backend_status=backend_status,
        )
        return context

    @classmethod
    def reset(cls) -> None:
        """Reset the singleton instance."""
        cls._instance = None
