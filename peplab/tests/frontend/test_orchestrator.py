# peplab/tests/frontend/test_orchestrator.py
"""
Tests for the ApplicationOrchestrator class.
"""

# Standard Library Imports
import pytest
from typing import List

# Third Party Imports
from pytest import MonkeyPatch

# Internal Imports
from peplab.frontend.src.infrastructure.orchestrator import (
    ApplicationOrchestrator,
    ApplicationContext,
)
from peplab.frontend.src.infrastructure.states.analysis_state import AnalysisState
from peplab.frontend.src.core.types.analysis_types import AnalysisType
from peplab.frontend.src.infrastructure.states.design_state import DesignState
from peplab.frontend.src.core.types.design_types import DesignType
from peplab.frontend.src.infrastructure.states.home_state import HomeState
from peplab.frontend.src.infrastructure.states.initialization_state import (
    InitializationState,
    InitializationPhase,
)
from peplab.frontend.src.infrastructure.states.dashboard_state import DashboardState
from peplab.frontend.src.infrastructure.services.backend_service import BackendService
from peplab.frontend.src.infrastructure.managers.state_manager import StateManager
from peplab.tests.frontend.mock_backend_service import MockBackendService


@pytest.fixture(autouse=True)
def reset_singletons() -> None:
    """Reset singleton instances before each test."""
    StateManager.reset()
    ApplicationOrchestrator.reset()


@pytest.fixture
def mock_backend_service(monkeypatch: MonkeyPatch) -> None:
    """Mock backend service to always return connected.

    Args:
        monkeypatch: The monkeypatch object to use.
    """

    # Create an instance of the mock service
    mock_service = MockBackendService()

    # Patch both the class and instance methods
    monkeypatch.setattr(BackendService, "__new__", lambda cls: mock_service)
    monkeypatch.setattr(
        BackendService, "verify_connection", mock_service.verify_connection
    )
    monkeypatch.setattr(BackendService, "is_connected", MockBackendService.is_connected)


@pytest.fixture
def orchestrator(
    mock_backend_service: MockBackendService, reset_singletons: None
) -> ApplicationOrchestrator:
    """Fixture to provide a fresh orchestrator instance for each test.

    Args:
        mock_backend_service: The mock backend service to use.
        reset_singletons: The reset singletons function to use.

    Returns:
        ApplicationOrchestrator: A fresh orchestrator instance.
    """
    orchestrator = ApplicationOrchestrator()
    orchestrator.initialize_application()
    return orchestrator


def test_singleton_pattern() -> None:
    """Test that orchestrator follows singleton pattern.

    Args:
        monkeypatch: The monkeypatch object to use.
    """
    orchestrator1 = ApplicationOrchestrator()
    orchestrator2 = ApplicationOrchestrator()
    assert orchestrator1 is orchestrator2


def test_initialization_flow(monkeypatch: MonkeyPatch) -> None:
    """Test the complete initialization flow."""
    orchestrator = ApplicationOrchestrator()

    # Initial state
    context: ApplicationContext = orchestrator.get_application_context()
    assert context.current_route == "/"
    assert context.current_state == "initialization"
    assert context.backend_status == "inactive"
    assert isinstance(orchestrator.state_manager.current_state, InitializationState)
    assert (
        orchestrator.state_manager.current_state.substate == InitializationPhase.START
    )
    assert orchestrator.state_manager.current_state.substate.value == "start"

    # Mock backend connection
    mock_service = MockBackendService()
    monkeypatch.setattr(BackendService, "__new__", lambda cls: mock_service)
    monkeypatch.setattr(
        BackendService, "verify_connection", mock_service.verify_connection
    )
    monkeypatch.setattr(BackendService, "is_connected", MockBackendService.is_connected)

    # After successful initialization
    orchestrator.initialize_application()
    context = orchestrator.get_application_context()
    assert context.current_route == "/"
    assert isinstance(orchestrator.state_manager.current_state, HomeState)
    assert context.current_state == "home"
    assert context.backend_status == "active"


def test_initialization_phases(monkeypatch: MonkeyPatch) -> None:
    """Test that initialization goes through all phases correctly."""
    orchestrator = ApplicationOrchestrator()

    # Should start in START phase
    assert isinstance(orchestrator.state_manager.current_state, InitializationState)
    assert (
        orchestrator.state_manager.current_state.substate == InitializationPhase.START
    )

    # Mock backend connection
    mock_service = MockBackendService()
    monkeypatch.setattr(BackendService, "__new__", lambda cls: mock_service)
    monkeypatch.setattr(
        BackendService, "verify_connection", mock_service.verify_connection
    )
    monkeypatch.setattr(BackendService, "is_connected", MockBackendService.is_connected)

    # Initialize application
    orchestrator.initialize_application()

    # Should end up in home state after successful initialization
    assert isinstance(orchestrator.state_manager.current_state, HomeState)


def test_home_to_dashboard_flow(orchestrator: ApplicationOrchestrator) -> None:
    """Test the flow from home to dashboard."""
    # Verify initial home state
    context: ApplicationContext = orchestrator.get_application_context()
    assert context.current_route == "/"
    assert context.current_state == "home"

    # Change to dashboard
    orchestrator.handle_state_change("DASHBOARD")
    context = orchestrator.get_application_context()
    assert context.current_route == "/dashboard"
    assert isinstance(orchestrator.state_manager.current_state, DashboardState)
    assert context.current_state == "dashboard"


def test_dashboard_to_design_flow(orchestrator: ApplicationOrchestrator) -> None:
    """Test the flow from dashboard to design state."""
    # Start at dashboard
    orchestrator.handle_state_change("DASHBOARD")
    context: ApplicationContext = orchestrator.get_application_context()
    assert context.current_route == "/dashboard"
    assert context.current_state == "dashboard"

    # Change to design hub
    orchestrator.handle_state_change("DESIGN")
    context = orchestrator.get_application_context()
    assert context.current_route == "/design"
    assert isinstance(orchestrator.state_manager.current_state, DesignState)
    assert context.current_state == "design"

    # Verify route history
    history = orchestrator.router.history
    assert "/dashboard" in history
    assert "/design" in history


def test_dashboard_to_analysis_flow(orchestrator: ApplicationOrchestrator) -> None:
    """Test the flow from dashboard to analysis state."""
    # Start at dashboard
    orchestrator.handle_state_change("DASHBOARD")
    context: ApplicationContext = orchestrator.get_application_context()
    assert context.current_route == "/dashboard"
    assert context.current_state == "dashboard"

    # Change to analysis hub
    orchestrator.handle_state_change("ANALYSIS")
    context = orchestrator.get_application_context()
    assert context.current_route == "/analysis"
    assert isinstance(orchestrator.state_manager.current_state, AnalysisState)
    assert context.current_state == "analysis"

    # Change to specific analysis type
    orchestrator.handle_state_change("CHEMINFORMATIC")
    context = orchestrator.get_application_context()
    assert context.current_route == "/analysis/cheminformatic"
    assert isinstance(orchestrator.state_manager.current_state, AnalysisState)
    assert context.current_state == "analysis"
    assert (
        orchestrator.state_manager.current_state.substate == AnalysisType.CHEMINFORMATIC
    )


def test_return_to_dashboard(orchestrator: ApplicationOrchestrator) -> None:
    """Test returning to dashboard from various states."""
    # Go to design state
    orchestrator.handle_state_change("DASHBOARD")
    orchestrator.handle_state_change("DESIGN")
    orchestrator.handle_state_change("COMBINATORIC")

    # Return to dashboard
    orchestrator.handle_state_change("DASHBOARD")
    context: ApplicationContext = orchestrator.get_application_context()
    assert context.current_route == "/dashboard"
    assert isinstance(orchestrator.state_manager.current_state, DashboardState)
    assert context.current_state == "dashboard"

    # Go to analysis state
    orchestrator.handle_state_change("ANALYSIS")
    orchestrator.handle_state_change("CHEMINFORMATIC")

    # Return to dashboard again
    orchestrator.handle_state_change("DASHBOARD")
    context = orchestrator.get_application_context()
    assert context.current_route == "/dashboard"
    assert isinstance(orchestrator.state_manager.current_state, DashboardState)
    assert context.current_state == "dashboard"


def test_direct_state_access_redirects_through_dashboard(
    orchestrator: ApplicationOrchestrator,
) -> None:
    """Test that accessing states directly properly routes through dashboard."""
    # Try to go directly to design type
    orchestrator.handle_state_change("COMBINATORIC")
    context: ApplicationContext = orchestrator.get_application_context()
    assert context.current_route == "/design/combinatoric"
    assert isinstance(orchestrator.state_manager.current_state, DesignState)
    assert context.current_state == "design"

    # Verify route history shows proper transition path
    history = orchestrator.router.history
    assert any(route == "/dashboard" for route in history)
    assert any(route == "/design" for route in history)
    assert history[-1] == "/design/combinatoric"


@pytest.mark.parametrize(
    "design_type",
    [design_type.value for design_type in DesignType],
)
def test_handle_design_state_changes(
    orchestrator: ApplicationOrchestrator, design_type: str
) -> None:
    """Test handling different design state changes.

    Args:
        orchestrator: The orchestrator to use.
        design_type: The design type to use.
    """
    orchestrator.handle_state_change(design_type)

    context: ApplicationContext = orchestrator.get_application_context()
    assert context.current_route == f"/design/{design_type.lower()}"
    assert isinstance(orchestrator.state_manager.current_state, DesignState)
    assert orchestrator.state_manager.current_state.substate == DesignType(
        design_type.lower()
    )


@pytest.mark.parametrize(
    "analysis_type",
    [analysis_type.value for analysis_type in AnalysisType],
)
def test_handle_analysis_state_changes(
    orchestrator: ApplicationOrchestrator, analysis_type: str
) -> None:
    """Test handling different analysis state changes.

    Args:
        orchestrator: The orchestrator to use.
        analysis_type: The analysis type to use.
    """
    orchestrator.handle_state_change(analysis_type)
    context: ApplicationContext = orchestrator.get_application_context()
    assert context.current_route == f"/analysis/{analysis_type.lower()}"
    assert isinstance(orchestrator.state_manager.current_state, AnalysisState)
    assert orchestrator.state_manager.current_state.substate == AnalysisType(
        analysis_type.lower()
    )


def test_application_context(orchestrator: ApplicationOrchestrator) -> None:
    """Test that application context is correctly updated.

    Args:
        orchestrator: The orchestrator to use.
    """
    # Test initial context (should be home state after initialization)
    context: ApplicationContext = orchestrator.get_application_context()
    assert context.current_route == "/"
    assert context.current_state == "home"
    assert context.backend_status == "active"

    # Test context after state change
    # First go to design hub
    orchestrator.handle_state_change("DESIGN")
    context = orchestrator.get_application_context()
    assert context.current_route == "/design"
    assert context.current_state == "design"
    assert context.backend_status == "active"

    # Then go to specific design type
    orchestrator.handle_state_change("COMBINATORIC")
    context = orchestrator.get_application_context()
    assert context.current_route == "/design/combinatoric"
    assert context.current_state == "design"
    assert context.backend_status == "active"


def test_invalid_state_change(orchestrator: ApplicationOrchestrator) -> None:
    """Test handling of invalid state changes.

    Args:
        orchestrator: The orchestrator to use.
    """
    initial_context: ApplicationContext = orchestrator.get_application_context()

    # Try to navigate to non-existent state
    orchestrator.handle_state_change("INVALID_STATE")

    # Context should remain unchanged
    current_context: ApplicationContext = orchestrator.get_application_context()
    assert current_context.current_route == initial_context.current_route
    assert current_context.current_state == initial_context.current_state


def test_state_manager_initialization(orchestrator: ApplicationOrchestrator) -> None:
    """Test that state manager is properly initialized.

    Args:
        orchestrator: The orchestrator to use.
    """
    assert orchestrator.state_manager is not None
    assert orchestrator.state_manager.current_state is not None
    assert isinstance(
        orchestrator.state_manager.current_state, HomeState
    )  # Now expects HomeState


def test_router_initialization(orchestrator: ApplicationOrchestrator) -> None:
    """Test that router is properly initialized.

    Args:
        orchestrator: The orchestrator to use.
    """
    assert orchestrator.router is not None
    assert orchestrator.router.current_route == "/"


def test_multiple_state_changes(orchestrator: ApplicationOrchestrator) -> None:
    """Test multiple sequential state changes.

    Args:
        orchestrator: The orchestrator to use.
    """
    # First go to design hub
    orchestrator.handle_state_change("DESIGN")
    context: ApplicationContext = orchestrator.get_application_context()
    assert context.current_route == "/design"
    assert context.current_state == "design"

    # Then test different design types
    states: List[str] = [design_type.value for design_type in DesignType]

    for state in states:
        orchestrator.handle_state_change(state)
        context: ApplicationContext = orchestrator.get_application_context()
        assert context.current_route == f"/design/{state.lower()}"
        assert isinstance(orchestrator.state_manager.current_state, DesignState)
        assert orchestrator.state_manager.current_state.substate == DesignType(
            state.lower()
        )


def test_home_to_design_flow(orchestrator: ApplicationOrchestrator) -> None:
    """Test the flow from home to design state.

    Args:
        orchestrator: The orchestrator to use.
    """
    # Verify initial home state
    context: ApplicationContext = orchestrator.get_application_context()
    assert context.current_route == "/"
    assert context.current_state == "home"

    # Change to design hub
    orchestrator.handle_state_change("DESIGN")
    context = orchestrator.get_application_context()
    assert context.current_route == "/design"
    assert isinstance(orchestrator.state_manager.current_state, DesignState)
    assert context.current_state == "design"

    # Change to specific design type
    orchestrator.handle_state_change("COMBINATORIC")
    context = orchestrator.get_application_context()
    assert context.current_route == "/design/combinatoric"
    assert isinstance(orchestrator.state_manager.current_state, DesignState)
    assert context.current_state == "design"


def test_home_to_analysis_flow(orchestrator: ApplicationOrchestrator) -> None:
    """Test the flow from home to analysis state.

    Args:
        orchestrator: The orchestrator to use.
    """
    # Verify initial home state
    context: ApplicationContext = orchestrator.get_application_context()
    assert context.current_route == "/"
    assert context.current_state == "home"

    # Change to analysis hub
    orchestrator.handle_state_change("ANALYSIS")
    context = orchestrator.get_application_context()
    assert context.current_route == "/analysis"
    assert isinstance(orchestrator.state_manager.current_state, AnalysisState)
    assert context.current_state == "analysis"

    # Change to specific analysis type
    orchestrator.handle_state_change("CHEMINFORMATIC")
    context = orchestrator.get_application_context()
    assert context.current_route == "/analysis/cheminformatic"
    assert isinstance(orchestrator.state_manager.current_state, AnalysisState)
    assert context.current_state == "analysis"


def test_multiple_analysis_state_changes(orchestrator: ApplicationOrchestrator) -> None:
    """Test multiple sequential analysis state changes."""
    # First go to analysis hub
    orchestrator.handle_state_change("ANALYSIS")
    context = orchestrator.get_application_context()
    assert context.current_route == "/analysis"
    assert context.current_state == "analysis"

    # Then test different analysis types
    analysis_states: List[str] = [
        analysis_type.value for analysis_type in AnalysisType
    ][
        :3
    ]  # Test first 3 analysis types

    for state in analysis_states:
        orchestrator.handle_state_change(state)
        context = orchestrator.get_application_context()
        assert context.current_route == f"/analysis/{state.lower()}"
        assert isinstance(orchestrator.state_manager.current_state, AnalysisState)
        assert orchestrator.state_manager.current_state.substate == AnalysisType(
            state.lower()
        )


def test_switching_between_design_and_analysis(
    orchestrator: ApplicationOrchestrator,
) -> None:
    """Test switching between design and analysis states."""
    # Start with design hub and type
    orchestrator.handle_state_change("DESIGN")
    orchestrator.handle_state_change("COMBINATORIC")
    context = orchestrator.get_application_context()
    assert context.current_route == "/design/combinatoric"
    assert isinstance(orchestrator.state_manager.current_state, DesignState)
    assert context.current_state == "design"

    # Switch to analysis hub and type
    orchestrator.handle_state_change("ANALYSIS")
    orchestrator.handle_state_change("MACHINE_LEARNING")
    context = orchestrator.get_application_context()
    assert context.current_route == "/analysis/machine_learning"
    assert isinstance(orchestrator.state_manager.current_state, AnalysisState)
    assert context.current_state == "analysis"

    # Back to design hub and type
    orchestrator.handle_state_change("DESIGN")
    orchestrator.handle_state_change("GENETIC")
    context = orchestrator.get_application_context()
    assert context.current_route == "/design/genetic"
    assert isinstance(orchestrator.state_manager.current_state, DesignState)
    assert context.current_state == "design"
