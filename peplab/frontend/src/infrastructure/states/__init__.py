# peplab/frontend/src/infrastructure/states/__init__.py
"""
This module contains the states of the frontend application.

States act as the interface between the user-side interactions with the frontend
and the backend. Each state corresponds to a discrete page or set of pages relating
to a set of backend processes.
"""
# External import
from typing import Any, List

# Internal imports
from peplab.frontend.src.infrastructure.states.analysis_state import AnalysisState
from peplab.frontend.src.infrastructure.states.dashboard_state import DashboardState
from peplab.frontend.src.infrastructure.states.design_state import DesignState
from peplab.frontend.src.infrastructure.states.home_state import HomeState
from peplab.frontend.src.infrastructure.states.initialization_state import (
    InitializationState,
)

__all__: List[Any] = [
    "AnalysisState",
    "DashboardState",
    "DesignState",
    "HomeState",
    "InitializationState",
]
