# peplab/frontend/src/infrastructure/managers/__init__.py
"""
This module contains the managers of the infrastructure of the frontend application.
"""

# External imports
from typing import Any, List

# Internal imports
from .state_manager import StateManager

__all__: List[Any] = [
    "StateManager",
]
