# peplab/frontend/src/application/states/design_state.py
"""
This module is responsible for handling the design state of the application.

The design state is the state of the application when the user is designing the library.
Classes:
    DesignState: The design state of the application.
"""

# External imports
from enum import Enum

# Internal imports
from ..interfaces.state import State
from core.types.design_types import DesignType


class DesignState(State):
    """The design state of the application.

    Methods:
        handle: Handle the design state.
    """

    def __init__(self) -> None:
        """Initialize the design state.

        This method is responsible for initializing the design state.
        """
        super().__init__()

    def handle(self) -> None:
        """Handle the design state.

        This method is responsible for handling the design state.
        """
        ...

    def __validate_design_type(self, design_type: DesignType) -> bool:
        """Validate the design type.

        This method is responsible for validating the design type.
        """
        return design_type in DesignType
