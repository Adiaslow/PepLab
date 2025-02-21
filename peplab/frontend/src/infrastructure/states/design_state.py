# peplab/frontend/src/application/states/design_state.py
"""
This module is responsible for handling the design state of the application.

The design state is the state of the application when the user is designing the library.
Classes:
    DesignState: The design state of the application.
"""

# Standard Library Imports
from enum import Enum
from typing import List

# Internal Imports
from ..interfaces.state import State
from core.types.design_types import DesignType
from infrastructure.managers.state_manager import StateManager


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
        self.substates: List[DesignType] = [substate for substate in DesignType]

    def handle(self) -> None:
        """Handle the design state.

        This method is responsible for handling the design state.
        """
        match self.substate:
            case DesignType.COMBINATORIC:
                # Initialize combinatoric backend scripts
                ...
            case DesignType.GENERATIVE:
                # Initialize generative backend scripts
                ...
            case DesignType.GENETIC:
                # Initialize genetic backend scripts
                ...
            case DesignType.MCMC:
                # Initialize MCMC backend scripts
                ...
            case DesignType.RANDOM:
                # Initialize random backend scripts
                ...
            case _:
                raise ValueError(f"Invalid design type: {self.substate}")


state_manager = StateManager(DesignState())
