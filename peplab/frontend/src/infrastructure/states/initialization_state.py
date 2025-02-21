# peplab/frontend/src/application/states/initialization_state.py
"""
This module is responsible for handling the initialization state of the application.

The initialization state is the first state of the application.
It is responsible for initializing the infrastructure of the application.

Classes:
    InitializationState: The initialization state of the application.
"""

# Internal imports
from ..interfaces.state import State


class InitializationState(State):
    """The initialization state of the application.

    Methods:
        handle: Handle the initialization state.
    """

    def __init__(self) -> None:
        """Initialize the initialization state.

        This method is responsible for initializing the initialization state.
        """
        super().__init__()

    def handle(self) -> None:
        """Handle the initialization state.

        This method is responsible for handling the initialization state.
        """
        ...
