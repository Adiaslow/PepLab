# peplab/frontend/src/application/interfaces/state.py
"""
This module is responsible for defining the State interface.

Classes:
    State: The State interface.
"""

# External imports
from abc import ABC, abstractmethod


class State(ABC):
    """The State interface.

    Methods:
        handle: Handle the state.
    """

    @abstractmethod
    def handle(self) -> None:
        """Handle the state.

        This method is responsible for handling the state.
        """
        ...
