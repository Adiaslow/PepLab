# peplab/frontend/src/core/exceptions/state_exceptions.py
"""
This module contains the state exceptions.
"""

# Standard Library Imports
from abc import ABC


# State Exceptions
class StateException(ABC, Exception):
    """The state exception."""

    def __init__(self, message: str) -> None:
        self.message: str = message
        super().__init__(self.message)


class InvalidStateError(StateException):
    """The invalid state error."""

    def __init__(self, message: str) -> None:
        super().__init__(message)


class InvalidSubstateError(StateException):
    """The invalid substate error."""

    def __init__(self, message: str) -> None:
        super().__init__(message)


class NoSubstateError(StateException):
    """The no substate error."""

    def __init__(self, message: str) -> None:
        super().__init__(message)
