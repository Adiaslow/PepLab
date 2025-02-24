# peplab/frontend/src/core/exceptions/backend_communication_exceptions.py
"""
This module contains the exceptions for the backend communication.
"""


class BackendConnectionError(Exception):
    """The backend connection error."""

    def __init__(self, message: str) -> None:
        self.message: str = message
        super().__init__(self.message)
