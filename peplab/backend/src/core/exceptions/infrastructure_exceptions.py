# peplab/backend/src/core/exceptions/infrustructure_exceptions.py

# External imports
from typing import List


class InfrastructureException(Exception):
    """The infrastructure exception."""

    def __init__(self, message: str) -> None:
        self.message: str = message
        super().__init__(self.message)


class EmptyBuildingBlockRepositoryError(InfrastructureException):
    """The empty building block repository error."""

    def __init__(self, message: str) -> None:
        super().__init__(message)
