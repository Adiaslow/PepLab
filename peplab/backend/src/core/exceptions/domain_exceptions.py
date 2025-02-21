# peplab/backend/src/core/exceptions/domain_exceptions.py
"""
This module contains the domain exceptions.

Classes:
    DomainException: The base domain exception.
    InvalidDirectionError: The invalid direction error.
    InvalidBuildingBlockError: The invalid building block error.
    InvalidPeptideError: The invalid peptide error.
"""

# External imports
from abc import ABC
from typing import Any


# Base domain exception
class DomainException(ABC, Exception):
    """The base domain exception."""

    def __init__(self, message: str) -> None:
        self.message: str = message
        super().__init__(self.message)


# Building block exceptions
class InvalidBuildingBlockError(DomainException):
    """The invalid building block error."""

    def __init__(self, message: str) -> None:
        super().__init__(message)


class InvalidBuildingBlockIndexError(DomainException):
    """The invalid building block index error."""

    def __init__(self, message: str) -> None:
        super().__init__(message)


# Peptide exceptions
class InvalidPeptideError(DomainException):
    """The invalid peptide error."""

    def __init__(self, message: str) -> None:
        super().__init__(message)


class InvalidDirectionError(DomainException):
    """The invalid direction error."""

    def __init__(self, message: str) -> None:
        super().__init__(message)


class NoBuildingBlocksError(DomainException):
    """The no building blocks error."""

    def __init__(self, message: str) -> None:
        super().__init__(message)


class NoPropertyOrMetadataError(DomainException):
    """The no property or metadata error."""

    def __init__(self, message: str) -> None:
        super().__init__(message)


# Library exceptions
class InvalidLibraryError(DomainException):
    """The invalid library error."""

    def __init__(self, message: str) -> None:
        super().__init__(message)


class InvalidLibraryIndexError(DomainException):
    """The invalid library index error."""

    def __init__(self, message: str) -> None:
        super().__init__(message)


class InvalidPeptideNameOrIDError(DomainException):
    """The invalid peptide name or ID error."""

    def __init__(self, message: str) -> None:
        super().__init__(message)


class InvalidBuildingBlocksError(DomainException):
    """The invalid building blocks error."""

    def __init__(self, message: str) -> None:
        super().__init__(message)


class InvalidPeptideUUIDError(DomainException):
    """The invalid peptide UUID error."""

    def __init__(self, message: str) -> None:
        super().__init__(message)
