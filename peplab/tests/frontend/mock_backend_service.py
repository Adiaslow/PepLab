# peplab/tests/frontend/mock_backend_service.py
"""
This module contains the MockBackendService class.

Classes:
    MockBackendService: A mock backend service that always returns connected.
"""


class MockBackendService:
    """Mock backend service to always return connected.

    Attributes:
        _connection_status: The connection status of the backend service.

    Methods:
        __init__: Initialize the mock backend service.
        verify_connection: Verify the connection to the backend service.
        is_connected: Check if the backend service is connected.
    """

    def __init__(self) -> None:
        """Initialize the mock backend service."""
        self._connection_status = True  # Initialize as connected

    def verify_connection(self) -> bool:
        """Verify the connection to the backend service.

        Returns:
            bool: Always returns True for testing.
        """
        return True

    @property
    def is_connected(self) -> bool:
        """Check if the backend service is connected.

        Returns:
            bool: Always returns True for testing.
        """
        return True
