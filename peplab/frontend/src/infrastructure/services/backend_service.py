"""
Service for handling backend communication.
"""

from typing import Optional
import requests
from peplab.frontend.src.core.exceptions.backend_communication_exceptions import (
    BackendConnectionError,
)


class BackendService:
    """Service for managing backend communication and status."""

    def __init__(self):
        """Initialize backend service with default configuration."""
        self._backend_url = "http://localhost:5000"  # Configure as needed
        self._connection_status: Optional[bool] = None

    def verify_connection(self) -> bool:
        """
        Verify connection to the backend service.

        Returns:
            bool: True if connection is successful, False otherwise

        Raises:
            BackendConnectionError: If connection attempt fails
        """
        try:
            response = requests.get(f"{self._backend_url}/health")
            self._connection_status = response.status_code == 200
            return self._connection_status
        except requests.RequestException as e:
            self._connection_status = False
            raise BackendConnectionError(f"Failed to connect to backend: {str(e)}")

    @property
    def is_connected(self) -> bool:
        """Return current backend connection status."""
        if self._connection_status is None:
            try:
                return self.verify_connection()
            except BackendConnectionError:
                return False
        return self._connection_status
