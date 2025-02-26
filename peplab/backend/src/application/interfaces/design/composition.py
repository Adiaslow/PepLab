# peplab/backend/src/application/interfaces/design/composition.py
"""
This module defines the Composition interface for combinatorial and group-theoretic compositions.

Classes:
    Composition: Abstract base class for combinatorial and group-theoretic compositions.
"""

# External imports
from abc import ABC, abstractmethod
from typing import List, Any


class Composition(ABC):
    """Abstract base class for combinatorial and group-theoretic compositions.

    Attributes:
        elements: List[Any] - The elements of the composition.
        filtered_elements: List[Any] - The filtered elements of the composition.
        unique_elements: List[Any] - The unique elements of the composition.
        ordered_elements: List[Any] - The ordered elements of the composition.
    """

    def generate_composition(self, *args, **kwargs) -> List[Any]:
        """Template method to generate compositions.

        Calls abstract `_generate_elements`, and optionally filters, deduplicates, and orders results.

        Returns:
            List[Any] - The generated composition.
        """
        elements: List[Any] = self._generate_elements(*args, **kwargs)
        filtered_elements: List[Any] = self._apply_filter(elements, *args, **kwargs)
        unique_elements: List[Any] = self._deduplicate(
            filtered_elements, *args, **kwargs
        )
        return self._order_elements(unique_elements, *args, **kwargs)

    @abstractmethod
    def _generate_elements(self, *args, **kwargs) -> List[Any]:
        """Abstract method to be implemented by subclasses to generate elements.

        Returns:
            List[Any] - The generated elements.
        """
        pass

    def _apply_filter(self, elements: List[Any], *args, **kwargs) -> List[Any]:
        """Optional: Filter elements.

        Returns:
            List[Any] - The filtered elements.
        """
        return elements  # Default: no filtering

    def _deduplicate(self, elements: List[Any], *args, **kwargs) -> List[Any]:
        """Optional: Deduplicate elements.

        Returns:
            List[Any] - The deduplicated elements.
        """
        return elements  # Default: no deduplication

    def _order_elements(self, elements: List[Any], *args, **kwargs) -> List[Any]:
        """Optional: Order elements.

        Returns:
            List[Any] - The ordered elements.
        """
        return elements  # Default: no ordering
