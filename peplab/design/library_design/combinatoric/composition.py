from abc import ABC, abstractmethod

class Composition(ABC):
    """
    Abstract base class for combinatorial and group-theoretic compositions.
    Uses the Template Pattern to define a structure for generating, filtering, 
    and deduplicating compositions.
    """

    def generate_composition(self):
        """Template method to generate compositions."""
        elements = self._generate_elements()
        filtered_elements = self._apply_filter(elements)
        unique_elements = self._deduplicate(filtered_elements)
        return self._order_elements(unique_elements)

    @abstractmethod
    def _generate_elements(self):
        """Abstract method to be implemented by subclasses to generate elements."""
        pass

    def _apply_filter(self, elements):
        """Hook method for filtering elements (optional, can be overridden)."""
        return elements  # Default: no filtering

    def _deduplicate(self, elements):
        """Hook method for deduplication (optional, can be overridden)."""
        return elements  # Default: no deduplication

    def _order_elements(self, elements):
        """Hook method for ordering elements (optional, can be overridden)."""
        return elements  # Default: no specific ordering
