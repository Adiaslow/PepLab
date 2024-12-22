# peplab/chemistry/reactivity/reactive_group_pattern_registry.py

from typing import Dict, List, Optional
from peplab.chemistry.reactivity.reactive_group_pattern import ReactiveGroupPattern
from peplab.chemistry.reactivity.reactive_group_patterns import (
    carbon_patterns,
    nitrogen_patterns,
    oxygen_patterns,
    sulfur_patterns
)

class ReactiveGroupPatternRegistry:
    """A registry for reactive group patterns."""
    _instance = None
    _initialized = False

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance

    def __init__(self):
        if not self._initialized:
            # Initialize separate dictionaries for each element
            self.patterns_by_element: Dict[str, List[ReactiveGroupPattern]] = {
                'C': [],
                'N': [],
                'O': [],
                'S': []
            }
            self._load_patterns()
            self.__class__._initialized = True

    def _load_patterns(self) -> None:
        """Load patterns into element-specific lists."""
        pattern_modules = [
            carbon_patterns,
            nitrogen_patterns,
            oxygen_patterns,
            sulfur_patterns
        ]
        for module in pattern_modules:
            for name, obj in vars(module).items():
                if isinstance(obj, ReactiveGroupPattern):
                    # Get the element from the first node criteria
                    element = obj.node_criteria[0].get('element')
                    if element:
                        self.patterns_by_element[element].append(obj)

    def get_patterns_for_element(self, element: str) -> List[ReactiveGroupPattern]:
        """Get all patterns for a specific element."""
        return self.patterns_by_element.get(element, [])

    def get_all_patterns(self) -> List[ReactiveGroupPattern]:
        """Get all registered patterns."""
        all_patterns = []
        for patterns in self.patterns_by_element.values():
            all_patterns.extend(patterns)
        return all_patterns

# Create and expose the singleton instance
pattern_registry = ReactiveGroupPatternRegistry()
