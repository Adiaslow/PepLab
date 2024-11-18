from abc import ABC, abstractmethod

class GroupTheoreticComposition(ABC):
    """
    Abstract base class for group-theoretic compositions.
    """

    def cyclic_permutations(self, items):
        """
        Generates all cyclic permutations of a list of items.
        """
        n = len(items)
        return [items[i:] + items[:i] for i in range(n)]

    @abstractmethod
    def generate_composition(self, items):
        """
        Abstract method to generate compositions.
        Must be implemented by subclasses.
        """
        pass


class CyclicComposition(GroupTheoreticComposition):
    """
    Generates cyclic permutations of a list of items.
    """

    def generate_composition(self, items):
        return self.cyclic_permutations(items)

class CyclicComposition(GroupTheoreticComposition):
    """
    Generates cyclic permutations of a list of items.
    """

    def generate_composition(self, items):
        return self.cyclic_permutations(items)
