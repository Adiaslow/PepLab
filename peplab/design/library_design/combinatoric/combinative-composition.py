from abc import ABC, abstractmethod
from composition import Composition

class CombinatorialComposition(Composition):
    @abstractmethod
    def generate_composition(self):
        raise NotImplementedError("This method must be implemented in a subclass")
