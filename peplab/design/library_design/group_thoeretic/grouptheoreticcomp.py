from abc import ABC, abstractmethod
from composition import Composition

class GroupTheoreticComposition(Composition):
    
    @abstractmethod
    def generate_composition(self):
        raise NotImplementedError("This method must be implemented in a subclass")
        
class CyclicPermutationComposition(GroupTheoreticComposition):
    def generate_composition(self):
        pass
        # Logic for cyclic perms
        
class DihedralPermutationComposition(GroupTheoreticComposition):
    def generate_composition(self):
        pass
        # Logic for dihedral perms