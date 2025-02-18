# backend/core/composition/building_block_registry.py

from rdkit import Chem
from rdkit.Chem import Descriptors
from backend.core.properties.property_calculator_factory import PropertyCalculatorFactory

class BuildingBlock:
    """Represents a peptide building block with its sequence and calculated properties."""
    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence
        self.properties = self.calculate_properties()

    def calculate_properties(self):
        """Calculates molecular properties using the registered property calculator."""
        calculator = PropertyCalculatorFactory.get_calculator("rdkit")
        if calculator:
            return calculator.calculate(self)
        return {}

class BuildingBlockRegistry:
    """Stores and manages building block peptides."""
    def __init__(self):
        self.registry = {}

    def add_block(self, name, sequence):
        """Adds a new building block to the registry."""
        if name in self.registry:
            print(f"Building block '{name}' already exists.")
            return
        self.registry[name] = BuildingBlock(name, sequence)

    def get_block(self, name):
        """Retrieves a building block by name."""
        return self.registry.get(name, None)

    def list_blocks(self):
        """Returns all registered building blocks."""
        return {name: block.sequence for name, block in self.registry.items()}

# Example usage:
if __name__ == "__main__":
    registry = BuildingBlockRegistry()
    registry.add_block("Block_A", "ARND")
    registry.add_block("Block_B", "GAVL")

    print("Registered Blocks:", registry.list_blocks())
    print("Block_A Properties:", registry.get_block("Block_A").properties)
