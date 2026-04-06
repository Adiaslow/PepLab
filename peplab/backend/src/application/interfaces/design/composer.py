# peplab/backend/src/application/interfaces/design/composer.py
"""
This module defines the Composer class for integrating and using composition strategies.

Classes:
    Composer: A class to integrate and use composition strategies for generating libraries.
"""

# Standard library imports
import csv
import sys
import os
# External imports
from typing import List, Any
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

# Internal imports
from peplab.backend.src.application.interfaces.design.composition import Composition


class Composer:
    """A class to integrate and use composition strategies for generating libraries.

    Attributes:
        composition_strategy: An instance of a Composition subclass.
    """

    def __init__(self, composition_strategy: Composition) -> None:
        """Initialize with a specific composition strategy.

        Args:
            composition_strategy: An instance of a Composition subclass.
        """
        self.composition_strategy: Composition = composition_strategy

    def generate_library(self, *args, **kwargs) -> List[Any]:
        """
        Generates a library using the composition strategy.

        Returns:
        - List of generated compositions.
        """
        return self.composition_strategy.generate_composition(*args, **kwargs)

    def export_to_csv(self, library: List[Any], filename: str) -> None:
        """Exports the generated library to a CSV file.

        Args:
            library: The library to export (list of lists or strings).
            filename: The name of the output CSV file.
        """
        with open(filename, mode="w", newline="") as file:
            writer = csv.writer(file)
            writer.writerow(["Generated Compositions"])
            for item in library:
                writer.writerow(["".join(item) if isinstance(item, list) else item])
        print(f"Library successfully exported to {filename}")
