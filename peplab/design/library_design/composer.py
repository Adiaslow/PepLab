import csv
from composition import Composition


class Composer:
    """
    A class to integrate and use composition strategies for generating libraries.
    """

    def __init__(self, composition_strategy: Composition):
        """
        Initialize with a specific composition strategy.

        Parameters:
        - composition_strategy: An instance of a Composition subclass.
        """
        self.composition_strategy = composition_strategy

    def generate_library(self, *args, **kwargs):
        """
        Generates a library using the composition strategy.

        Returns:
        - List of generated compositions.
        """
        return self.composition_strategy.generate_composition(*args, **kwargs)

    def export_to_csv(self, library, filename):
        """
        Exports the generated library to a CSV file.

        Parameters:
        - library: The library to export (list of lists or strings).
        - filename: The name of the output CSV file.
        """
        with open(filename, mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(["Generated Compositions"])
            for item in library:
                writer.writerow([''.join(item) if isinstance(item, list) else item])
        print(f"Library successfully exported to {filename}")
