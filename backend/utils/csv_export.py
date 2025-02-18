import csv
from typing import List


def export_to_csv(library: List[List[str]], filename: str) -> None:
    """
    Exports a given library to a CSV file.

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
