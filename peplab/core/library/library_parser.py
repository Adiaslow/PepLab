# @title Library Parser

"""Parsers for peptide library input files.

This module handles parsing of CSV and JSON input files into
standardized LibraryInfo objects.
"""

import json
import pandas as pd
from typing import Union
from pathlib import Path

from .library import LibraryInfo
from ..molecule.residue import ResidueInfo
from ..molecule.position import PositionInfo


class LibraryParser:
    """Parser for peptide library input files."""

    @staticmethod
    def parse(file_path: Union[str, Path]) -> LibraryInfo:
        """Parses a library file into a LibraryInfo object.

        Args:
            file_path: Path to the input file (CSV or JSON).

        Returns:
            LibraryInfo object containing the parsed data.

        Raises:
            ValueError: If the file format is not supported.
        """
        file_path = Path(file_path)
        if file_path.suffix.lower() == '.csv':
            return LibraryParser._parse_csv(file_path)
        elif file_path.suffix.lower() == '.json':
            return LibraryParser._parse_json(file_path)
        else:
            raise ValueError(
                f"Unsupported file format: {file_path.suffix}. "
                "Please provide either CSV or JSON file."
            )

    @staticmethod
    def _parse_csv(file_path: Path) -> LibraryInfo:
        """Parses a CSV file into a LibraryInfo object.

        Args:
            file_path: Path to the CSV file.

        Returns:
            LibraryInfo object containing the parsed data.
        """
        df = pd.read_csv(file_path)
        positions = {}

        # Group by position
        for pos_idx in df['position_index'].unique():
            pos_data = df[df['position_index'] == pos_idx]
            residues = {}

            # Create residues for this position
            for _, row in pos_data.iterrows():
                residue_key = f"residue_{row['residue_index'] + 1}"
                residues[residue_key] = ResidueInfo(
                    name=row['name'],
                    smiles=row['smiles'],
                    nucleophile=row['nucleophile'],
                    electrophile=row['electrophile']
                )

            position_key = f"position_{pos_idx + 1}"
            positions[position_key] = PositionInfo(residues=residues)

        # Default properties
        properties = {
            'alogp': True,
            'exact_mass': True,
            'rotatable_bonds': True,
            'hbd_count': True,
            'hba_count': True
        }

        return LibraryInfo(positions=positions, properties=properties)

    @staticmethod
    def _parse_json(file_path: Path) -> LibraryInfo:
        """Parses a JSON file into a LibraryInfo object.

        Args:
            file_path: Path to the JSON file.

        Returns:
            LibraryInfo object containing the parsed data.
        """
        with open(file_path, 'r') as f:
            data = json.load(f)

        # Get first library if multiple are present
        library_key = next(iter(data.keys()))
        library_data = data[library_key]

        return LibraryInfo.from_dict(library_data)
