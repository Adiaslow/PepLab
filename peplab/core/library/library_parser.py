# peplab/builder/parser.py

import json
import pandas as pd
from typing import Union
from pathlib import Path
from uuid import uuid4
from .library import LibraryInfo
from ..molecule.residue import ResidueInfo
from ..molecule.position import PositionInfo

class LibraryParser:
    """Parser for peptide library input files."""

    @staticmethod
    def parse(file_path: Union[str, Path]) -> LibraryInfo:
        """Parse library file into LibraryInfo object."""
        file_path = Path(file_path)
        if file_path.suffix.lower() == '.csv':
            return LibraryParser._parse_csv(file_path)
        else:
            raise ValueError("Please provide a CSV file")

    @staticmethod
    def _parse_csv(file_path: Path) -> LibraryInfo:
        """Parse CSV file with specified format."""
        try:
            # Read CSV
            df = pd.read_csv(file_path)
            required_columns = [
                'library_index', 'position_index', 'residue_index',
                'name', 'smiles', 'nucleophile', 'electrophile'
            ]

            # Validate columns
            missing_cols = [col for col in required_columns if col not in df.columns]
            if missing_cols:
                raise ValueError(f"Missing required columns: {missing_cols}")

            # Get first library
            library_idx = df['library_index'].iloc[0]
            df = df[df['library_index'] == library_idx]
            positions = {}

            # Group by position
            for pos_idx in df['position_index'].unique():
                pos_data = df[df['position_index'] == pos_idx]
                residues = {}

                # Create residues for this position
                for _, row in pos_data.iterrows():
                    residue_key = f"residue_{row['residue_index'] + 1}"
                    residues[residue_key] = ResidueInfo(
                        id=str(uuid4()),
                        name=row['name'],
                        smiles=row['smiles'],
                        nucleophile=row['nucleophile'],
                        electrophile=row['electrophile']
                    )
                    print(residues[residue_key] )
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

        except Exception as e:
            raise ValueError(f"Error parsing CSV file: {str(e)}")

    @staticmethod
    def _parse_json(file_path: Path) -> LibraryInfo:
        """Parses a JSON file into a LibraryInfo object."""
        with open(file_path, 'r') as f:
            data = json.load(f)
        # Get first library if multiple are present
        library_key = next(iter(data.keys()))
        library_data = data[library_key]
        return LibraryInfo.from_dict(library_data)
