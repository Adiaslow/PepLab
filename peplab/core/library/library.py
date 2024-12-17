# peplab/builder/library.py
from dataclasses import dataclass
from typing import Dict
from ..molecule.residue import ResidueInfo
from ..molecule.position import PositionInfo

@dataclass
class LibraryInfo:
    """Complete information about a peptide library."""
    positions: Dict[str, PositionInfo]
    properties: Dict[str, bool]

    DEFAULT_PROPERTIES = {
        'alogp': True,
        'exact_mass': True,
        'rotatable_bonds': True,
        'hbd_count': True,
        'hba_count': True
    }

    @classmethod
    def from_dict(cls, data: Dict) -> 'LibraryInfo':
        """Creates a LibraryInfo instance from a dictionary."""
        positions = {}
        for pos_key, pos_data in data.items():
            if pos_key not in ('properties', 'temp'):
                residues = {}
                for res_key, res_data in pos_data.items():
                    if res_key != 'temp':
                        # Create ResidueInfo with UUID
                        residues[res_key] = ResidueInfo.from_dict(res_data)
                positions[pos_key] = PositionInfo(
                    residues=residues,
                    temperature=pos_data.get('temp', 20.0)
                )

        # Use provided properties or defaults
        properties = data.get('properties', cls.DEFAULT_PROPERTIES)
        return cls(positions=positions, properties=properties)

    def to_dict(self) -> Dict:
        """Converts the LibraryInfo instance to a dictionary."""
        result = {}
        for pos_key, pos_info in self.positions.items():
            result[pos_key] = {}
            for res_key, res_info in pos_info.residues.items():
                result[pos_key][res_key] = res_info.to_dict()
            result[pos_key]['temp'] = pos_info.temperature
        result['properties'] = self.properties
        return result
