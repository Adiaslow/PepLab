from dataclasses import dataclass
from typing import Dict

from ..molecule.residue import ResidueInfo
from ..molecule.position import PositionInfo


@dataclass
class LibraryInfo:
    """Complete information about a peptide library.

    Attributes:
        positions: Dictionary mapping position keys to PositionInfo objects.
        properties: Dictionary of property calculation flags.
    """
    positions: Dict[str, PositionInfo]
    properties: Dict[str, bool]

    @classmethod
    def from_dict(cls, data: Dict) -> 'LibraryInfo':
        """Creates a LibraryInfo instance from a dictionary.

        Args:
            data: Dictionary containing library information.

        Returns:
            A new LibraryInfo instance.
        """
        positions = {}
        for pos_key, pos_data in data.items():
            if pos_key not in ('properties', 'temp'):
                residues = {}
                for res_key, res_data in pos_data.items():
                    if res_key != 'temp':
                        residues[res_key] = ResidueInfo(
                            name=res_data['id'],
                            smiles=res_data['smiles'],
                            nucleophile=res_data['nuc'],
                            electrophile=res_data['elec']
                        )
                positions[pos_key] = PositionInfo(
                    residues=residues,
                    temperature=pos_data.get('temp', 20.0)
                )

        properties = data.get('properties', {
            'alogp': True,
            'exact_mass': True,
            'rotatable_bonds': True,
            'hbd_count': True,
            'hba_count': True
        })

        return cls(positions=positions, properties=properties)

    def to_dict(self) -> Dict:
        """Converts the LibraryInfo instance to a dictionary.

        Returns:
            Dictionary representation of the library info.
        """
        result = {}

        for pos_key, pos_info in self.positions.items():
            result[pos_key] = {}
            for res_key, res_info in pos_info.residues.items():
                result[pos_key][res_key] = {
                    'id': res_info.name,
                    'smiles': res_info.smiles,
                    'nuc': res_info.nucleophile,
                    'elec': res_info.electrophile
                }
            result[pos_key]['temp'] = pos_info.temperature

        result['properties'] = self.properties
        return result
