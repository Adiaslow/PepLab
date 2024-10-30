from typing import Tuple, Optional, Union
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, Draw

from peplab.utils.rdkit_utils import graph_dict_to_mol

class RDKitVisualizer:
    """Visualizes molecules using RDKit's 2D drawing capabilities."""

    @staticmethod
    def create_2d_depiction(
        graph_dict: dict,
        size: Tuple[int, int] = (400, 400),
        title: Optional[str] = None,
        fmt: str = 'svg'
    ) -> Tuple[Union[str, bytes], Optional[str]]:
        """Creates a 2D depiction of a molecule."""
        mol = graph_dict_to_mol(graph_dict)
        AllChem.Compute2DCoords(mol)

        if fmt.lower() == 'svg':
            drawer = Draw.rdMolDraw2D.MolDraw2DSVG(size[0], size[1])
        elif fmt.lower() == 'png':
            drawer = Draw.rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
        else:
            raise ValueError(f"Unsupported format: {fmt}. Use 'svg' or 'png'.")

        RDKitVisualizer._configure_drawer(drawer)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()

        return drawer.GetDrawingText(), title

    @staticmethod
    def _configure_drawer(drawer: Draw.rdMolDraw2D.MolDraw2D) -> None:
        """Configures drawing options for a drawer."""
        opts = drawer.drawOptions()
        opts.clearBackground = False
        opts.includeAtomTags = False
        opts.additionalAtomLabelPadding = 0.15
        opts.explicitMethyl = True
        opts.bondLineWidth = 2
        opts.scaleBondWidth = True
        opts.fixedBondLength = 20
        opts.fixedScale = 20
