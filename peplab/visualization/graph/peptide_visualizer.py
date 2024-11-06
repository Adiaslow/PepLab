from typing import Tuple, Optional, Union
from pathlib import Path

from ..rdkit_visualizer import RDKitVisualizer
from ...core.molecule.peptide import PeptideInfo

class PeptideVisualizer:
    """Combines different visualization methods for peptides."""

    def __init__(
        self,
        size: Tuple[int, int] = (400, 400),
        fmt: str = 'svg',
        **kwargs  # Accept additional kwargs for flexibility
    ):
        """Initializes the visualizer.

        Args:
            size: Size for molecular depictions in pixels.
            fmt: Output format for molecular depictions.
            **kwargs: Additional configuration options.
        """
        self.size = size
        self.fmt = fmt
        # Store graph size from kwargs or use default
        self.graph_size = kwargs.get('graph_size', (10, 10))
        self.visualizer = RDKitVisualizer()

    def visualize_peptide(
        self,
        peptide_info: PeptideInfo,
        save_path: Optional[Union[str, Path]] = None
    ) -> None:
        """Visualizes a peptide structure.

        Args:
            peptide_info: PeptideInfo object containing peptide data.
            save_path: Optional path to save the visualization.
        """
        # Create title
        title = (f"{'->'.join(peptide_info.sequence)} "
                f"({peptide_info.peptide_type})")

        # Generate visualization
        img_data, _ = self.visualizer.create_2d_depiction(
            peptide_info.graph,
            size=self.size,
            title=title,
            fmt=self.fmt
        )

        # Save if path provided
        if save_path:
            mode = 'wb' if self.fmt == 'png' else 'w'
            save_path = Path(save_path)
            mol_path = save_path.with_stem(f"{save_path.stem}_mol")
            mol_path = mol_path.with_suffix(f".{self.fmt}")

            with open(mol_path, mode) as f:
                f.write(img_data)
