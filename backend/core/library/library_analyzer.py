from typing import Dict, List, Optional
import pandas as pd

from ...core.base.property_store import PropertyConfig
from ...core.molecule.peptide import PeptideInfo
from ..peptide.peptide_analyzer import PeptideAnalyzer

class LibraryAnalysis:
    """High-level interface for peptide library analysis."""

    def __init__(self, property_config: Optional[Dict[str, bool]] = None):
        """Initializes library analysis.

        Args:
            property_config: Optional dictionary of property calculation flags.
        """
        if property_config is None:
            property_config = {
                'alogp': True,
                'exact_mass': True,
                'rotatable_bonds': True,
                'hbd_count': True,
                'hba_count': True
            }

        self.config = PropertyConfig.from_dict(property_config)
        self.analyzer = PeptideAnalyzer(self.config)

    def analyze(
        self,
        peptides: List[PeptideInfo],
        output_path: Optional[str] = None
    ) -> pd.DataFrame:
        """Analyzes peptide library and optionally saves results.

        Args:
            peptides: List of PeptideInfo objects to analyze.
            output_path: Optional path to save results CSV.

        Returns:
            DataFrame containing analysis results.
        """
        df = self.analyzer.analyze_library(peptides)

        if output_path:
            df.to_csv(output_path, index=False)

        return df
