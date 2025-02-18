from pathlib import Path
import logging
from typing import Dict, Optional, Union, List
import json

from .library_parser import LibraryParser
from .peptide_builder import PeptideBuilder
from ...visualization.graph.peptide_visualizer import PeptideVisualizer
from ...analysis.library.library_analyzer import LibraryAnalysis
from ..molecule.peptide import PeptideInfo
from ...analysis.property.property_analyzer import PropertyAnalysis

class PeptideLibraryGenerator:
    """Main interface for peptide library generation and analysis."""

    def __init__(
        self,
        input_path: Optional[Union[str, Path]] = None,
        output_dir: Optional[Union[str, Path]] = None,
        property_config: Optional[Dict[str, bool]] = None,
        visualization_config: Optional[Dict] = None,
        logger: Optional[logging.Logger] = None
    ):
        """Initializes the peptide library generator.

        Args:
            input_path: Path to input library file (CSV or JSON).
            output_dir: Directory for output files.
            property_config: Configuration for property calculations.
            visualization_config: Configuration for visualizations.
            logger: Optional logger instance.
        """
        self.logger = logger or self._setup_logger()

        self.output_dir = Path(output_dir) if output_dir else Path.cwd()
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Initialize components
        self.parser = LibraryParser()
        self.builder = PeptideBuilder()
        self.visualizer = PeptideVisualizer(**(visualization_config or {}))
        self.analyzer = LibraryAnalysis(property_config)

        # Load library if path provided
        self.library_info = None
        if input_path:
            self.load_library(input_path)

    def load_library(self, input_path: Union[str, Path]) -> None:
        """Loads library from file.

        Args:
            input_path: Path to input library file.
        """
        try:
            self.library_info = self.parser.parse(input_path)
            self.logger.info(f"Successfully loaded library from {input_path}")
        except Exception as e:
            self.logger.error(f"Error loading library: {str(e)}")
            raise

    def generate_peptides(
        self,
        cyclize: bool = False,
        save_intermediates: bool = False
    ) -> List[PeptideInfo]:
        """Generates peptides from loaded library.

        Args:
            cyclize: Whether to generate cyclic peptides.
            save_intermediates: Whether to save intermediate structures.

        Returns:
            List of PeptideInfo objects.
        """
        if not self.library_info:
            raise ValueError("No library loaded")

        try:
            peptides = self.builder.enumerate_library(
                self.library_info,
                cyclize=cyclize
            )

            self.logger.info(
                f"Generated {len(peptides)} "
                f"{'cyclic' if cyclize else 'linear'} peptides"
            )

            if save_intermediates:
                self._save_intermediates(peptides)

            return peptides

        except Exception as e:
            self.logger.error(f"Error generating peptides: {str(e)}")
            raise

    def analyze_peptides(
        self,
        peptides: List[PeptideInfo],
        save_results: bool = True
    ) -> Dict:
        """Analyzes generated peptides.

        Args:
            peptides: List of PeptideInfo objects.
            save_results: Whether to save analysis results.

        Returns:
            Dictionary containing analysis results.
        """
        try:
            # Run analysis
            df = self.analyzer.analyze(peptides)
            stats = PropertyAnalysis.calculate_statistics(df)
            summary = PropertyAnalysis.generate_summary(df)

            results = {
                'dataframe': df,
                'statistics': stats,
                'summary': summary
            }

            if save_results:
                self._save_analysis_results(results)

            self.logger.info("Analysis completed successfully")
            return results

        except Exception as e:
            self.logger.error(f"Error analyzing peptides: {str(e)}")
            raise

    def visualize_peptides(
        self,
        peptides: List[PeptideInfo],
        save_images: bool = True
    ) -> None:
        """Generates visualizations for peptides.

        Args:
            peptides: List of PeptideInfo objects.
            save_images: Whether to save visualization files.
        """
        try:
            for i, peptide in enumerate(peptides):
                # Create visualization filename
                base_name = f"peptide_{i+1}"
                if save_images:
                    save_path = self.output_dir / base_name
                else:
                    save_path = None

                # Generate visualizations
                self.visualizer.visualize_peptide(peptide, save_path)

            self.logger.info(
                f"Generated visualizations for {len(peptides)} peptides"
            )

        except Exception as e:
            self.logger.error(f"Error generating visualizations: {str(e)}")
            raise

    def _save_intermediates(self, peptides: List[PeptideInfo]) -> None:
        """Saves intermediate peptide structures."""
        intermediates_dir = self.output_dir / "intermediates"
        intermediates_dir.mkdir(exist_ok=True)

        for i, peptide in enumerate(peptides):
            peptide_data = {
                'sequence': peptide.sequence,
                'type': peptide.peptide_type,
                'graph': peptide.graph
            }

            output_path = intermediates_dir / f"peptide_{i+1}.json"
            with open(output_path, 'w') as f:
                json.dump(peptide_data, f, indent=2)

    def _save_analysis_results(self, results: Dict) -> None:
        """Saves analysis results to files."""
        # Save DataFrame
        results['dataframe'].to_csv(
            self.output_dir / "analysis_results.csv",
            index=False
        )

        # Save statistics
        with open(self.output_dir / "statistics.json", 'w') as f:
            json.dump(results['statistics'], f, indent=2)

        # Save summary
        with open(self.output_dir / "summary.txt", 'w') as f:
            f.write(results['summary'])

    @staticmethod
    def _setup_logger() -> logging.Logger:
        """Sets up logging configuration."""
        logger = logging.getLogger('PeptideLibraryGenerator')
        logger.setLevel(logging.WARNING)

        handler = logging.StreamHandler()
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )
        handler.setFormatter(formatter)
        logger.addHandler(handler)

        return logger
