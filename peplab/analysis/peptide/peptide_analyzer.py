from typing import Any, Dict, List
from rdkit import Chem
import pandas as pd

from ...core.base.property_store import PropertyConfig
from ...core.molecule.peptide import PeptideInfo
from ...utils.smiles_generator import SMILESGenerator
from ..property.property_calculator import PropertyCalculator
from ..property.property_calculator_factory import PropertyCalculatorFactory

class PeptideAnalyzer:
    """Analyzes peptide structures and computes properties."""

    def __init__(self, property_config: PropertyConfig):
        """Initializes analyzer with property configuration.

        Args:
            property_config: Configuration specifying which properties to compute.
        """
        self.property_config = property_config
        self.smiles_generator = SMILESGenerator()
        self.calculators = {}

        # Initialize calculators based on config
        for prop_name, enabled in property_config.__dict__.items():
            if enabled:
                self.calculators[prop_name] = (
                    PropertyCalculatorFactory.create(prop_name)
                )

    def analyze_peptide(
        self,
        peptide: PeptideInfo
    ) -> Dict[str, Any]:
        """Analyzes a single peptide.

        Args:
            peptide: PeptideInfo object to analyze.

        Returns:
            Dictionary containing analysis results.
        """
        results = {
            'sequence': peptide.sequence,
            'type': peptide.peptide_type
        }

        # Generate SMILES
        print(peptide)
        smiles, success = self.smiles_generator.generate(peptide.graph)
        results['smiles'] = smiles if success else None

        if success and smiles:
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                # Calculate properties
                for prop_name, calculator in self.calculators.items():
                    try:
                        results[prop_name] = calculator.calculate(mol)
                    except Exception as e:
                        print(f"Error calculating {prop_name}: {str(e)}")
                        results[prop_name] = None
            else:
                # Set all properties to None if molecule conversion failed
                for prop_name in self.calculators:
                    results[prop_name] = None

        return results

    def analyze_library(
        self,
        peptides: List[PeptideInfo]
    ) -> pd.DataFrame:
        """Analyzes a list of peptides.

        Args:
            peptides: List of PeptideInfo objects to analyze.

        Returns:
            DataFrame containing analysis results.
        """
        results = []
        for peptide in peptides:
            analysis = self.analyze_peptide(peptide)
            results.append(analysis)

        return pd.DataFrame(results)
