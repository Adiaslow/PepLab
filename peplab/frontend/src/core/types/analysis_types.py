# peplab/frontend/src/core/types/analysis_types.py
"""
This module contains the types of analysis.
"""

# External imports
from enum import Enum
from typing import List


class AnalysisType(Enum):
    """The type of analysis."""

    CHEMINFORMATIC = "cheminformatic"
    DATA_ANALYSIS = "data_analysis"
    MACHINE_LEARNING = "machine_learning"
    SIMULATION = "simulation"

    @classmethod
    def get_all_analysis_types(cls) -> List[str]:
        """
        This method returns all the analysis types.
        """
        return [analysis_type.value for analysis_type in cls]
