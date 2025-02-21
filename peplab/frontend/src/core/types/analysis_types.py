# peplab/frontend/src/core/types/analysis_types.py
"""
This module contains the types of analysis.
"""

# External imports
from enum import Enum


class AnalysisType(Enum):
    """The type of analysis."""

    CHEMINFORMATIC = "cheminformatic"
    DATA_ANALYSIS = "data_analysis"
    MACHINE_LEARNING = "machine_learning"
