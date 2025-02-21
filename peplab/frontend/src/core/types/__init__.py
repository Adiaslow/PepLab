# peplab/frontend/src/core/types/__init__.py
"""
This module contains the types used in the frontend application.
"""

# External imports
from typing import Any, List

# Internal imports
from peplab.frontend.src.core.types.analysis_types import AnalysisType
from peplab.frontend.src.core.types.design_types import DesignType

__all__: List[Any] = [
    "AnalysisType",
    "DesignType",
]
