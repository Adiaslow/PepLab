# frontend/src/core/types/modeling_types.py
"""
This module defines the ModelingTypes enum, which is responsible for
managing the types of modeling available in the application.
"""

# External Imports
from enum import Enum


class ModelingType(Enum):
    """
    Enum for modeling types.
    """

    DOCKING = "docking"
    FORCE_FIELD = "force_field"
    MACHINE_LEARNING = "machine_learning"
    MOLECULAR_DYNAMICS = "molecular_dynamics"
