# peplab/frontend/src/core/types/design_types.py
"""
This module is responsible for defining the types of designs.

Classes:
    DesignType: The type of design to perform.
"""

# External imports
from enum import Enum


class DesignType(Enum):
    """The type of design to perform."""

    COMBINATORIC = "combinatoric"
    GENERATIVE = "generative"
    GENETIC = "genetic"
    MCMC = "mcmc"
    RANDOM = "random"
