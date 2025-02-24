# peplab/frontend/src/core/types/design_types.py
"""
This module is responsible for defining the types of designs.

Classes:
    DesignType: The type of design to perform.
"""

# External imports
from enum import Enum
from typing import List


class DesignType(Enum):
    """The type of design to perform."""

    COMBINATORIC = "combinatoric"
    GENERATIVE = "generative"
    GENETIC = "genetic"
    MCMC = "mcmc"
    FRACTAL = "fractal"
    RANDOM = "random"

    @classmethod
    def get_all_design_types(cls) -> List[str]:
        """
        This method returns all the design types.
        """
        return [design_type.value for design_type in cls]
