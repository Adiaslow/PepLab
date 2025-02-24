# frontend/src/infrastructure/states/optimization_types.py
"""
This module defines the OptimizationTypes enum, which is responsible for
managing the types of optimizations available in the application.
"""

from enum import Enum


class OptimizationType(Enum):
    """
    Enum for optimization types.
    """

    GENETIC = "genetic"
    MCMC = "mcmc"
    MACHINE_LEARNING = "machine_learning"
