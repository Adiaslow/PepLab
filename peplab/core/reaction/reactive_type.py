# peplab/core/reaction/reactive_type.py

"""
Module defining the ReactiveType enumeration.

This module defines the ReactiveType enumeration, which specifies different types of reactive patterns.
"""

from enum import Enum

class ReactiveType(Enum):
    """
    Enumeration of different types of reactive patterns.
    """
    NH2 = 'NH2'
    NH = 'NH'
    COOH = 'COOH'
    AZIDE = 'AZIDE'
    ALKYNE = 'ALKYNE'
    # Add additional reactive types as needed

    def __str__(self):
        return self.value
