"""
PepLab: A Python package for peptide library design and analysis.

This package provides tools for generating, optimizing, and analyzing
in silico peptide libraries with a focus on chemical reactivity and
structural properties.
"""

from .core import (
    MolecularEntity,
    PropertyStore,
    Atom,
    Bond,
    Residue,
    Peptide,
    MoleculeGraph,
    ReactiveSite,
    ReactionMechanism,
    ReactionPathway,
    PeptideLibrary,
    LibraryParser
)

from .analysis import (
    ThermodynamicParameters,
    ReactivityProfile,
    StructureGenerator
)

from .utils import (
    ReactivityType,
    validate_smiles
)

__version__ = '0.1.0'
__author__ = 'Adam Murray'
__email__ = 'admmurra@ucsc.edu'

__all__ = [
    # Core classes
    'MolecularEntity',
    'PropertyStore',
    'Atom',
    'Bond',
    'Residue',
    'Peptide',
    'MoleculeGraph',
    'ReactiveSite',
    'ReactionMechanism',
    'ReactionPathway',
    'PeptideLibrary',
    'LibraryParser',

    # Analysis classes
    'ThermodynamicParameters',
    'ReactivityProfile',
    'StructureGenerator',

    # Utilities
    'ReactivityType',
    'validate_smiles'
]
