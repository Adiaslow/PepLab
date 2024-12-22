# peplab/chemistry/reactivity/reactive_group_patterns/__init__.py

from peplab.chemistry.reactivity.reactive_group_pattern import ReactiveGroupPattern

# Import all patterns from submodules
from .carbon_patterns import *
from .nitrogen_patterns import *
from .oxygen_patterns import *
from .sulfur_patterns import *

# These are already defined in their respective modules but let's make them
# available at the package level
__all__ = [
    # Carbon patterns
    "METHYL_PATTERN",
    "METHYLENE_PATTERN",
    "METHINE_PATTERN",
    "AROMATIC_CARBON_PATTERN",
    "ALKENE_PATTERN",
    "ALKYNE_PATTERN",
    # Nitrogen patterns
    "PRIMARY_AMINE_PATTERN",
    "SECONDARY_AMINE_PATTERN",
    "TERTIARY_AMINE_PATTERN",
    "QUATERNARY_AMMONIUM_PATTERN",
    "IMINE_PATTERN",
    "AMIDE_PATTERN",
    "NITRO_PATTERN",
    # Oxygen patterns
    "HYDROXYL_PATTERN",
    "CARBONYL_PATTERN",
    "CARBOXYL_PATTERN",
    "ESTER_PATTERN",
    "ETHER_PATTERN",
    # Sulfur patterns
    "THIOL_PATTERN",
    "SULFIDE_PATTERN",
    "DISULFIDE_PATTERN",
    "SULFOXIDE_PATTERN",
    "SULFONE_PATTERN"
]
