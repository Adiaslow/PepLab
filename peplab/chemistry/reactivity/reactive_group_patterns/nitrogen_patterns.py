# peplab/chemistry/reactivity/reactive_group_patterns/nitrogen_patterns.py

"""
Pattern definitions for nitrogen-containing reactive groups.
"""

from .. import ReactiveGroupPattern

# Primary Amine Pattern: NH2
PRIMARY_AMINE_PATTERN = ReactiveGroupPattern(
    name='primary_amine',
    node_criteria=[
        {'element': 'N', 'neighbors': 3, 'hydrogen_count': 2}
    ],
    edge_criteria=[
        {'bond_type': 'SINGLE', 'from_element': 'N', 'to_element': 'H'},
        {'bond_type': 'SINGLE', 'from_element': 'N', 'to_element': 'H'}
    ]
)

# Secondary Amine Pattern: NH
SECONDARY_AMINE_PATTERN = ReactiveGroupPattern(
    name='secondary_amine',
    node_criteria=[
        {'element': 'N', 'neighbors': 3, 'hydrogen_count': 1}
    ],
    edge_criteria=[
        {'bond_type': 'SINGLE', 'from_element': 'N', 'to_element': 'H'}
    ]
)

# Tertiary Amine Pattern: N
TERTIARY_AMINE_PATTERN = ReactiveGroupPattern(
    name='tertiary_amine',
    node_criteria=[
        {'element': 'N', 'neighbors': 3, 'hydrogen_count': 0}
    ],
    edge_criteria=[]
)

# Quaternary Ammonium Pattern: N+
QUATERNARY_AMMONIUM_PATTERN = ReactiveGroupPattern(
    name='quaternary_ammonium',
    node_criteria=[
        {'element': 'N', 'neighbors': 4, 'charge': 1}
    ],
    edge_criteria=[]
)

# Imine Pattern: C=N
IMINE_PATTERN = ReactiveGroupPattern(
    name='imine',
    node_criteria=[
        {'element': 'N', 'neighbors': 2, 'double_bond': 'C'}
    ],
    edge_criteria=[
        {'bond_type': 'DOUBLE', 'from_element': 'N', 'to_element': 'C'}
    ]
)

# Amide Pattern: C(=O)-N
AMIDE_PATTERN = ReactiveGroupPattern(
    name='amide',
    node_criteria=[
        {'element': 'N', 'neighbors': 3, 'single_bond': 'C', 'carbonyl_bond': True}
    ],
    edge_criteria=[
        {'bond_type': 'SINGLE', 'from_element': 'N', 'to_element': 'C'},
        {'bond_type': 'DOUBLE', 'from_element': 'C', 'to_element': 'O'}
    ]
)

# Nitro Group Pattern: NO2
NITRO_PATTERN = ReactiveGroupPattern(
    name='nitro',
    node_criteria=[
        {'element': 'N', 'neighbors': 3, 'double_bond': 'O', 'single_bond': 'O'}
    ],
    edge_criteria=[
        {'bond_type': 'DOUBLE', 'from_element': 'N', 'to_element': 'O'},
        {'bond_type': 'SINGLE', 'from_element': 'N', 'to_element': 'O'}
    ]
)
