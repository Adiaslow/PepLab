# peplab/chemistry/reactivity/reactive_group_patterns/sulfur_patterns.py

"""
Pattern definitions for sulfur-containing reactive groups.
"""

from .. import ReactiveGroupPattern

# Thiol Pattern: SH
THIOL_PATTERN = ReactiveGroupPattern(
    name='thiol',
    node_criteria=[
        {'element': 'S', 'neighbors': 1, 'single_bond': 'H'}
    ],
    edge_criteria=[
        {'bond_type': 'SINGLE', 'from_element': 'S', 'to_element': 'H'}
    ]
)

# Sulfide Pattern: RSR
SULFIDE_PATTERN = ReactiveGroupPattern(
    name='sulfide',
    node_criteria=[
        {'element': 'S', 'neighbors': 2}
    ],
    edge_criteria=[
        {'bond_type': 'SINGLE', 'from_element': 'S', 'to_element': 'R'},
        {'bond_type': 'SINGLE', 'from_element': 'S', 'to_element': 'R'}
    ]
)

# Disulfide Pattern: RSSR
DISULFIDE_PATTERN = ReactiveGroupPattern(
    name='disulfide',
    node_criteria=[
        {'element': 'S', 'neighbors': 2, 'single_bond': 'S'}
    ],
    edge_criteria=[
        {'bond_type': 'SINGLE', 'from_element': 'S', 'to_element': 'S'}
    ]
)

# Sulfoxide Pattern: S=O
SULFOXIDE_PATTERN = ReactiveGroupPattern(
    name='sulfoxide',
    node_criteria=[
        {'element': 'S', 'neighbors': 3, 'double_bond': 'O'}
    ],
    edge_criteria=[
        {'bond_type': 'DOUBLE', 'from_element': 'S', 'to_element': 'O'}
    ]
)

# Sulfone Pattern: SO2
SULFONE_PATTERN = ReactiveGroupPattern(
    name='sulfone',
    node_criteria=[
        {'element': 'S', 'neighbors': 4, 'double_bond': 'O', 'single_bond': 'O'}
    ],
    edge_criteria=[
        {'bond_type': 'DOUBLE', 'from_element': 'S', 'to_element': 'O'},
        {'bond_type': 'SINGLE', 'from_element': 'S', 'to_element': 'O'}
    ]
)
