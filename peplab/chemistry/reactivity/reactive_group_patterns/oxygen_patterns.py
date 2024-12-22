# peplab/chemistry/reactivity/reactive_group_patterns/oxygen_patterns.py

"""
Pattern definitions for oxygen-containing reactive groups.
"""

from .. import ReactiveGroupPattern

# Hydroxyl Pattern: OH
HYDROXYL_PATTERN = ReactiveGroupPattern(
    name='hydroxyl',
    node_criteria=[
        {'element': 'O', 'neighbors': 1, 'single_bond': 'H'}
    ],
    edge_criteria=[
        {'bond_type': 'SINGLE', 'from_element': 'O', 'to_element': 'H'}
    ]
)

# Carbonyl Pattern: C=O
CARBONYL_PATTERN = ReactiveGroupPattern(
    name='carbonyl',
    node_criteria=[
        {'element': 'O', 'neighbors': 1, 'double_bond': 'C'}
    ],
    edge_criteria=[
        {'bond_type': 'DOUBLE', 'from_element': 'C', 'to_element': 'O'}
    ]
)

# Carboxyl Pattern: COOH
CARBOXYL_PATTERN = ReactiveGroupPattern(
    name='carboxyl',
    node_criteria=[
        {'element': 'C', 'neighbors': 3, 'double_bond': 'O'},
        {'element': 'O', 'single_bond': 'H'}
    ],
    edge_criteria=[
        {'bond_type': 'DOUBLE', 'from_element': 'C', 'to_element': 'O'},
        {'bond_type': 'SINGLE', 'from_element': 'O', 'to_element': 'H'}
    ]
)

# Ester Pattern: COOR
ESTER_PATTERN = ReactiveGroupPattern(
    name='ester',
    node_criteria=[
        {'element': 'C', 'neighbors': 3, 'double_bond': 'O', 'single_bond': 'O'}
    ],
    edge_criteria=[
        {'bond_type': 'DOUBLE', 'from_element': 'C', 'to_element': 'O'},
        {'bond_type': 'SINGLE', 'from_element': 'C', 'to_element': 'O'}
    ]
)

# Ether Pattern: ROR
ETHER_PATTERN = ReactiveGroupPattern(
    name='ether',
    node_criteria=[
        {'element': 'O', 'neighbors': 2}
    ],
    edge_criteria=[
        {'bond_type': 'SINGLE', 'from_element': 'O', 'to_element': 'R'},
        {'bond_type': 'SINGLE', 'from_element': 'O', 'to_element': 'R'}
    ]
)
