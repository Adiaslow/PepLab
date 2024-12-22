# peplab/chemistry/reactivity/reactive_group_patterns/carbon_patterns.py

"""
Pattern definitions for carbon-containing reactive groups.
"""

from ..reactive_group_pattern import ReactiveGroupPattern

# Methyl Pattern: CH3
METHYL_PATTERN = ReactiveGroupPattern(
    name='methyl',
    node_criteria=[
        {'element': 'C', 'neighbors': 4, 'hydrogen_count': 3}
    ],
    edge_criteria=[
        {'bond_type': 'SINGLE', 'from_element': 'C', 'to_element': 'H'},
        {'bond_type': 'SINGLE', 'from_element': 'C', 'to_element': 'H'},
        {'bond_type': 'SINGLE', 'from_element': 'C', 'to_element': 'H'}
    ]
)

# Methylene Pattern: CH2
METHYLENE_PATTERN = ReactiveGroupPattern(
    name='methylene',
    node_criteria=[
        {'element': 'C', 'neighbors': 4, 'hydrogen_count': 2}
    ],
    edge_criteria=[
        {'bond_type': 'SINGLE', 'from_element': 'C', 'to_element': 'H'},
        {'bond_type': 'SINGLE', 'from_element': 'C', 'to_element': 'H'}
    ]
)

# Methine Pattern: CH
METHINE_PATTERN = ReactiveGroupPattern(
    name='methine',
    node_criteria=[
        {'element': 'C', 'neighbors': 4, 'hydrogen_count': 1}
    ],
    edge_criteria=[
        {'bond_type': 'SINGLE', 'from_element': 'C', 'to_element': 'H'}
    ]
)

# Aromatic Carbon Pattern: Ar-C
AROMATIC_CARBON_PATTERN = ReactiveGroupPattern(
    name='aromatic_carbon',
    node_criteria=[
        {
            'element': 'C',
            'neighbors': 3,
            'check_aromatic_system': True
        }
    ],
    edge_criteria=[
        {'bond_type': ['SINGLE', 'DOUBLE'], 'count': 3}
    ]
)

# Alkene Pattern: C=C
ALKENE_PATTERN = ReactiveGroupPattern(
    name='alkene',
    node_criteria=[
        {'element': 'C', 'neighbors': 3, 'double_bond': 'C'}
    ],
    edge_criteria=[
        {'bond_type': 'DOUBLE', 'from_element': 'C', 'to_element': 'C'}
    ]
)

# Alkyne Pattern: C≡C
ALKYNE_PATTERN = ReactiveGroupPattern(
    name='alkyne',
    node_criteria=[
        {'element': 'C', 'neighbors': 2, 'triple_bond': 'C'}
    ],
    edge_criteria=[
        {'bond_type': 'TRIPLE', 'from_element': 'C', 'to_element': 'C'}
    ]
)
