# Python Codebase Summary

Generated on: 2024-12-20 16:35:19

## Summary Statistics
- Total Python files: 105
- Total functions: 122

---



## Directory: peplab

### __init__.py
**File Statistics:**
- Total lines: 61
- Non-empty lines: 53
- Number of functions: 0
**File Description:**
PepLab: A Python package for peptide library design and analysis.

This package provides tools for generating, optimizing, and analyzing
in silico peptide libraries with a focus on chemical reactivity and
structural properties.
---

### setup.py
**File Statistics:**
- Total lines: 26
- Non-empty lines: 24
- Number of functions: 0
---


## Directory: peplab/visualization

### __init__.py
**File Statistics:**
- Total lines: 9
- Non-empty lines: 6
- Number of functions: 0
---

### cpk_colors.py
**File Statistics:**
- Total lines: 119
- Non-empty lines: 117
- Number of functions: 1
**Functions:**
```python
def get_color
```
---

### rdkit_visualizer.py
**File Statistics:**
- Total lines: 47
- Non-empty lines: 39
- Number of functions: 2
**Functions:**
```python
def create_2d_depiction
def _configure_drawer
```
---


## Directory: peplab/visualization/graph

### __init__.py
**File Statistics:**
- Total lines: 8
- Non-empty lines: 6
- Number of functions: 0
---

### graph_visualizer.py
**File Statistics:**
- Total lines: 320
- Non-empty lines: 295
- Number of functions: 4
**Functions:**
```python
def create_graph_plot
def _draw_nodes
def _draw_edges
def _add_legend
```
---

### peptide_visualizer.py
**File Statistics:**
- Total lines: 61
- Non-empty lines: 51
- Number of functions: 2
**Functions:**
```python
def __init__
def visualize_peptide
```
---


## Directory: peplab/design

### __init__.py
**File Statistics:**
- Total lines: 2
- Non-empty lines: 0
- Number of functions: 0
---


## Directory: peplab/design/library_design

### __init__.py
**File Statistics:**
- Total lines: 2
- Non-empty lines: 0
- Number of functions: 0
---

### composer.py
**File Statistics:**
- Total lines: 42
- Non-empty lines: 33
- Number of functions: 3
**Functions:**
```python
def __init__
def generate_library
def export_to_csv
```
---


## Directory: peplab/design/library_design/reaction

### __init__.py
**File Statistics:**
- Total lines: 1
- Non-empty lines: 0
- Number of functions: 0
---

### reaction_planner.py
**File Statistics:**
- Total lines: 1
- Non-empty lines: 0
- Number of functions: 0
---


## Directory: peplab/design/library_design/generative/genetic

### atom_based.py
**File Statistics:**
- Total lines: 1
- Non-empty lines: 0
- Number of functions: 0
---

### building_block_based.py
**File Statistics:**
- Total lines: 1
- Non-empty lines: 0
- Number of functions: 0
---

### residue_based.py
**File Statistics:**
- Total lines: 1
- Non-empty lines: 0
- Number of functions: 0
---


## Directory: peplab/design/library_design/generative/markov_chain_monte_carlo

### atom_based.py
**File Statistics:**
- Total lines: 1
- Non-empty lines: 0
- Number of functions: 0
---

### building_block_based.py
**File Statistics:**
- Total lines: 1
- Non-empty lines: 0
- Number of functions: 0
---

### residue_based.py
**File Statistics:**
- Total lines: 1
- Non-empty lines: 0
- Number of functions: 0
---


## Directory: peplab/design/library_design/generative/neural

### atom_based.py
**File Statistics:**
- Total lines: 1
- Non-empty lines: 0
- Number of functions: 0
---

### building_block_based.py
**File Statistics:**
- Total lines: 1
- Non-empty lines: 0
- Number of functions: 0
---

### residue_based.py
**File Statistics:**
- Total lines: 1
- Non-empty lines: 0
- Number of functions: 0
---


## Directory: peplab/design/library_design/combinatoric

### __init__.py
**File Statistics:**
- Total lines: 1
- Non-empty lines: 0
- Number of functions: 0
---

### combinative-composition.py
**File Statistics:**
- Total lines: 91
- Non-empty lines: 68
- Number of functions: 5
**Functions:**
```python
def generate_composition
def generate_composition
def generate_composition
def generate_composition
def generate_composition
```
---

### combinative.py
**File Statistics:**
- Total lines: 27
- Non-empty lines: 22
- Number of functions: 1
**Functions:**
```python
def generate_combinations
```
---

### combinative_composition.py
**File Statistics:**
- Total lines: 64
- Non-empty lines: 47
- Number of functions: 5
**Functions:**
```python
def _generate_elements
def _generate_elements
def _generate_elements
def _generate_elements
def _generate_elements
```
---

### permutative.py
**File Statistics:**
- Total lines: 25
- Non-empty lines: 21
- Number of functions: 1
**Functions:**
```python
def generate_permutations
```
---


## Directory: peplab/design/library_design/group_thoeretic

### cyclic_permutative.py
**File Statistics:**
- Total lines: 63
- Non-empty lines: 52
- Number of functions: 2
**Functions:**
```python
def generate_cyclic_permutations
def _deduplicate
```
---

### dihedral_permutative.py
**File Statistics:**
- Total lines: 32
- Non-empty lines: 25
- Number of functions: 1
**Functions:**
```python
def generate_dihedral_permutations
```
---

### permutative.py
**File Statistics:**
- Total lines: 2
- Non-empty lines: 0
- Number of functions: 0
---


## Directory: peplab/core

### __init__.py
**File Statistics:**
- Total lines: 31
- Non-empty lines: 27
- Number of functions: 0
---


## Directory: peplab/core/reaction

### __init__.py
**File Statistics:**
- Total lines: 9
- Non-empty lines: 7
- Number of functions: 0
---

### reaction_conditions.py
**File Statistics:**
- Total lines: 8
- Non-empty lines: 5
- Number of functions: 1
**Functions:**
```python
def __init__
```
---

### reaction_mechanism.py
**File Statistics:**
- Total lines: 7
- Non-empty lines: 5
- Number of functions: 2
**Functions:**
```python
def __init__
def get_type
```
---

### reaction_pathway.py
**File Statistics:**
- Total lines: 4
- Non-empty lines: 3
- Number of functions: 1
**Functions:**
```python
def __init__
```
---

### reactive_pattern.py
**File Statistics:**
- Total lines: 14
- Non-empty lines: 11
- Number of functions: 0
---

### reactive_site.py
**File Statistics:**
- Total lines: 16
- Non-empty lines: 12
- Number of functions: 1
**Functions:**
```python
def __hash__
```
---

### reactive_type.py
**File Statistics:**
- Total lines: 9
- Non-empty lines: 7
- Number of functions: 0
---


## Directory: peplab/core/graph

### __init__.py
**File Statistics:**
- Total lines: 4
- Non-empty lines: 2
- Number of functions: 0
---

### molecule_graph.py
**File Statistics:**
- Total lines: 264
- Non-empty lines: 226
- Number of functions: 9
**Functions:**
```python
def __init__
def from_smiles
def get_neighbors
def find_reactive_sites
def _find_amine_patterns
def _find_carboxyl_patterns
def _find_azide_patterns
def _find_alkyne_patterns
def to_dict
```
---


## Directory: peplab/core/library

### __init__.py
**File Statistics:**
- Total lines: 10
- Non-empty lines: 8
- Number of functions: 0
---

### library.py
**File Statistics:**
- Total lines: 51
- Non-empty lines: 45
- Number of functions: 2
**Functions:**
```python
def from_dict
def to_dict
```
---

### library_parser.py
**File Statistics:**
- Total lines: 87
- Non-empty lines: 74
- Number of functions: 3
**Functions:**
```python
def parse
def _parse_csv
def _parse_json
```
---

### peptide_builder.py
**File Statistics:**
- Total lines: 898
- Non-empty lines: 760
- Number of functions: 23
**Functions:**
```python
def __init__
def _process_combination
def _progress_tracker
def enumerate_library
def build_linear_peptide
def _form_peptide_bond
def _get_amine_type
def _remove_h_from_nh
def _remove_oh_from_cooh
def _has_click_pairs
def click_cyclize_peptide
def _get_azide_chain
def _get_alkyne_chain
def _is_azide_nitrogen
def _is_alkyne_carbon
def validate_click_chemistry_pair
def _get_connected_atoms
def _get_azide_nitrogens
def cyclize_peptide
def _reindex_graph
def count_atoms_by_element
def check_ionic_azide
def check_covalent_azide
```
---

### peptide_library_generator.py
**File Statistics:**
- Total lines: 215
- Non-empty lines: 178
- Number of functions: 8
**Functions:**
```python
def __init__
def load_library
def generate_peptides
def analyze_peptides
def visualize_peptides
def _save_intermediates
def _save_analysis_results
def _setup_logger
```
---


## Directory: peplab/core/molecule

### __init__.py
**File Statistics:**
- Total lines: 7
- Non-empty lines: 5
- Number of functions: 0
---

### atom.py
**File Statistics:**
- Total lines: 53
- Non-empty lines: 49
- Number of functions: 2
**Functions:**
```python
def to_dict
def __hash__
```
---

### bond.py
**File Statistics:**
- Total lines: 31
- Non-empty lines: 27
- Number of functions: 2
**Functions:**
```python
def to_dict
def __hash__
```
---

### peptide.py
**File Statistics:**
- Total lines: 12
- Non-empty lines: 10
- Number of functions: 0
---

### position.py
**File Statistics:**
- Total lines: 16
- Non-empty lines: 12
- Number of functions: 0
---

### residue.py
**File Statistics:**
- Total lines: 31
- Non-empty lines: 27
- Number of functions: 2
**Functions:**
```python
def from_dict
def to_dict
```
---


## Directory: peplab/core/base

### __init__.py
**File Statistics:**
- Total lines: 5
- Non-empty lines: 3
- Number of functions: 0
---

### molecular_entity.py
**File Statistics:**
- Total lines: 20
- Non-empty lines: 14
- Number of functions: 3
**Functions:**
```python
def __init__
def get_type
def __hash__
```
---

### property_store.py
**File Statistics:**
- Total lines: 48
- Non-empty lines: 39
- Number of functions: 5
**Functions:**
```python
def __init__
def set_property
def get_property
def __hash__
def from_dict
```
---


## Directory: peplab/analysis

### __init__.py
**File Statistics:**
- Total lines: 12
- Non-empty lines: 10
- Number of functions: 0
---


## Directory: peplab/analysis/peptide

### __init__.py
**File Statistics:**
- Total lines: 6
- Non-empty lines: 4
- Number of functions: 0
---

### peptide_analyzer.py
**File Statistics:**
- Total lines: 89
- Non-empty lines: 72
- Number of functions: 3
**Functions:**
```python
def __init__
def analyze_peptide
def analyze_library
```
---


## Directory: peplab/analysis/peptide/prediction/bioactivity

### bioactivity_prediction_manager.py
**File Statistics:**
- Total lines: 1
- Non-empty lines: 0
- Number of functions: 0
---


## Directory: peplab/analysis/peptide/prediction/permeability

### permeability_prediction_manager.py
**File Statistics:**
- Total lines: 1
- Non-empty lines: 0
- Number of functions: 0
---


## Directory: peplab/analysis/peptide/structure

### structure_manager.py
**File Statistics:**
- Total lines: 1
- Non-empty lines: 0
- Number of functions: 0
---


## Directory: peplab/analysis/peptide/structure/prediction

### alphafold.py
**File Statistics:**
- Total lines: 1
- Non-empty lines: 0
- Number of functions: 0
---

### chai1.py
**File Statistics:**
- Total lines: 1
- Non-empty lines: 0
- Number of functions: 0
---

### esmfold.py
**File Statistics:**
- Total lines: 1
- Non-empty lines: 0
- Number of functions: 0
---

### rosettafold.py
**File Statistics:**
- Total lines: 1
- Non-empty lines: 0
- Number of functions: 0
---


## Directory: peplab/analysis/peptide/structure/force_field

### amber.py
**File Statistics:**
- Total lines: 1
- Non-empty lines: 0
- Number of functions: 0
---

### mmff.py
**File Statistics:**
- Total lines: 1
- Non-empty lines: 0
- Number of functions: 0
---

### uff.py
**File Statistics:**
- Total lines: 1
- Non-empty lines: 0
- Number of functions: 0
---


## Directory: peplab/analysis/property

### __init__.py
**File Statistics:**
- Total lines: 5
- Non-empty lines: 3
- Number of functions: 0
---

### alogp_calculator.py
**File Statistics:**
- Total lines: 11
- Non-empty lines: 7
- Number of functions: 1
**Functions:**
```python
def calculate
```
---

### exact_mass_calculator.py
**File Statistics:**
- Total lines: 11
- Non-empty lines: 7
- Number of functions: 1
**Functions:**
```python
def calculate
```
---

### hba_count_calculator.py
**File Statistics:**
- Total lines: 11
- Non-empty lines: 7
- Number of functions: 1
**Functions:**
```python
def calculate
```
---

### hbd_count_calculator.py
**File Statistics:**
- Total lines: 11
- Non-empty lines: 7
- Number of functions: 1
**Functions:**
```python
def calculate
```
---

### property_analyzer.py
**File Statistics:**
- Total lines: 75
- Non-empty lines: 61
- Number of functions: 2
**Functions:**
```python
def calculate_statistics
def generate_summary
```
---

### property_calculator.py
**File Statistics:**
- Total lines: 11
- Non-empty lines: 8
- Number of functions: 1
**Functions:**
```python
def calculate
```
---

### property_calculator_factory.py
**File Statistics:**
- Total lines: 37
- Non-empty lines: 29
- Number of functions: 1
**Functions:**
```python
def create
```
---

### reactivity_profile.py
**File Statistics:**
- Total lines: 25
- Non-empty lines: 19
- Number of functions: 3
**Functions:**
```python
def __init__
def get_reactivity_score
def _calculate_rate_constant
```
---

### rotatable_bonds_calculator.py
**File Statistics:**
- Total lines: 11
- Non-empty lines: 7
- Number of functions: 1
**Functions:**
```python
def calculate
```
---

### thermodynamics.py
**File Statistics:**
- Total lines: 14
- Non-empty lines: 11
- Number of functions: 1
**Functions:**
```python
def get_free_energy
```
---


## Directory: peplab/analysis/prediction/bioactivity

### bioactivity_prediction_manager.py
**File Statistics:**
- Total lines: 2
- Non-empty lines: 0
- Number of functions: 0
---


## Directory: peplab/analysis/prediction/permeability

### permeability_prediction_manager.py
**File Statistics:**
- Total lines: 2
- Non-empty lines: 0
- Number of functions: 0
---


## Directory: peplab/analysis/library

### library_analyzer.py
**File Statistics:**
- Total lines: 49
- Non-empty lines: 38
- Number of functions: 2
**Functions:**
```python
def __init__
def analyze
```
---


## Directory: peplab/analysis/docking

### docking_manager.py
**File Statistics:**
- Total lines: 1
- Non-empty lines: 0
- Number of functions: 0
---

### rosetta_dock.py
**File Statistics:**
- Total lines: 2
- Non-empty lines: 0
- Number of functions: 0
---


## Directory: peplab/analysis/docking/chai

### chai_docking.py
**File Statistics:**
- Total lines: 1
- Non-empty lines: 0
- Number of functions: 0
---


## Directory: peplab/analysis/docking/rosetta

### rosetta_dock.py
**File Statistics:**
- Total lines: 1
- Non-empty lines: 0
- Number of functions: 0
---


## Directory: peplab/analysis/structure

### __init__.py
**File Statistics:**
- Total lines: 1
- Non-empty lines: 0
- Number of functions: 0
---

### structure_generator.py
**File Statistics:**
- Total lines: 8
- Non-empty lines: 5
- Number of functions: 2
**Functions:**
```python
def __init__
def generate
```
---


## Directory: peplab/analysis/structure/prediction

### alphafold.py
**File Statistics:**
- Total lines: 2
- Non-empty lines: 0
- Number of functions: 0
---

### chai1.py
**File Statistics:**
- Total lines: 2
- Non-empty lines: 0
- Number of functions: 0
---

### esmfold.py
**File Statistics:**
- Total lines: 2
- Non-empty lines: 0
- Number of functions: 0
---

### rosettafold.py
**File Statistics:**
- Total lines: 2
- Non-empty lines: 0
- Number of functions: 0
---


## Directory: peplab/analysis/structure/force_field

### amber.py
**File Statistics:**
- Total lines: 2
- Non-empty lines: 0
- Number of functions: 0
---

### mmff.py
**File Statistics:**
- Total lines: 2
- Non-empty lines: 0
- Number of functions: 0
---

### uff.py
**File Statistics:**
- Total lines: 2
- Non-empty lines: 0
- Number of functions: 0
---


## Directory: peplab/test

### test_analysis.py
**File Statistics:**
- Total lines: 1
- Non-empty lines: 0
- Number of functions: 0
---

### test_core.py
**File Statistics:**
- Total lines: 1
- Non-empty lines: 0
- Number of functions: 0
---

### test_design.py
**File Statistics:**
- Total lines: 1
- Non-empty lines: 0
- Number of functions: 0
---

### test_peptide_library_generator.py
**File Statistics:**
- Total lines: 50
- Non-empty lines: 39
- Number of functions: 0
---


## Directory: peplab/utils

### __init__.py
**File Statistics:**
- Total lines: 8
- Non-empty lines: 6
- Number of functions: 0
---

### constants.py
**File Statistics:**
- Total lines: 2
- Non-empty lines: 0
- Number of functions: 0
---

### ge- Number of functions: 1
**Functions:**
```python
def main
```
---

### rdkit_utils.py
**File Statistics:**
- Total lines: 188
- Non-empty lines: 153
- Number of functions: 4
**File Description:**
RDKit utility functions for molecular operations.

This module provides helper functions for working with RDKit molecules
and extracting chemical information.
**Functions:**
```python
def get_atom_info
def get_bond_info
def mol_to_graph_dict
def graph_dict_to_mol
```
---

### scrape_cyclic_pepedia.py
**File Statistics:**
- Total lines: 95
- Non-empty lines: 78
- Number of functions: 1
**Functions:**
```python
def scrape_cyclicpepedia
```
---

### smiles_generator.py
**File Statistics:**
- Total lines: 42
- Non-empty lines: 31
- Number of functions: 1
**Functions:**
```python
def generate
```
---


## Directory: peplab/cli

### generate_peptide_library.py
**File Statistics:**
- Total lines: 100
- Non-empty lines: 85
- Number of functions: 4
**Functions:**
```python
def setup_logging
def parse_args
def load_config
def main
```
---

