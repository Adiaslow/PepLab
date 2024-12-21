# Python Codebase Summary

Generated on: 2024-12-20 16:29:09


## Root Directory

### __init__.py
**File Statistics:**
- Total lines: 31
- Non-empty lines: 27
- Number of functions: 0
---


## Directory: reaction

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


## Directory: graph

### __init__.py
**File Statistics:**
- Total lines: 4
- Non-empty lines: 2
- Number of functions: 0
---

### molecule_graph.py
**File Statistics:**
- Total lines: 266
- Non-empty lines: 227
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


## Directory: library

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
- Total lines: 217
- Non-empty lines: 179
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


## Directory: molecule

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


## Directory: base

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

