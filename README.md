
# PepLab: Peptide Library Generator

**PepLab** is a modular and extensible tool for generating peptide libraries using various combinatorial and group-theoretic composition strategies. This tool enables researchers to efficiently create libraries of peptide sequences, export results to CSV files, and integrate them into downstream workflows.

---

## **Features**

- **Modular Composition Framework**:
  - Combinatorial strategies: Combinations, permutations, Cartesian products, and k-fold Cartesian products.
  - Group-theoretic strategies: Cyclic and dihedral permutations.

- **Flexible Export Options**:
  - Export generated libraries to CSV files for easy sharing and downstream analysis.

- **Extendable Design**:
  - Built using the Template and Strategy patterns to facilitate the addition of new composition methods.

- **Testing and Validation**:
  - Comprehensive unit tests for each component ensure reliability and accuracy.

---

## **Workflow**

### **1. Define Input Data**
Start with a list of elements (e.g., amino acids or fragments) to combine into sequences.
```python
amino_acids = ['A', 'R', 'N', 'D']
```

---

### **2. Choose a Composition Strategy**
Select a strategy for combining the input:
```python
# Combinative Strategies
# - Combinations: Select subsets of a specific size (order doesn’t matter).
# - Permutations: Generate all possible orderings (order matters).
# - Cartesian Products: Combine elements from multiple sets.
# - K-fold Cartesian Products: Repeated combinations of a single set.

# Group-Theoretic Strategies
# - Cyclic Permutations: Rotate elements within a list.
# - Dihedral Permutations: Include both rotations and reflections.
```

---

### **3. Use the Composer**
Pass the input data and selected strategy to the `Composer` class to generate the library:
```python
from combinative_composition import CombinationComposition
from composer import Composer

# Initialize the composition strategy and composer
combination_strategy = CombinationComposition()
composer = Composer(combination_strategy)

# Generate the peptide library
peptide_library = composer.generate_library(amino_acids, r=3)
```

---

### **4. Export the Results**
Save the generated library to a CSV file using the export functionality:
```python
composer.export_to_csv(peptide_library, "peptide_library.csv")
```

