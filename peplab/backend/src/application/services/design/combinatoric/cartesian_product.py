# peplab/backend/src/application/services/design/combinatoric/cartesian_product.py

import itertools
from typing import List, Any
import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from peplab.backend.src.application.interfaces.design.composition import Composition

class CartesianProduct:
    """N-ary Cartesian Product algorithm for combinatorial synthesis."""
    
    @staticmethod
    def generate_composition(items: List[Any], length: int = 1, **kwargs) -> List[List[Any]]:
        """
        Generates Cartesian product of a single set with itself length times.
        """
        # Calculate permutations with replacement (e.g., A/B with r=2 -> AA, AB, BA, BB)
        prod = list(itertools.product(items, repeat=length))
        return [list(p) for p in prod]
