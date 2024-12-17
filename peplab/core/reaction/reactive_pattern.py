from dataclasses import dataclass
from typing import List

from ...core.molecule.atom import GraphNode
from ...core.molecule.bond import GraphEdge
from .reactive_type import ReactiveType

@dataclass
class ReactivePattern:
    type: ReactiveType
    atoms: List[GraphNode]
    bonds: List[GraphEdge]
    pattern_id: str
