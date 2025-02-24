# peplab/frontend/src/core/types/workflow_types.py
"""
This module contains the types for the workflow.
"""

# External Imports
from enum import Enum
from typing import List


class WorkflowType(Enum):
    """
    This enum contains the types of workflows.
    """

    DESIGN = "design"
    ANALYSIS = "analysis"
    MODELING = "modeling"
    OPTIMIZATION = "optimization"
    SETTINGS = "settings"

    @classmethod
    def get_all_workflow_types(cls) -> List[str]:
        """
        This method returns all the workflow types.
        """
        return [workflow_type.value for workflow_type in cls]
