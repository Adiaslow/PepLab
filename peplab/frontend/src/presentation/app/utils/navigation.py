"""Navigation utilities for the application."""

from typing import Dict

# Define the page hierarchy
PAGE_HIERARCHY: Dict[str, str] = {
    # Design section
    "design.handle_design": "design.design",  # All design method pages go back to design hub
    "blocks.explore": "design.design",  # Building blocks explorer goes back to design
    "blocks.view_database": "blocks.explore",  # Database view goes back to explorer
    "blocks.view_library": "blocks.explore",  # Library view goes back to explorer
    # Analysis section
    "analysis.handle_analysis": "analysis.analysis",
    # Modeling section
    "modeling.handle_modeling": "modeling.modeling",
    # Optimization section
    "optimization.handle_optimization": "optimization.optimization",
}


def get_parent_route(current_endpoint: str) -> str:
    """Get the parent route for the current endpoint.

    Args:
        current_endpoint: Current route endpoint

    Returns:
        Parent route endpoint or dashboard if not found
    """
    return PAGE_HIERARCHY.get(current_endpoint, "dashboard.index")
