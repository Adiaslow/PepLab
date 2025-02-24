# presentation/app/routes/optimization.py
"""Optimization routes for the PepLab application."""

# Standard Library Imports
from typing import Any, Dict

# External Imports
from flask import Blueprint, render_template

# Internal Imports
from peplab.frontend.src.infrastructure.managers.state_manager import StateManager
from peplab.frontend.src.infrastructure.states.optimization_state import (
    OptimizationState,
)
from peplab.frontend.src.core.types.optimization_types import OptimizationType

optimization_bp = Blueprint("optimization", __name__)

# Optimization method descriptions
METHOD_DESCRIPTIONS: Dict[str, str] = {
    "genetic": "Evolve peptide sequences using genetic algorithms",
    "mcmc": "Sample peptide space using Markov Chain Monte Carlo methods",
    "machine_learning": "Optimize sequences using ML models",
}

# Display names for methods
DISPLAY_NAMES: Dict[str, str] = {
    "genetic": "Genetic Algorithm",
    "mcmc": "Markov Chain Monte Carlo",
    "machine_learning": "Machine Learning",
}


@optimization_bp.route("/")
def optimization() -> Any:
    """Render the optimization hub page.

    Returns:
        Rendered optimization hub page
    """
    state_manager = StateManager()
    if not isinstance(state_manager.current_state, OptimizationState):
        state_manager.set_state(OptimizationState())

    return render_template(
        "optimization/optimization.html",
        method_descriptions=METHOD_DESCRIPTIONS,
        display_names=DISPLAY_NAMES,
    )


@optimization_bp.route("/<optimization_type>")
def handle_optimization(optimization_type: str) -> Any:
    """Handle specific optimization type routes.

    Args:
        optimization_type: The type of optimization to handle

    Returns:
        Rendered optimization type page
    """
    try:
        optimization_enum = OptimizationType(optimization_type.lower())
        state_manager = StateManager()

        if not isinstance(state_manager.current_state, OptimizationState):
            state = OptimizationState()
            state.substate = optimization_enum
            state_manager.set_state(state)
        else:
            state_manager.current_state.substate = optimization_enum

        return render_template(
            f"optimization/{optimization_type}.html",
            optimization_type=optimization_type,
        )
    except ValueError:
        return "Invalid optimization type", 404
