# presentation/app/routes/modeling.py
"""Modeling routes for the PepLab application."""

# Standard Library Imports
from typing import Any, Dict

# External Imports
from flask import Blueprint, render_template

# Internal Imports
from peplab.frontend.src.infrastructure.managers.state_manager import StateManager
from peplab.frontend.src.infrastructure.states.modeling_state import ModelingState
from peplab.frontend.src.core.types.modeling_types import ModelingType

modeling_bp = Blueprint("modeling", __name__)

# Modeling method descriptions
METHOD_DESCRIPTIONS: Dict[str, str] = {
    "docking": "Simulate peptide-protein docking interactions",
    "force_field": "Calculate molecular mechanics using force fields",
    "machine_learning": "Apply ML models for structure prediction",
    "molecular_dynamics": "Simulate peptide dynamics and conformations",
}

# Display names for methods
DISPLAY_NAMES: Dict[str, str] = {
    "docking": "Molecular Docking",
    "force_field": "Force Field",
    "machine_learning": "Machine Learning",
    "molecular_dynamics": "Molecular Dynamics",
}


@modeling_bp.route("/")
def modeling() -> Any:
    """Render the modeling hub page.

    Returns:
        Rendered modeling hub page or coming soon page
    """
    return render_template("coming_soon.html")  # Instead of modeling.html


@modeling_bp.route("/<modeling_type>")
def handle_modeling(modeling_type: str) -> Any:
    """Handle specific modeling type routes.

    Args:
        modeling_type: The type of modeling to handle

    Returns:
        Rendered modeling type page
    """
    try:
        modeling_enum = ModelingType(modeling_type.lower())
        state_manager = StateManager()

        if not isinstance(state_manager.current_state, ModelingState):
            state = ModelingState()
            state.substate = modeling_enum
            state_manager.set_state(state)
        else:
            state_manager.current_state.substate = modeling_enum

        return render_template(
            f"modeling/{modeling_type}.html", modeling_type=modeling_type
        )
    except ValueError:
        return "Invalid modeling type", 404
