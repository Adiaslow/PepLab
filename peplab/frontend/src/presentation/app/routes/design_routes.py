"""Module for design-related routes."""

# Standard Library Imports
from typing import Any, Dict

# External Imports
from flask import Blueprint, render_template

# Internal Imports
from peplab.frontend.src.infrastructure.states.design_state import DesignState
from peplab.frontend.src.infrastructure.managers.state_manager import StateManager
from peplab.frontend.src.core.types.design_types import DesignType
from peplab.frontend.src.infrastructure.interfaces.state import State

design_bp: Blueprint = Blueprint("design", __name__)


@design_bp.route("/")
def design() -> Any:
    """Render the main design page.

    Returns:
        The main design page.
    """
    # Access the state manager singleton
    state_manager: StateManager = StateManager()

    # Set the state to design state if it's not already
    if not isinstance(state_manager.state, DesignState):
        design_state = DesignState()
        state_manager.set_state(design_state)

    # Get current state for the template
    current_state: State | None = state_manager.get_state()

    # Prepare data for template
    state_data: Dict[str, Any] = {
        "state_type": "design",
        "substate": (
            current_state.substate
            if current_state and hasattr(current_state, "substate")
            else None
        ),
    }

    return render_template("design/design.html", state_data=state_data)


@design_bp.route("/combinatoric")
def combinatoric() -> Any:
    """Render the combinatoric design page."""
    return render_template("design/combinatoric.html")


@design_bp.route("/generative")
def generative() -> Any:
    """Render the generative design page."""
    return render_template("design/generative.html")


@design_bp.route("/genetic")
def genetic() -> Any:
    """Render the genetic design page."""
    return render_template("design/genetic.html")


@design_bp.route("/mcmc")
def mcmc() -> Any:
    """Render the MCMC design page."""
    return render_template("design/mcmc.html")


@design_bp.route("/random")
def random() -> Any:
    """Render the random design page."""
    return render_template("design/random.html")
