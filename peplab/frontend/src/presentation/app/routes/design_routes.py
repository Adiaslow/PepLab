# presentation/app/routes/design.py
"""
Module for design-related routes.

This module contains the routes for the design page.

Classes:
    DesignBlueprint: Blueprint for design routes

Functions:
    design: Render the design hub page
    handle_design: Handle specific design type routes
    upload_blocks: Handle building blocks file upload
    explore_blocks: Render the building blocks explorer page
    save_blocks: Handle building blocks export
"""

# Standard Library Imports
from typing import Any, Dict, Set

# External Imports
from flask import Blueprint, render_template, request, flash, redirect, url_for
from werkzeug.utils import secure_filename

# Internal Imports
from peplab.frontend.src.infrastructure.states.design_state import DesignState
from peplab.frontend.src.infrastructure.managers.state_manager import StateManager
from peplab.frontend.src.core.types.design_types import DesignType
from peplab.frontend.src.infrastructure.interfaces.state import State
from peplab.frontend.src.infrastructure.orchestrator import (
    ApplicationOrchestrator,
    ApplicationContext,
)

design_bp = Blueprint("design", __name__)

# Design method descriptions
METHOD_DESCRIPTIONS: Dict[str, str] = {
    "combinatoric": "Generate a peptide library using combinatorial methods",
    "generative": "Generate a peptide library using generative AI models",
    "genetic": "Generate a peptide library through simulated evolution using genetic algorithms",
    "mcmc": "Generate a peptide library using Markov Chain Monte Carlo sampling",
    "fractal": "Generate a peptide library using fractal-based patterns",
    "random": "Generate a peptide library using random sampling",
}

# Display names for methods
DISPLAY_NAMES: Dict[str, str] = {
    "combinatoric": "Combinatorial",
    "generative": "Generative AI",
    "genetic": "Genetic Algorithm",
    "mcmc": "Markov Chain Monte Carlo",
    "fractal": "Fractal-Based",
    "random": "Random Sampling",
}


@design_bp.route("/")
def design() -> Any:
    """Render the design hub page.

    Returns:
        Rendered design hub page
    """
    state_manager = StateManager()
    if not isinstance(state_manager.current_state, DesignState):
        state_manager.set_state(DesignState())

    return render_template(
        "design/design.html",
        method_descriptions=METHOD_DESCRIPTIONS,
        display_names=DISPLAY_NAMES,
    )


@design_bp.route("/<design_type>")
def handle_design(design_type: str) -> Any:
    """Handle specific design type routes.

    Args:
        design_type: The type of design to handle

    Returns:
        Rendered design type page
    """
    try:
        design_enum = DesignType(design_type.lower())
        state_manager = StateManager()

        if not isinstance(state_manager.current_state, DesignState):
            state = DesignState()
            state.substate = design_enum
            state_manager.set_state(state)
        else:
            state_manager.current_state.substate = design_enum

        # Return the appropriate template based on design type
        if design_type == "combinatoric":
            return render_template("design/combinatoric.html")
        elif design_type == "genetic":
            return render_template("design/genetic.html")
        elif design_type == "mcmc":
            return render_template("design/mcmc.html")
        elif design_type == "fractal":
            return render_template("design/fractal.html")
        elif design_type == "generative":
            return render_template("design/generative.html")
        return render_template("coming_soon.html")
    except ValueError:
        return "Invalid design type", 404


@design_bp.route("/<design_type>/<method>")
def handle_design_method(design_type: str, method: str) -> Any:
    """Handle specific design method routes.

    Args:
        design_type: The type of design to handle
        method: The specific method to use

    Returns:
        Rendered method page
    """
    try:
        if design_type == "combinatoric":
            # For combinatoric methods, use method_base.html
            return render_template(
                "design/method_base.html",
                method_type=method,
                display_names=DISPLAY_NAMES,
                method_descriptions=METHOD_DESCRIPTIONS,
            )
        return render_template("coming_soon.html")
    except ValueError:
        return "Invalid method", 404


@design_bp.route("/upload-blocks", methods=["POST"])
def upload_blocks() -> Any:
    """Handle building blocks file upload.

    Returns:
        Redirect to design page with status message
    """
    if "blocks" not in request.files:
        flash("No file selected", "error")
        return redirect(url_for("design.index"))

    file: Any = request.files["blocks"]
    if file.filename == "":
        flash("No file selected", "error")
        return redirect(url_for("design.index"))

    def allowed_file(filename: str) -> bool:
        """Check if the file has an allowed extension.

        Args:
            filename: The name of the file to check.

        Returns:
            True if the file has an allowed extension, False otherwise.
        """
        ALLOWED_EXTENSIONS: Set[str] = {"csv", "xlsx", "txt"}
        return (
            "." in filename and filename.rsplit(".", 1)[1].lower() in ALLOWED_EXTENSIONS
        )

    if file and allowed_file(file.filename):
        filename: str = secure_filename(file.filename)
        # TODO: Process the building blocks file
        # This is where you would add the logic to parse and store the building blocks
        flash(f"Successfully loaded building blocks from {filename}", "success")
    else:
        flash("Invalid file type. Please upload a CSV, XLSX, or TXT file.", "error")

    return redirect(url_for("design.index"))


@design_bp.route("/explore-blocks")
def explore_blocks() -> Any:
    """Render the building blocks explorer page.

    Returns:
        Rendered building blocks explorer page
    """
    state_manager = StateManager()
    if not isinstance(state_manager.current_state, DesignState):
        state_manager.set_state(DesignState())

    # TODO: Get actual building blocks data
    building_blocks = []  # This will be replaced with real data

    return render_template(
        "design/explore_blocks.html", building_blocks=building_blocks
    )


@design_bp.route("/save-blocks")
def save_blocks() -> Any:
    """Handle building blocks export.

    Returns:
        File download response
    """
    # TODO: Implement actual building blocks export
    flash("Building blocks export not yet implemented", "error")
    return redirect(url_for("design.index"))


@design_bp.route("/generative")
def generative():
    return render_template("coming_soon.html", title="Generative Design")


@design_bp.route("/force_field")
def force_field():
    return render_template("coming_soon.html", title="Force Field Modeling")
