# presentation/app/routes/analysis.py
"""
Module for analysis-related routes.

This module contains the routes for the analysis page.
"""

# Standard Library Imports
from typing import Any, Dict

# External Imports
from flask import Blueprint, render_template

# Internal Imports
from peplab.frontend.src.infrastructure.managers.state_manager import StateManager
from peplab.frontend.src.infrastructure.states.analysis_state import AnalysisState
from peplab.frontend.src.core.types.analysis_types import AnalysisType

analysis_bp = Blueprint("analysis", __name__)

# Analysis method descriptions
METHOD_DESCRIPTIONS: Dict[str, str] = {
    "cheminformatic": "Analyze chemical properties of peptides",
    "data_analysis": "Analyze experimental results from assays, DEL sequencing, etc.",
    "machine_learning": "Analyze peptides using machine learning models",
    "simulation": "Analyze simulation results",
}

# Display names for methods
DISPLAY_NAMES: Dict[str, str] = {
    "cheminformatic": "Cheminformatics",
    "data_analysis": "Data Analysis",
    "machine_learning": "Machine Learning",
    "simulation": "Simulation Analysis",
}


@analysis_bp.route("/")
def analysis() -> Any:
    """Render the analysis hub page.

    Returns:
        Rendered analysis hub page
    """
    state_manager = StateManager()
    if not isinstance(state_manager.current_state, AnalysisState):
        state_manager.set_state(AnalysisState())

    return render_template(
        "analysis/analysis.html",
        method_descriptions=METHOD_DESCRIPTIONS,
        display_names=DISPLAY_NAMES,
    )


@analysis_bp.route("/<analysis_type>")
def handle_analysis(analysis_type: str) -> Any:
    """Handle specific analysis type routes.

    Args:
        analysis_type: The type of analysis to handle

    Returns:
        Rendered analysis type page or coming soon page
    """
    try:
        analysis_enum = AnalysisType(analysis_type.lower())
        state_manager = StateManager()

        if not isinstance(state_manager.current_state, AnalysisState):
            state = AnalysisState()
            state.substate = analysis_enum
            state_manager.set_state(state)
        else:
            state_manager.current_state.substate = analysis_enum

        return render_template("coming_soon.html")
    except ValueError:
        return "Invalid analysis type", 404
