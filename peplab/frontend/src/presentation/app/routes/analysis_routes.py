# presentation/app/routes/analysis.py
"""
Module for analysis-related routes.

This module contains the routes for the analysis page.
"""

# Standard Library Imports
from typing import Any

# External Imports
from flask import Blueprint, render_template

# Internal Imports
from infrastructure.managers.state_manager import StateManager

analysis_bp: Blueprint = Blueprint("analysis", __name__)


@analysis_bp.route("/")
def analysis() -> Any:
    """Render the analysis page."""
    return render_template("analysis/analysis.html")


@analysis_bp.route("/cheminformatics")
def cheminformatics() -> Any:
    """Render the cheminformatics page."""
    return render_template("analysis/cheminformatics.html")


@analysis_bp.route("/data_analysis")
def data_analysis() -> Any:
    """Render the data analysis page."""
    return render_template("analysis/data_analysis.html")


@analysis_bp.route("/machine_learning")
def machine_learning() -> Any:
    """Render the machine learning page."""
    return render_template("analysis/machine_learning.html")
