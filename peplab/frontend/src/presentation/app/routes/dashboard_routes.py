# presentation/app/routes/dashboard_routes.py
"""Module for dashboard-related routes."""

# Standard Library Imports
from typing import Any

# External Imports
from flask import Blueprint, render_template

# Internal Imports
from peplab.frontend.src.infrastructure.states.dashboard_state import DashboardState

dashboard_bp: Blueprint = Blueprint("dashboard", __name__)


@dashboard_bp.route("/")
def dashboard() -> Any:
    """Render the dashboard page."""
    return render_template("dashboard/dashboard.html")


@dashboard_bp.route("/design")
def design() -> Any:
    """Render the design page."""
    return render_template("dashboard/design.html")


@dashboard_bp.route("/analysis")
def analysis() -> Any:
    """Render the analysis page."""
    return render_template("dashboard/analysis.html")


@dashboard_bp.route("/settings")
def settings() -> Any:
    """Render the settings page."""
    return render_template("dashboard/settings.html")
