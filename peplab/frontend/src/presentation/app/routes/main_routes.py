# routes/main_routes.py
"""Main routes for the application.

Classes:
    MainBlueprint: Blueprint for main routes

Functions:
    index: Render the main index (splash) page
    dashboard: Render the dashboard page
    coming_soon: Render the coming soon page
"""

# Standard Library Imports
from typing import Any, Dict

# External Imports
from flask import Blueprint, render_template, redirect, url_for

main_bp = Blueprint("main", __name__)


@main_bp.route("/")
def index() -> Any:
    """Render the main index (splash) page.

    Returns:
        Rendered splash page
    """
    return render_template("index.html")


@main_bp.route("/dashboard")
def dashboard() -> Any:
    """Render the dashboard page.

    Returns:
        Rendered dashboard page
    """
    # TODO: Get actual statistics from your data store
    stats: Dict[str, int] = {
        "peptides_count": 0,  # Total number of peptides in library
        "building_blocks_count": 0,  # Number of unique building blocks
        "properties_count": 0,  # Number of calculated/measured properties
        "recent_projects_count": 0,  # Number of recent projects
    }

    # This is the hub that branches to all other workflows
    return render_template("dashboard/dashboard.html", stats=stats)


@main_bp.route("/coming-soon")
def coming_soon() -> Any:
    """Render the coming soon page.

    Returns:
        Rendered coming soon page
    """
    return render_template("coming_soon.html")
