# routes/main_routes.py
"""Module for main application routes."""
from typing import Any

# External Imports
from flask import Blueprint, render_template

main_bp: Blueprint = Blueprint("main", __name__)


@main_bp.route("/")
def index() -> Any:
    """Render the index page."""
    return render_template("index.html")


@main_bp.route("/dashboard")
def dashboard() -> Any:
    """Render the dashboard page."""
    return render_template("dashboard.html")
