# presentation/app/routes/settings.py
"""
Module for settings-related routes.

This module contains the routes for the settings page.
"""

# Standard Library Imports
from typing import Any

# External Imports
from flask import Blueprint, render_template

# Internal Imports
from infrastructure.managers.state_manager import StateManager

settings_bp = Blueprint("settings", __name__)


@settings_bp.route("/")
def settings() -> Any:
    """Render the settings page."""
    return render_template("settings.html")
