# presentation/app/routes/settings.py
"""
Module for settings-related routes.

This module contains the routes for the settings page.

Classes:
    SettingsBlueprint: Blueprint for settings routes

Functions:
    settings: Render the settings page
    save_general_settings: Save general settings
    save_compute_settings: Save compute settings
    save_appearance_settings: Save appearance settings
    clear_data: Clear application data
"""

# Standard Library Imports
from typing import Any

# External Imports
from flask import Blueprint, render_template, request, flash, redirect, url_for
import multiprocessing
import torch

# Internal Imports
from peplab.frontend.src.infrastructure.managers.state_manager import StateManager

settings_bp = Blueprint("settings", __name__)


@settings_bp.route("/")
def settings() -> Any:
    """Render the settings page.

    Returns:
        Rendered settings page
    """
    # Get number of CPUs
    num_cpus = multiprocessing.cpu_count()
    # Check for GPU availability
    gpu_available = (
        torch.cuda.is_available() if hasattr(torch.cuda, "is_available") else False
    )

    return render_template(
        "settings/settings.html", num_cpus=num_cpus, gpu_available=gpu_available
    )


@settings_bp.route("/save-general", methods=["POST"])
def save_general() -> Any:
    """Save general settings."""
    try:
        # TODO: Save settings to configuration
        default_format = request.form.get("default_format")
        auto_save = request.form.get("auto_save")
        flash("General settings saved successfully", "success")
    except Exception as e:
        flash(f"Error saving settings: {str(e)}", "error")
    return redirect(url_for("settings.index"))


@settings_bp.route("/save-compute", methods=["POST"])
def save_compute() -> Any:
    """Save computation settings."""
    try:
        # TODO: Save settings to configuration
        max_threads = request.form.get("max_threads")
        gpu_enabled = request.form.get("gpu_enabled") == "on"
        flash("Computation settings saved successfully", "success")
    except Exception as e:
        flash(f"Error saving settings: {str(e)}", "error")
    return redirect(url_for("settings.index"))


@settings_bp.route("/save-appearance", methods=["POST"])
def save_appearance() -> Any:
    """Save appearance settings."""
    try:
        # TODO: Save settings to configuration
        theme = request.form.get("theme")
        compact_view = request.form.get("compact_view") == "on"
        flash("Appearance settings saved successfully", "success")
    except Exception as e:
        flash(f"Error saving settings: {str(e)}", "error")
    return redirect(url_for("settings.index"))


@settings_bp.route("/clear-data", methods=["POST"])
def clear_data() -> Any:
    """Clear application data based on action.

    Returns:
        Redirect to settings page with status message
    """
    action = request.form.get("action")
    try:
        if action == "clear_cache":
            # TODO: Implement cache clearing
            flash("Cache cleared successfully", "success")
        elif action == "reset_settings":
            # TODO: Implement settings reset
            flash("Settings reset to defaults", "success")
    except Exception as e:
        flash(f"Error clearing data: {str(e)}", "error")
    return redirect(url_for("settings.index"))
