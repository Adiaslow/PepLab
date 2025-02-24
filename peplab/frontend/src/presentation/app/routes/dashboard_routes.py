# presentation/app/routes/dashboard_routes.py
"""
Module for dashboard-related routes.

This module contains the routes for the dashboard page.

Classes:
    DashboardBlueprint: Blueprint for dashboard routes

Functions:
    dashboard: Render the dashboard page
    design: Render the design dashboard page
    analysis: Render the analysis dashboard page
    modeling: Render the modeling dashboard page
    optimization: Render the optimization dashboard page
"""
# Standard Library Imports
from typing import Any, Set
import os

# External Imports
from werkzeug.utils import secure_filename
from werkzeug.datastructures import FileStorage
from flask import (
    Blueprint,
    render_template,
    request,
    flash,
    redirect,
    url_for,
    send_file,
)
import io
import csv

# Internal Imports
from peplab.frontend.src.infrastructure.managers.state_manager import StateManager
from peplab.frontend.src.infrastructure.states.dashboard_state import DashboardState

dashboard_bp: Blueprint = Blueprint("dashboard", __name__)

ALLOWED_EXTENSIONS: Set[str] = {"csv", "xlsx", "txt"}


def allowed_file(filename: str) -> bool:
    """Check if the file extension is allowed.

    Args:
        filename: Name of the file to check

    Returns:
        bool: True if file extension is allowed
    """
    return "." in filename and filename.rsplit(".", 1)[1].lower() in ALLOWED_EXTENSIONS


@dashboard_bp.route("/")
def dashboard() -> Any:
    """Render the dashboard hub page.

    Returns:
        Rendered dashboard hub page
    """
    state_manager = StateManager()
    if not isinstance(state_manager.current_state, DashboardState):
        state_manager.set_state(DashboardState())

    # TODO: Get actual statistics from your data store
    stats = {
        "peptides_count": 0,  # Total number of peptides in library
        "building_blocks_count": 0,  # Number of unique building blocks
        "properties_count": 0,  # Number of calculated/measured properties
        "recent_projects": [],  # List of recent projects
    }

    return render_template("dashboard/dashboard.html", stats=stats)


@dashboard_bp.route("/design")
def design() -> Any:
    """Render the design dashboard page.

    Returns:
        Rendered design dashboard page
    """
    return render_template("dashboard/design.html")


@dashboard_bp.route("/analysis")
def analysis() -> Any:
    """Render the analysis dashboard page.

    Returns:
        Rendered analysis dashboard page
    """
    return render_template("dashboard/analysis.html")


@dashboard_bp.route("/modeling")
def modeling() -> Any:
    """Render the modeling dashboard page.

    Returns:
        Rendered modeling dashboard page
    """
    return render_template("dashboard/modeling.html")


@dashboard_bp.route("/optimization")
def optimization() -> Any:
    """Render the optimization dashboard page.

    Returns:
        Rendered optimization dashboard page
    """
    return render_template("dashboard/optimization.html")


@dashboard_bp.route("/upload-library", methods=["POST"])
def upload_library() -> Any:
    """Handle library file upload.

    Returns:
        Redirect to dashboard with status message
    """
    if "library" not in request.files:
        flash("No file selected", "error")
        return redirect(url_for("dashboard.index"))

    file: FileStorage = request.files["library"]
    if not file.filename:  # Type check for None
        flash("No file selected", "error")
        return redirect(url_for("dashboard.index"))

    if allowed_file(file.filename):  # Now we know filename is not None
        filename: str = secure_filename(file.filename)
        # TODO: Process the library file
        # This is where you would add the logic to parse and process the library file
        flash(f"Successfully loaded library from {filename}", "success")
    else:
        flash("Invalid file type. Please upload a CSV, XLSX, or TXT file.", "error")

    return redirect(url_for("dashboard.index"))


@dashboard_bp.route("/save-library")
def save_library() -> Any:
    """Handle library file download.

    Returns:
        File download response
    """
    try:
        # TODO: Get the current library data from your state/database
        # This is a placeholder that creates a sample CSV
        output = io.StringIO()
        writer = csv.writer(output)
        writer.writerow(["Sequence", "MW", "LogP"])  # Headers
        writer.writerow(["ACDEFGH", "823.32", "1.23"])  # Sample data

        # Create the response
        output.seek(0)
        return send_file(
            io.BytesIO(output.getvalue().encode("utf-8")),
            mimetype="text/csv",
            as_attachment=True,
            download_name="peptide_library.csv",
        )
    except Exception as e:
        flash(f"Error saving library: {str(e)}", "error")
        return redirect(url_for("dashboard.index"))
