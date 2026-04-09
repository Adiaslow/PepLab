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

    from peplab import db
    from peplab.backend.src.infrastructure.database.models import BuildingBlockModel, PeptideModel, LibraryModel

    stats = {
        "peptides_count": db.session.query(PeptideModel).count(),
        "blocks_count": db.session.query(BuildingBlockModel).count(),
        "properties_count": db.session.query(LibraryModel).count(),  # Switched to Library count for utility
        "recent_projects": [],  
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
        return redirect(url_for("dashboard.dashboard"))

    file: FileStorage = request.files["library"]
    if not file.filename:  # Type check for None
        flash("No file selected", "error")
        return redirect(url_for("dashboard.dashboard"))

    if allowed_file(file.filename):  # Now we know filename is not None
        filename: str = secure_filename(file.filename)
        try:
            import csv
            import io
            from peplab import db
            from peplab.backend.src.infrastructure.database.models import LibraryModel, PeptideModel
            import uuid

            stream = io.StringIO(file.read().decode("utf-8", errors="ignore"), newline=None)
            reader = csv.DictReader(stream, skipinitialspace=True)
            
            # Use 'Sequence' if exists, otherwise fallback to first column
            fieldnames = reader.fieldnames if reader.fieldnames else []
            seq_col = "Sequence" if "Sequence" in fieldnames else (fieldnames[0] if fieldnames else None)
            
            if not seq_col:
                raise ValueError("Could not determine sequence column from CSV.")
                
            new_lib = LibraryModel(
                id=str(uuid.uuid4()),
                name=f"Imported Library: {filename}",
                description="Uploaded from dashboard"
            )
            db.session.add(new_lib)
            
            count = 0
            for row in reader:
                seq_val = row.get(seq_col, "").strip()
                if not seq_val:
                    continue
                pep = PeptideModel(
                    id=str(uuid.uuid4()),
                    name=seq_val,
                    library_id=new_lib.id
                )
                db.session.add(pep)
                count += 1
                
            db.session.commit()
            flash(f"Successfully loaded {count} peptides into new database library from {filename}", "success")
        except Exception as e:
            flash(f"Error processing library CSV: {str(e)}", "error")
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
        from peplab import db
        from peplab.backend.src.infrastructure.database.models import LibraryModel
        from peplab.backend.src.application.services.analysis.cheminformatics_service import CheminformaticsService
        import csv
        import io

        latest_lib = db.session.query(LibraryModel).order_by(LibraryModel.id.desc()).first()
        if not latest_lib:
            flash("No libraries found to export", "error")
            return redirect(url_for("dashboard.dashboard"))

        # Reconstruct sequences: in MVP we map unique building blocks from the peptide model, or just use peptide name
        sequences = []
        for peptide in latest_lib.peptides:
            # name is formatted like "Alanine Glycine"
            seq = peptide.name.split(" ") if peptide.name else []
            if seq:
                sequences.append(seq)
                
        # Generate properties
        service = CheminformaticsService()
        results = service.analyze_sequences(sequences)

        output = io.StringIO()
        writer = csv.writer(output)
        writer.writerow(["Sequence", "Mol. Wt (Da)", "Exact Mass", "LogP", "TPSA", "H-Donors", "H-Acceptors", "Status"])
        
        for res in results:
            seq_str = " ".join(res["sequence"])
            if res["status"] == "Success":
                writer.writerow([
                    seq_str,
                    res["molecular_weight"],
                    res["exact_mass"],
                    res["log_p"],
                    res["tpsa"],
                    res["h_donors"],
                    res["h_acceptors"],
                    res["status"]
                ])
            else:
                 writer.writerow([seq_str, "-", "-", "-", "-", "-", "-", res["status"]])

        output.seek(0)
        return send_file(
            io.BytesIO(output.getvalue().encode("utf-8")),
            mimetype="text/csv",
            as_attachment=True,
            download_name=f"{latest_lib.name.replace(' ', '_')}_library.csv",
        )
    except Exception as e:
        flash(f"Error saving library: {str(e)}", "error")
        return redirect(url_for("dashboard.dashboard"))

