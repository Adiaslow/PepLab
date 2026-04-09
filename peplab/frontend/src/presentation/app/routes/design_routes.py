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
from flask import Blueprint, render_template, request, flash, redirect, url_for, jsonify
from werkzeug.utils import secure_filename
import requests
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
}

# Display names for methods
DISPLAY_NAMES: Dict[str, str] = {
    "combinatoric": "Combinatorial",
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
        return redirect(url_for("design.design"))

    file: Any = request.files["blocks"]
    if file.filename == "":
        flash("No file selected", "error")
        return redirect(url_for("design.design"))

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
        try:
            import csv
            import io
            from peplab.backend.src.domain.models.building_block import BuildingBlock
            from peplab.backend.src.infrastructure.repositories.building_block_repository import BuildingBlockRepository
            
            repo = BuildingBlockRepository()
            stream = io.StringIO(file.read().decode("utf-8", errors="ignore"), newline=None)
            reader = csv.DictReader(stream, skipinitialspace=True)
            
            count = 0
            for row in reader:
                name = row.get("name", "").strip()
                if not name:
                    continue
                
                bb = BuildingBlock(
                    name=name,
                    properties={
                        "alt_name1": row.get("alt_name1", "").strip(),
                        "alt_name2": row.get("alt_name2", "").strip(),
                        "position": row.get("position", "").strip()
                    },
                    metadata={
                        "smiles": row.get("smiles", "").strip()
                    }
                )
                
                try:
                    repo.get_building_block_by_name(name)
                except Exception:
                    # Not found, add it
                    repo.add_building_block(bb)
                    count += 1
                    
            flash(f"Successfully loaded {count} new building blocks from {filename}", "success")
        except Exception as e:
            flash(f"Error processing CSV: {str(e)}", "error")
    else:
        flash("Invalid file type. Please upload a CSV, XLSX, or TXT file.", "error")

    return redirect(url_for("design.design"))


@design_bp.route("/explore-blocks")
def explore_blocks() -> Any:
    """Render the building blocks explorer page.

    Returns:
        Rendered building blocks explorer page
    """
    state_manager = StateManager()
    if not isinstance(state_manager.current_state, DesignState):
        state_manager.set_state(DesignState())

    # Redirect to the dedicated blocks explorer (blocks_routes.py) — it already
    # has the full toggle/search/preview UI and pulls from the DB correctly.
    return redirect(url_for("blocks.explore"))


@design_bp.route("/save-blocks")
def save_blocks() -> Any:
    """Handle building blocks export.

    Returns:
        File download response
    """
    try:
        from peplab.backend.src.infrastructure.repositories.building_block_repository import BuildingBlockRepository
        from flask import Response
        import csv
        import io
        
        repo = BuildingBlockRepository()
        
        # Adding simple fallback empty array if get_all_building_blocks errors
        try:
            blocks = repo.get_all_building_blocks()
        except:
            blocks = []

        si = io.StringIO()
        writer = csv.writer(si)
        writer.writerow(["name", "alt_name1", "alt_name2", "position", "smiles"])
        
        for b in blocks:
            props = b.properties or {}
            metadata = b.metadata or {}
            writer.writerow([
                b.name,
                props.get("alt_name1", ""),
                props.get("alt_name2", ""),
                props.get("position", ""),
                metadata.get("smiles", "")
            ])
            
        return Response(
            si.getvalue(),
            mimetype="text/csv",
            headers={"Content-Disposition": "attachment; filename=building_blocks.csv"}
        )
    except Exception as e:
        flash(f"Failed to export building blocks: {str(e)}", "error")
        return redirect(url_for("design.design"))






    
@design_bp.route("/design/combinatoric/combination")
def combination_design():
    """Render the Combination design page."""
    return render_template("combination.html")