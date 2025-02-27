"""
Module for building blocks-related routes.

This module contains the routes for managing and exploring building blocks.
"""

from typing import Any, Dict, List
from flask import Blueprint, render_template, jsonify, redirect, url_for
from rdkit import Chem
from rdkit.Chem import Draw
import pandas as pd
import base64
import io

blocks_bp = Blueprint("blocks", __name__)

# Example database categories
BLOCK_DATABASES: Dict[str, str] = {
    "natural": "Natural Amino Acids",
    "unnatural": "Unnatural Amino Acids",
    "peptoids": "Peptoid Building Blocks",
    "linkers": "Linker Molecules",
    "specialty": "Specialty Building Blocks",
}


@blocks_bp.route("/api/building-blocks/<set_name>")
def get_building_blocks(set_name: str) -> Any:
    """Get building blocks for a specific set.

    Args:
        set_name: Name of the building block set to load

    Returns:
        JSON response with building blocks data
    """
    try:
        # Map set names to file paths (you can expand this)
        set_files = {
            "canonical": "test_building_blocks.csv",
            # Add more mappings as needed
        }

        if set_name not in set_files:
            return jsonify({"error": f"Unknown set: {set_name}"}), 404

        # Read the CSV file
        df = pd.read_csv(set_files[set_name])
        building_blocks = []

        for _, row in df.iterrows():
            # Create RDKit molecule from SMILES
            mol = Chem.MolFromSmiles(row["smiles"])

            # Generate 2D depiction
            img = Draw.MolToImage(mol)

            # Convert image to base64
            img_buffer = io.BytesIO()
            img.save(img_buffer, format="PNG")
            img_str = base64.b64encode(img_buffer.getvalue()).decode()

            # Create building block object
            building_block = {
                "name": row["name"],
                "code": row["alt_name1"],
                "alt_code": row["alt_name2"],
                "smiles": row["smiles"],
                "position": row["position"],
                "image": f"data:image/png;base64,{img_str}",
                "set": set_name,
            }
            building_blocks.append(building_block)

        return jsonify({"building_blocks": building_blocks})
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@blocks_bp.route("/explore")
def explore() -> Any:
    """Render the building blocks explorer page.

    Returns:
        Rendered building blocks explorer page
    """
    return render_template("blocks/explore.html", databases=BLOCK_DATABASES)


@blocks_bp.route("/database/<db_type>")
def view_database(db_type: str) -> Any:
    """View specific database of building blocks.

    Args:
        db_type: Type of database to view

    Returns:
        Rendered database view page
    """
    if db_type not in BLOCK_DATABASES:
        return "Invalid database type", 404

    return render_template(
        "blocks/database.html", db_type=db_type, db_name=BLOCK_DATABASES[db_type]
    )


@blocks_bp.route("/library")
def view_library() -> Any:
    """View user's loaded building blocks.

    Returns:
        Rendered library view page
    """
    # TODO: Get actual library data
    return render_template("blocks/library.html")


@blocks_bp.route("/upload-blocks", methods=["POST"])
def upload_blocks() -> Any:
    """Handle building blocks upload and redirect to explorer.

    Returns:
        Redirect to explorer page
    """
    try:
        # Handle file upload logic here
        # ... existing upload logic ...

        # Redirect to explore page with query parameter
        return redirect(url_for("blocks.explore", show_user_library=True))
    except Exception as e:
        return jsonify({"error": str(e)}), 500
