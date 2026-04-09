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
        from peplab.backend.src.infrastructure.repositories.building_block_repository import BuildingBlockRepository
        repo = BuildingBlockRepository()
        
        try:
            db_blocks = repo.get_all_building_blocks()
            
            # Simple simulation of "sets" for MVP filter capability
            filtered = []
            canonical_names = ["Alanine", "Arginine", "Asparagine", "Aspartic Acid", "Cysteine", "Glutamic Acid", "Glutamine", "Glycine", "Histidine", "Isoleucine", "Leucine", "Lysine", "Methionine", "Phenylalanine", "Proline", "Serine", "Threonine", "Tryptophan", "Tyrosine", "Valine"]
            
            for b in db_blocks:
                props = b.properties or {}
                b_set = props.get("set")
                if not b_set:
                    # Fallback determination
                    if b.name in canonical_names:
                        b_set = "canonical"
                    else:
                        b_set = "user"
                
                # 'user' view shows everything they explicitly uploaded recently or non canonical
                if set_name == "user" and b_set != "canonical":
                    filtered.append(b)
                elif set_name == b_set:
                    filtered.append(b)
                # Fallback: if 'user' is selected, also show everything custom
                elif set_name == "user" and "set" not in props and b.name not in canonical_names:
                    filtered.append(b)
                    
            db_blocks = filtered
        except Exception:
            db_blocks = []

        building_blocks = []

        for row in db_blocks:
            smiles = row.metadata.get("smiles", "") if row.metadata else ""
            img_str = ""
            if smiles:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    # Generate 2D depiction
                    img = Draw.MolToImage(mol)
                    # Convert image to base64
                    img_buffer = io.BytesIO()
                    img.save(img_buffer, format="PNG")
                    img_str = base64.b64encode(img_buffer.getvalue()).decode()

            props = row.properties or {}
            # Create building block object
            building_block = {
                "name": row.name,
                "code": props.get("alt_name1", ""),
                "alt_code": props.get("alt_name2", ""),
                "smiles": smiles,
                "position": props.get("position", ""),
                "image": f"data:image/png;base64,{img_str}" if img_str else "",
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
    try:
        from peplab.backend.src.infrastructure.repositories.building_block_repository import BuildingBlockRepository
        repo = BuildingBlockRepository()
        blocks = repo.get_all_building_blocks()
    except Exception:
        blocks = []
    return render_template("blocks/library.html", blocks=blocks)


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
