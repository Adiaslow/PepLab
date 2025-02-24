"""
Module for building blocks-related routes.

This module contains the routes for managing and exploring building blocks.
"""

from typing import Any, Dict, List
from flask import Blueprint, render_template

blocks_bp = Blueprint("blocks", __name__)

# Example database categories
BLOCK_DATABASES: Dict[str, str] = {
    "natural": "Natural Amino Acids",
    "unnatural": "Unnatural Amino Acids",
    "peptoids": "Peptoid Building Blocks",
    "linkers": "Linker Molecules",
    "specialty": "Specialty Building Blocks",
}


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
