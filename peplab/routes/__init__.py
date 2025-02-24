"""
Main routes for the PepLab application.
"""

from flask import Blueprint

main = Blueprint("main", __name__)


@main.route("/")
def index():
    return "Welcome to PepLab"
