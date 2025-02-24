"""Flask application entry point."""

# External Imports
from flask import Flask, Blueprint, render_template

# Internal Imports
from peplab.frontend.src.presentation.app import create_app

app: Flask = create_app()

# Main routes for the PepLab application.
main = Blueprint("main", __name__)


@main.route("/")
def index():
    return render_template("index.html")


@main.route("/dashboard")
def dashboard():
    return render_template("dashboard.html")


@main.route("/design")
def design():
    return render_template("design.html")


@main.route("/analysis")
def analysis():
    return render_template("analysis.html")


@main.route("/modeling")
def modeling():
    return render_template("modeling.html")


@main.route("/optimization")
def optimization():
    return render_template("optimization.html")


@main.route("/settings")
def settings():
    return render_template("settings.html")


# ... your other routes ...

if __name__ == "__main__":
    app.run(debug=True, port=5001)
