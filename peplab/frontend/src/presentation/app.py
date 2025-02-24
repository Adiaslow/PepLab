"""
This file is the main entry point for the application.
It is responsible for setting up the Flask application and routing the requests.
"""

# Standard Library Imports
from typing import Any
import os
import multiprocessing
import torch  # for GPU detection

# External Imports
from flask import Flask
from flask_assets import Environment, Bundle

# Internal Imports
from peplab.frontend.src.infrastructure.managers.state_manager import StateManager
from peplab.frontend.src.presentation.app.routes.main_routes import main_bp
from peplab.frontend.src.presentation.app.routes.design_routes import design_bp
from peplab.frontend.src.presentation.app.routes.analysis_routes import analysis_bp
from peplab.frontend.src.presentation.app.routes.modeling_routes import modeling_bp
from peplab.frontend.src.presentation.app.routes.optimization_routes import (
    optimization_bp,
)
from peplab.frontend.src.presentation.app.routes.settings_routes import settings_bp
from peplab.frontend.src.presentation.app.config import Config
from peplab.frontend.src.presentation.app.extensions import init_extensions


def create_app(config_class=Config) -> Flask:
    """Create and configure the Flask application."""
    app = Flask(__name__, template_folder="app/templates", static_folder="app/static")

    # Load configuration
    app.config.from_object(config_class)

    # Initialize extensions
    init_extensions(app)

    # Register blueprints
    app.register_blueprint(main_bp)
    app.register_blueprint(design_bp, url_prefix="/design")
    app.register_blueprint(analysis_bp, url_prefix="/analysis")
    app.register_blueprint(modeling_bp, url_prefix="/modeling")
    app.register_blueprint(optimization_bp, url_prefix="/optimization")
    app.register_blueprint(settings_bp, url_prefix="/settings")

    return app


if __name__ == "__main__":
    app = create_app()
    app.run(debug=True)
