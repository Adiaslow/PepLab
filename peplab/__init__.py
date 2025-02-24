"""
Initialize the PepLab Flask application.
"""

import os
from flask import Flask
from flask_cors import CORS
from flask_sqlalchemy import SQLAlchemy
from peplab.config import Config

db = SQLAlchemy()


def create_app(config_class=Config):
    """
    Create and configure the Flask application instance.

    Args:
        config_class: Configuration class to use (default: Config)

    Returns:
        Flask application instance
    """
    app = Flask(
        __name__,
        template_folder="frontend/src/presentation/app/templates",
        static_folder="frontend/src/presentation/app/static",
    )

    # Enable CORS for all routes
    CORS(app)

    # Basic configuration
    app.config["ENV"] = os.getenv("FLASK_ENV", "development")
    app.config["DEBUG"] = os.getenv("FLASK_DEBUG", "1") == "1"
    app.config["SECRET_KEY"] = os.getenv("SECRET_KEY", "527-436-828")

    # Security settings
    app.config["ALLOWED_HOSTS"] = os.getenv(
        "ALLOWED_HOSTS", "localhost,127.0.0.1"
    ).split(",")

    app.config.from_object(config_class)

    db.init_app(app)

    # Register all blueprints using the blueprint manager
    from peplab.frontend.src.presentation.app.blueprint_manager import (
        register_blueprints,
    )

    register_blueprints(app)

    return app


app = create_app()
