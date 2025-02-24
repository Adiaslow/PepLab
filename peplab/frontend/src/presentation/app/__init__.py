# frontend/src/presentation/app/__init__.py
"""
Flask application factory.
"""

from flask import Flask
from flask_cors import CORS
from peplab.frontend.src.presentation.app.routes import (
    main_bp,
    dashboard_bp,
    design_bp,
    analysis_bp,
    modeling_bp,
    optimization_bp,
    settings_bp,
    blocks_bp,
)
from peplab.frontend.src.presentation.app.utils.navigation import get_parent_route
from typing import Dict, Any
from flask_assets import Environment, Bundle
import os
from flask_wtf.csrf import CSRFProtect


def create_app(config=None) -> Flask:
    """Create and configure the Flask application.

    Args:
        config: Configuration object

    Returns:
        Flask application instance
    """
    app = Flask(__name__)
    CORS(app, supports_credentials=True)

    # Set secret key
    app.config["SECRET_KEY"] = os.environ.get(
        "SECRET_KEY", "dev-secret-key"
    )  # Use environment variable in production
    app.config["WTF_CSRF_CHECK_DEFAULT"] = False  # We'll check manually in our routes
    app.config["WTF_CSRF_TIME_LIMIT"] = None  # Optional: disable CSRF token timeout

    # Initialize CSRF protection
    csrf = CSRFProtect(app)

    # Initialize Flask-Assets
    assets = Environment(app)
    assets.load_path = [
        os.path.join(os.path.dirname(__file__), "static"),
        os.path.join(os.path.dirname(__file__), "static/css"),
    ]

    # Create CSS bundle
    css = Bundle(
        "base/_variables.css",
        "base/_reset.css",
        "base/_typography.css",
        "base/_layout.css",
        "components/_buttons.css",
        "components/_cards.css",
        "components/_forms.css",
        "components/_navigation.css",
        "components/_flash.css",
        "pages/_dashboard.css",
        "pages/_design.css",
        "pages/_splash.css",
        "pages/_coming_soon.css",
        "themes/_neuomorphic.css",
        "main.css",
        filters="cssmin",
        output="gen/packed.css",
    )
    assets.register("css_all", css)

    # Configure the app
    if config:
        app.config.from_object(config)

    # Register blueprints
    app.register_blueprint(main_bp, url_prefix="/")
    app.register_blueprint(dashboard_bp, url_prefix="/dashboard")
    app.register_blueprint(design_bp, url_prefix="/design")
    app.register_blueprint(analysis_bp, url_prefix="/analysis")
    app.register_blueprint(modeling_bp, url_prefix="/modeling")
    app.register_blueprint(optimization_bp, url_prefix="/optimization")
    app.register_blueprint(settings_bp, url_prefix="/settings")
    app.register_blueprint(blocks_bp, url_prefix="/blocks")

    # Add template context processors
    @app.context_processor
    def utility_processor() -> Dict[str, Any]:
        return {"get_parent_route": get_parent_route}

    return app
