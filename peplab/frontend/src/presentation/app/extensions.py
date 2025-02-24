# presentation/app/extensions.py
"""Flask extensions initialization.

This module initializes the Flask extensions for the application.

Methods:
    init_extensions(app: Flask) -> None: Initialize the Flask extensions.
"""

# External Imports
from flask import Flask

# Flask Extension Imports
from flask_postgresql import PostgreSQL
from flask_migrate import Migrate
from flask_wtf.csrf import CSRFProtect

# Initialize Flask extensions
db: PostgreSQL = PostgreSQL()
migrate: Migrate = Migrate()
csrf: CSRFProtect = CSRFProtect()


def init_extensions(app: Flask) -> None:
    """Initialize Flask extensions.

    Args:
        app: The Flask application.
    """
    db.init_app(app)
    migrate.init_app(app)
    csrf.init_app(app)
