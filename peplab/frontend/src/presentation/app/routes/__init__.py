# frontend/src/presentation/app/routes/__init__.py
"""
This module contains the routes for the application.
"""

# Standard Library Imports
from typing import List

# Internal Imports
from peplab.frontend.src.presentation.app.routes.main_routes import main_bp
from peplab.frontend.src.presentation.app.routes.dashboard_routes import dashboard_bp
from peplab.frontend.src.presentation.app.routes.design_routes import design_bp
from peplab.frontend.src.presentation.app.routes.analysis_routes import analysis_bp
from peplab.frontend.src.presentation.app.routes.modeling_routes import modeling_bp
from peplab.frontend.src.presentation.app.routes.optimization_routes import (
    optimization_bp,
)
from peplab.frontend.src.presentation.app.routes.settings_routes import settings_bp
from peplab.frontend.src.presentation.app.routes.blocks_routes import blocks_bp

__all__: List[str] = [
    "main_bp",
    "dashboard_bp",
    "design_bp",
    "analysis_bp",
    "modeling_bp",
    "optimization_bp",
    "settings_bp",
    "blocks_bp",
]
