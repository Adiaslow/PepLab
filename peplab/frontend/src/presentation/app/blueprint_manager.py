# peplab/frontend/src/presentation/app/blueprint_manager.py
"""
This module manages the registration and coordination of Flask blueprints
with the application state management system.
"""

# External imports
from flask import Flask
from typing import Dict, Any
from dataclasses import dataclass

# Internal imports
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
from peplab.frontend.src.infrastructure.orchestrator import (
    ApplicationOrchestrator,
    ApplicationContext,
)


@dataclass
class BlueprintConfig:
    blueprint: Any
    url_prefix: str


class BlueprintManager:
    """Manages Flask blueprints and their integration with the state system.

    Attributes:
        app: The Flask application.
        orchestrator: The application orchestrator.
        blueprints: A dictionary of blueprints and their URL prefixes.

    Methods:
        __init__: Initialize the BlueprintManager.
        register_blueprints: Register all blueprints with the Flask application.
        initialize_state_handlers: Initialize state handlers for all blueprints.
        _setup_state_handlers: Set up state handlers for a specific blueprint.
        before_request: Update application state before handling request.
        after_request: Handle any necessary state cleanup after request.
    """

    def __init__(self, app: Flask) -> None:
        """Initialize the BlueprintManager.

        Args:
            app: The Flask application.
        """
        self.app: Flask = app
        self.orchestrator: ApplicationOrchestrator = ApplicationOrchestrator()
        self.blueprints: Dict[str, BlueprintConfig] = {
            "main": BlueprintConfig(main_bp, "/"),
            "dashboard": BlueprintConfig(dashboard_bp, "/dashboard"),
            "design": BlueprintConfig(design_bp, "/design"),
            "analysis": BlueprintConfig(analysis_bp, "/analysis"),
            "modeling": BlueprintConfig(modeling_bp, "/modeling"),
            "optimization": BlueprintConfig(optimization_bp, "/optimization"),
            "settings": BlueprintConfig(settings_bp, "/settings"),
            "blocks": BlueprintConfig(blocks_bp, "/blocks"),
        }

    def register_blueprints(self) -> None:
        """Register all blueprints with the Flask application.

        Args:
            self: The BlueprintManager instance.
        """
        for config in self.blueprints.values():
            self.app.register_blueprint(config.blueprint, url_prefix=config.url_prefix)

    def initialize_state_handlers(self) -> None:
        """Initialize state handlers for all blueprints.

        Args:
            self: The BlueprintManager instance.
        """
        for blueprint_name, config in self.blueprints.items():
            self._setup_state_handlers(config.blueprint, blueprint_name)

    def _setup_state_handlers(self, blueprint: Any, blueprint_name: str) -> None:
        """Set up before/after request handlers for state management.

        Args:
            blueprint: The blueprint to set up state handlers for.
            blueprint_name: The name of the blueprint.
        """

        @blueprint.before_request
        def before_request() -> None:
            """Update application state before handling request.

            Args:
                self: The BlueprintManager instance.
            """
            context: ApplicationContext = self.orchestrator.get_application_context()
            # Store context in g for access in route handlers
            from flask import g

            g.app_context = context

        @blueprint.after_request
        def after_request(response: Any) -> Any:
            """Handle any necessary state cleanup after request.

            Args:
                response: The response to the request.
            """
            return response


def register_blueprints(app: Flask) -> None:
    """
    Register all blueprints with the Flask application.

    Args:
        app: Flask application instance
    """
    manager = BlueprintManager(app)
    manager.register_blueprints()
