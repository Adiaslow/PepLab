# peplab/run.py
"""
Development server script.

This script sets up the development environment and runs the Flask application.
It checks for required environment variables and sets sensible defaults if needed.
"""

# Standard Library Imports
import os
import sys
from pathlib import Path
from typing import Dict, List

# External Imports
from flask import Flask
from dotenv import load_dotenv

# Internal Imports
from peplab import create_app


def get_project_root() -> Path:
    """Get the project root directory.

    Returns:
        Path: The project root directory
    """
    return Path(__file__).parent.parent


def check_environment() -> Dict[str, str]:
    """
    Check and set up environment variables.

    Returns:
        Dict[str, str]: Dictionary of environment variables

    Raises:
        SystemExit: If critical environment variables are missing
    """
    # Load environment variables from .flaskenv file
    env_path: Path = get_project_root() / ".flaskenv"
    load_dotenv(env_path)

    required_vars: Dict[str, str] = {
        "FLASK_APP": "wsgi.py",
        "FLASK_ENV": "development",
        "FLASK_DEBUG": "1",
        "APP_PORT": "5001",
        "APP_HOST": "localhost",
    }

    env_vars: Dict[str, str] = {}
    missing_vars: List[str] = []

    for var, default in required_vars.items():
        value = os.environ.get(var)
        if value is None:
            if default:
                os.environ[var] = default
                print(f"Setting {var} to default: {default}")
                env_vars[var] = default
            else:
                missing_vars.append(var)
        else:
            env_vars[var] = value

    if missing_vars:
        print("Error: Missing required environment variables:", file=sys.stderr)
        for var in missing_vars:
            print(f"  - {var}", file=sys.stderr)
        sys.exit(1)

    return env_vars


def setup_development_environment() -> None:
    """Set up the development environment."""
    try:
        # Ensure we're in the correct directory
        os.chdir(get_project_root())

        # Check Python version
        python_version = sys.version_info
        if python_version.major < 3 or (
            python_version.major == 3 and python_version.minor < 8
        ):
            print("Error: Python 3.8 or higher is required", file=sys.stderr)
            sys.exit(1)

    except Exception as e:
        print(f"Error setting up development environment: {str(e)}", file=sys.stderr)
        sys.exit(1)


def main() -> None:
    """Main entry point for the application."""
    # Set up environment
    setup_development_environment()
    env_vars: Dict[str, str] = check_environment()

    # Create and configure the application
    app: Flask = create_app()

    # Run the application
    app.run(
        host=env_vars["APP_HOST"],
        port=int(env_vars["APP_PORT"]),
        debug=(env_vars["FLASK_DEBUG"] == "1"),
    )


if __name__ == "__main__":
    main()
