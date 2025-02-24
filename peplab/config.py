"""
Configuration settings for the PepLab application.
"""

import os
from pathlib import Path


class Config:
    """Base configuration class."""

    # Get the project root directory
    BASE_DIR = Path(__file__).resolve().parent

    # Flask settings
    SECRET_KEY = os.getenv("SECRET_KEY", "527-436-828")

    # Database settings
    SQLALCHEMY_DATABASE_URI = os.getenv(
        "DATABASE_URL", "sqlite:///" + str(BASE_DIR / "peplab.db")
    )
    SQLALCHEMY_TRACK_MODIFICATIONS = False

    # Security settings
    SESSION_COOKIE_SECURE = True
    REMEMBER_COOKIE_SECURE = True
    SESSION_COOKIE_HTTPONLY = True
    REMEMBER_COOKIE_HTTPONLY = True
