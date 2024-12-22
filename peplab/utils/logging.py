# peplab/utils/logging.py

"""
Logging Utility Module

This module provides a centralized logging utility for the application. It includes a logger class that implements the Singleton pattern to ensure that only one instance of the logger exists throughout the application.

Classes:
    Logger: A singleton logger that provides centralized logging with file and console handlers.

Usage:
    from peplab.utils.logging import Logger
    logger = Logger().get_logger()
    logger.info('This is an info message')
    logger.debug('This is a debug message')
"""

import logging
import os
from logging.handlers import RotatingFileHandler

class Logger:
    """Centralized Logger

    This class implements a centralized logger for the application, ensuring that all logging is done through a single instance. It provides both file and console logging with rotating files to manage log sizes.

    Design Patterns:
        Singleton: Ensures that only one instance of the logger exists throughout the application. This is important for maintaining a consistent logging configuration and avoiding the overhead of multiple logger instances.

    Attributes:
        logger: The logger instance used to log messages.
    """

    _instance = None

    def __new__(cls, *args, **kwargs):
        if not cls._instance:
            cls._instance = super(Logger, cls).__new__(cls, *args, **kwargs)
            cls._instance._initialize()
        return cls._instance

    def _initialize(self):
        log_dir = os.path.join(os.path.dirname(__file__), '..', '..', 'logs')
        os.makedirs(log_dir, exist_ok=True)
        log_file = os.path.join(log_dir, 'peplab.log')

        self.logger = logging.getLogger('peplab')
        self.logger.setLevel(logging.DEBUG)

        # Create a rotating file handler
        file_handler = RotatingFileHandler(log_file, maxBytes=5 * 1024 * 1024, backupCount=5)
        file_handler.setLevel(logging.DEBUG)

        # Create a console handler
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.INFO)

        # Create a formatter and set it for both handlers
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        file_handler.setFormatter(formatter)
        console_handler.setFormatter(formatter)

        # Add the handlers to the logger
        self.logger.addHandler(file_handler)
        self.logger.addHandler(console_handler)

    def get_logger(self):
        """Returns the logger instance."""
        return self.logger

# Usage example:
# logger = Logger().get_logger()
# logger.info('This is an info message')
# logger.debug('This is a debug message')
