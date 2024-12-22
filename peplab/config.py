# peplab/config.py

"""
Configuration settings for the PepLab application.

This module contains constants and configuration settings used throughout the PepLab application. All configurable settings are centralized here to ensure consistency and ease of maintenance.

Attributes:
    LOG_DIR (str): The directory where log files are stored.
    LOG_FILE (str): The path to the main log file.
"""

import os

# Directory for storing log files
LOG_DIR = os.path.join(os.path.dirname(__file__), '..', 'logs')

# Path to the main log file
LOG_FILE = os.path.join(LOG_DIR, 'peplab.log')

def create_log_directory():
    """
    Creates the log directory if it does not exist.

    This function ensures that the directory specified by LOG_DIR exists. If the directory does not exist, it is created. This is useful for ensuring that log files can be written to the specified location without encountering directory not found errors.
    """
    if not os.path.exists(LOG_DIR):
        os.makedirs(LOG_DIR)

# Create the log directory at module load time
create_log_directory()
