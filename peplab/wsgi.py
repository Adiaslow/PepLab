# peplab/wsgi.py
"""
WSGI server configuration for PepLab.

This module provides the entry point for running the application
with the correct host and port settings.
"""

import os
from peplab import create_app

app = create_app()

if __name__ == "__main__":
    port = int(os.getenv("APP_PORT", 5001))
    host = os.getenv("APP_HOST", "0.0.0.0")
    app.run(host=host, port=port, debug=True)
