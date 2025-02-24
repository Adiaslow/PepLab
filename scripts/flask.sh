#!/bin/bash
export FLASK_RUN_PORT=5001
export FLASK_RUN_HOST=0.0.0.0
flask "$@"
