#!/bin/bash

# peplab/scripts/flask.sh
# Script to manage Flask operations

case "$1" in
    "run")
        export FLASK_APP=peplab.wsgi
        export FLASK_ENV=development
        export FLASK_DEBUG=1
        export FLASK_RUN_PORT=5001
        flask run
        ;;
    *)
        echo "Usage: $0 run"
        exit 1
        ;;
esac 