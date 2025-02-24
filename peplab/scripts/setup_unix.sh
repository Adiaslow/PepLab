#!/bin/bash

# Setup script for Unix-like systems (macOS/Linux)

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Get the absolute path to the project root
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

echo -e "${YELLOW}Setting up PepLab development environment...${NC}"
echo -e "${YELLOW}Project root: ${PROJECT_ROOT}${NC}"

# Check if Python 3.8+ is installed
if ! command -v python3 &> /dev/null; then
    echo -e "${RED}Python 3 not found. Please install Python 3.8 or higher.${NC}"
    exit 1
fi

# Check for required packages
echo -e "${YELLOW}Checking required packages...${NC}"
if ! python3 -c "import flask" 2>/dev/null; then
    echo -e "${YELLOW}Installing Flask...${NC}"
    pip install flask
fi

if ! python3 -c "import dotenv" 2>/dev/null; then
    echo -e "${YELLOW}Installing python-dotenv...${NC}"
    pip install python-dotenv
fi

if ! python3 -c "import flask_cors" 2>/dev/null; then
    echo -e "${YELLOW}Installing Flask-CORS...${NC}"
    pip install flask-cors
fi

# Check if port 5000 is in use
if lsof -Pi :5000 -sTCP:LISTEN -t >/dev/null ; then
    echo "Notice: Port 5000 is in use. Using port 5001 instead."
    export FLASK_RUN_PORT=5001
fi

# Set environment variables
export FLASK_APP=peplab.wsgi
export FLASK_ENV=development
export FLASK_DEBUG=1
export PYTHONPATH="${PROJECT_ROOT}:${PYTHONPATH}"

echo -e "${GREEN}Environment setup complete!${NC}"
echo -e "${GREEN}You can now run:${NC}"
echo -e "${YELLOW}flask run${NC}" 