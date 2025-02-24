# Setup script for Windows PowerShell
$env:FLASK_APP="wsgi.py"
$env:FLASK_ENV="development"
$env:FLASK_DEBUG="1"
$env:APP_PORT="5001"
$env:APP_HOST="localhost"

Write-Host "Environment variables set. You can now run:"
Write-Host "flask run" 