from flask import Blueprint, request, jsonify
from peplab.backend.src.application.services.analysis.cheminformatics_service import CheminformaticsService

analysis_api_bp = Blueprint("analysis_api", __name__, url_prefix="/api/analysis")

@analysis_api_bp.route("/properties", methods=["POST"])
def get_properties():
    """
    Expects JSON payload: {"sequences": [["Ala", "Gly"], ...]}
    Returns the computed properties for each.
    """
    if not request.is_json:
        return jsonify({"error": "Request must be JSON"}), 400
        
    data = request.get_json()
    sequences = data.get("sequences", [])
    
    if not sequences:
        return jsonify({"error": "No sequences provided"}), 400
        
    service = CheminformaticsService()
    try:
        results = service.analyze_sequences(sequences)
        return jsonify({"results": results})
    except Exception as e:
        return jsonify({"error": str(e)}), 500
