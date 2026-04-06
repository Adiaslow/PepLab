import logging
from flask import Blueprint, request, jsonify

from peplab.backend.src.application.services.design.combinatoric.combination import Combination
from peplab.backend.src.application.services.design.combinatoric.permutation import Permutation
from peplab.backend.src.application.services.design.combinatoric.cyclic_permutation import CyclicPermutative
from peplab.backend.src.application.services.design.combinatoric.dihedral_permutation import DihedralPermutative
from peplab.backend.src.application.interfaces.design.composer import Composer

api_bp = Blueprint("api", __name__)
logger = logging.getLogger(__name__)

def validate_input(data):
    required_fields = ["strategy", "input_data"]
    for field in required_fields:
        if field not in data:
            return False, f"Missing required field: {field}"
    if not isinstance(data["input_data"], list):
        return False, "'input_data' must be a list of items"
    return True, None

@api_bp.route('/health', methods=['GET'])
def health_check():
    """Health check endpoint for the backend service."""
    return jsonify({"status": "ok"})

@api_bp.route('/design/generate', methods=['POST'])
def generate_composition():
    """
    API endpoint to generate peptide compositions based on selected strategy.
    """
    data = request.json
    logger.info(f"📩 Received API request: {data}")

    is_valid, error_message = validate_input(data)
    if not is_valid:
        return jsonify({"error": error_message}), 400

    strategy = data.get('strategy')
    input_data = data.get('input_data')
    r = data.get('r')

    strategy_map = {
        "combination": Combination,
        "permutation": Permutation,
        "cyclic": CyclicPermutative,
        "dihedral": DihedralPermutative,
    }

    if strategy not in strategy_map:
        return jsonify({"error": f"Invalid strategy '{strategy}'"}), 400

    try:
        composer = Composer(strategy_map[strategy]())
        result = composer.generate_library(input_data, length=r) if r is not None else composer.generate_library(input_data)

        return jsonify({
            "strategy": strategy,
            "input_data": input_data,
            "result": result
        })

    except Exception as e:
        logger.error(f"❌ Error in generate_composition: {str(e)}")
        return jsonify({"error": str(e)}), 500
