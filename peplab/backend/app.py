'''from flask import Flask, request, jsonify
import sys
import os
from flask_cors import CORS

# Ensure project root is in sys.path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from peplab.design.library_design.composer import Composer
from peplab.design.library_design.combinatoric.combinative_composition import (
    CombinationComposition, PermutationComposition, CartesianProductComposition
)
from peplab.design.library_design.group_theoretic.grouptheoreticcomp import (
    CyclicPermutationComposition, DihedralPermutationComposition
)

app = Flask(__name__)
CORS(app)

# Helper function to validate input
def validate_input(data):
    required_fields = ["strategy", "input_data"]
    for field in required_fields:
        if field not in data:
            return False, f"Missing required field: {field}"
    if not isinstance(data["input_data"], list):
        return False, "'input_data' must be a list of items"
    return True, None

@app.route('/api/generate', methods=['POST'])  # Updated route to fit API structure
def generate_composition():
    data = request.json
    app.logger.info(f"Received request with data: {request.json}")

    # Validate input
    is_valid, error_message = validate_input(data)
    if not is_valid:
        return jsonify({"error": error_message}), 400

    strategy = data.get('strategy')
    input_data = data.get('input_data')
    r = data.get('r')  

    strategy_map = {
        "combinative": CombinationComposition,
        "permutative": PermutationComposition,
        "cartesian": CartesianProductComposition,
        "cyclic": CyclicPermutationComposition,
        "dihedral": DihedralPermutationComposition,
    }

    if strategy not in strategy_map:
        return jsonify({"error": "Invalid strategy"}), 400

    try:
        composer = Composer(strategy_map[strategy]())
        result = composer.generate_library(input_data, r=r) if r else composer.generate_library(input_data)

        # Return structured JSON
        return jsonify({"strategy": strategy, "input_data": input_data, "result": result})

    except Exception as e:
        app.logger.error(f"Error generating composition: {str(e)}")
        return jsonify({"error": str(e)}), 500

if __name__ == "__main__":
    app.run(debug=True, port=5001)  # Ensure it's on 5001
'''