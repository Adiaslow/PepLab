from flask import Flask, request, jsonify
import sys
import os
from flask_cors import CORS  # Import Flask-CORS

# Add project root to sys.path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from peplab.design.library_design.composer import Composer
from peplab.design.library_design.combinatoric.combinative_composition import (
    CombinationComposition, PermutationComposition, CartesianProductComposition
)
from peplab.design.library_design.group_theoretic.grouptheoreticcomp import (
    CyclicPermutationComposition, DihedralPermutationComposition
)

app = Flask(__name__)
CORS(app, resources={r"/*": {"origins": ["http://127.0.0.1:8000", "http://localhost:8000"]}})

# Helper function to validate input
def validate_input(data):
    required_fields = ["strategy", "input_data"]
    for field in required_fields:
        if field not in data:
            return False, f"Missing required field: {field}"
    if not isinstance(data["input_data"], list):
        return False, "'input_data' must be a list of items"
    return True, None
@app.route('/generate', methods=['POST'])
def generate_composition():
    data = request.json
    app.logger.info(f"Received request with data: {request.json}")

    # Validate input
    is_valid, error_message = validate_input(data)
    if not is_valid:
        return jsonify({"testerror": error_message}), 400

    strategy = data.get('strategy')  # e.g., "combinative", "cyclic", etc.
    input_data = data.get('input_data')  # List of numbers or items
    r = data.get('r')  # Optional length for combinations/permutations

    # Map strategy to the corresponding composition class
    strategy_map = {
        "combinative": CombinationComposition,
        "permutative": PermutationComposition,
        "cartesian": CartesianProductComposition,
        "cyclic": CyclicPermutationComposition,
        "dihedral": DihedralPermutationComposition,
    }

    # Check if strategy is valid
    if strategy not in strategy_map:
        return jsonify({"error": "Invalid strategy"}), 400

    try:
        # Initialize the appropriate composition strategy
        composer = Composer(strategy_map[strategy]())

        # Generate the library
        if r:
            result = composer.generate_library(input_data, r=r)
        else:
            result = composer.generate_library(input_data)

        # Flatten result for combination strategies
        formatted_result = [tuple(seq) if isinstance(seq, list) else (seq,) for seq in result]

        return "\n".join([f"{pair[0]} {pair[1]}" for pair in formatted_result])


    except Exception as e:
        # Log and return the error
        app.logger.error(f"Error generating composition: {str(e)}")
        return jsonify({"testerror": str(e)}), 500

if __name__ == "__main__":
    app.run(debug=True)
