from flask import Blueprint, request, jsonify
from peplab.backend.src.application.services.design.combinatoric.permutation import generate_permutations

combinatorial_routes = Blueprint("combinatorial_routes", __name__)

@combinatorial_routes.route("/api/combinatorial/permutations", methods=["POST"])
def generate_permutation():
    data = request.json
    sequence = data.get("sequence")
    num_permutations = data.get("num_permutations")

    if not sequence or not num_permutations:
        return jsonify({"error": "Missing required parameters"}), 400

    result = generate_permutations(sequence, num_permutations)
    return jsonify({"permutations": result})