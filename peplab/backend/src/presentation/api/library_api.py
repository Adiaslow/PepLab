from flask import Blueprint, request, jsonify
from peplab import db
from peplab.backend.src.infrastructure.database.models import LibraryModel, PeptideModel, BuildingBlockModel
import uuid

library_api_bp = Blueprint("library_api", __name__, url_prefix="/api/library")

@library_api_bp.route("/save", methods=["POST"])
def save_library():
    """
    Expects JSON payload: { "name": "Generated Library", "sequences": [["Ala", "Gly"], ...] }
    """
    data = request.get_json()
    sequences = data.get("sequences", [])
    name = data.get("name", "Unnamed Library")
    
    if not sequences:
        return jsonify({"error": "No sequences provided."}), 400
        
    try:
        library = LibraryModel(name=name)
        db.session.add(library)
        
        # Get all unique block names in the sequences
        unique_block_names = set(block for seq in sequences for block in seq)
        blocks = db.session.query(BuildingBlockModel).filter(BuildingBlockModel.name.in_(unique_block_names)).all()
        block_map = {block.name: block for block in blocks}
        
        # Record keeping
        for seq in sequences:
            seq_name = " ".join(seq)
            peptide = PeptideModel(name=seq_name, library=library)
            
            # Map constituents uniquely to avoid primary key collisions on the assoc table
            unique_constituents = set(seq)
            for bb_name in unique_constituents:
                if bb_name in block_map:
                    peptide.building_blocks.append(block_map[bb_name])
            
            db.session.add(peptide)
            
        db.session.commit()
        return jsonify({"success": True, "library_id": library.id, "peptide_count": len(sequences)})
    except Exception as e:
        db.session.rollback()
        return jsonify({"error": str(e)}), 500
