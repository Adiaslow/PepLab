from flask import Blueprint, jsonify
from rdkit import Chem
from rdkit.Chem import Draw
import pandas as pd
import base64
import io

bp = Blueprint("building_blocks", __name__)


@bp.route("/api/building-blocks", methods=["GET"])
def get_building_blocks():
    try:
        # Read the CSV file
        df = pd.read_csv("test_building_blocks.csv")
        building_blocks = []

        for _, row in df.iterrows():
            # Create RDKit molecule from SMILES
            mol = Chem.MolFromSmiles(row["smiles"])

            # Generate 2D depiction
            img = Draw.MolToImage(mol)

            # Convert image to base64
            img_buffer = io.BytesIO()
            img.save(img_buffer, format="PNG")
            img_str = base64.b64encode(img_buffer.getvalue()).decode()

            # Create building block object
            building_block = {
                "name": row["name"],
                "alt_name1": row["alt_name1"],
                "alt_name2": row["alt_name2"],
                "smiles": row["smiles"],
                "position": row["position"],
                "image": f"data:image/png;base64,{img_str}",
            }
            building_blocks.append(building_block)

        return jsonify({"building_blocks": building_blocks})
    except Exception as e:
        return jsonify({"error": str(e)}), 500
