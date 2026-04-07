import csv
import os
from peplab import create_app, db
from peplab.backend.src.domain.models.building_block import BuildingBlock
from peplab.backend.src.infrastructure.repositories.building_block_repository import BuildingBlockRepository

def seed_database():
    app = create_app()
    with app.app_context():
        repo = BuildingBlockRepository()
        
        csv_path = os.path.abspath("test_building_blocks.csv")
        if not os.path.exists(csv_path):
            print(f"CSV not found at {csv_path}")
            return
            
        with open(csv_path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f, skipinitialspace=True)
            for row in reader:
                name = row.get("name", "").strip()
                if not name:
                    continue
                
                bb = BuildingBlock(
                    name=name,
                    properties={
                        "alt_name1": row.get("alt_name1", "").strip(),
                        "alt_name2": row.get("alt_name2", "").strip(),
                        "position": row.get("position", "").strip()
                    },
                    metadata={
                        "smiles": row.get("smiles", "").strip()
                    }
                )
                
                try:
                    # check if already exists
                    repo.get_building_block_by_name(name)
                    print(f"Skipping {name}, already exists in DB.")
                except Exception:
                    repo.add_building_block(bb)
                    print(f"Added {name} to DB.")

if __name__ == "__main__":
    seed_database()
