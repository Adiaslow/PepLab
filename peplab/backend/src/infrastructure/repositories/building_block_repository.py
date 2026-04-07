# peplab/backend/src/infrastructure/repositories/building_block_repository.py

# External imports
from typing import List, Dict

# Internal imports
from peplab.backend.src.domain.models.building_block import BuildingBlock
from peplab.backend.src.core.exceptions.domain_exceptions import (
    InvalidBuildingBlockError,
)
from peplab.backend.src.core.exceptions.infrastructure_exceptions import (
    EmptyBuildingBlockRepositoryError,
)
from peplab.backend.src.infrastructure.database.models import BuildingBlockModel
from peplab import db


class BuildingBlockRepository:
    """Repository for building blocks.
    
    Uses SQLAlchemy for persistence while returning Domain models.
    """
    
    def _to_domain(self, model: BuildingBlockModel) -> BuildingBlock:
        return BuildingBlock(
            id=model.id,
            name=model.name,
            embeddings=model.embeddings,
            encodings=model.encodings,
            properties=model.properties,
            metadata=model.metadata_json
        )

    def _to_model(self, domain: BuildingBlock) -> BuildingBlockModel:
        return BuildingBlockModel(
            id=str(domain.id) if domain.id else None,
            name=domain.name,
            embeddings=domain.embeddings,
            encodings=domain.encodings,
            properties=domain.properties,
            metadata_json=domain.metadata
        )

    def index_building_blocks(self) -> None:
        """Index the building blocks in the repository."""
        pass

    def add_building_block(self, building_block: BuildingBlock) -> None:
        """Add a building block to the repository."""
        model = self._to_model(building_block)
        db.session.add(model)
        db.session.commit()

    def get_building_block_by_id(self, id: str) -> BuildingBlock:
        """Get a building block from the repository by id."""
        model = db.session.get(BuildingBlockModel, id)
        if not model:
            raise InvalidBuildingBlockError(f"Building block with id {id} not found.")
        return self._to_domain(model)

    def get_building_block_by_name(self, name: str) -> BuildingBlock:
        """Get a building block from the repository by name."""
        model = db.session.query(BuildingBlockModel).filter_by(name=name).first()
        if not model:
            raise InvalidBuildingBlockError(f"Building block with name {name} not found.")
        return self._to_domain(model)

    def get_all_building_blocks(self) -> List[BuildingBlock]:
        """Get all building blocks from the repository."""
        models = db.session.query(BuildingBlockModel).all()
        if not models:
            raise EmptyBuildingBlockRepositoryError("No building blocks in the repository.")
        return [self._to_domain(m) for m in models]

    def remove_building_block_by_id(self, id: str) -> None:
        """Remove a building block from the repository by id."""
        model = db.session.get(BuildingBlockModel, id)
        if not model:
            raise InvalidBuildingBlockError(f"Building block with id {id} not found.")
        db.session.delete(model)
        db.session.commit()

    def remove_building_block_by_name(self, name: str) -> None:
        """Remove a building block from the repository by name."""
        model = db.session.query(BuildingBlockModel).filter_by(name=name).first()
        if not model:
            raise InvalidBuildingBlockError(f"Building block with name {name} not found.")
        db.session.delete(model)
        db.session.commit()
