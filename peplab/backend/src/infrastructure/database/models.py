from peplab import db
import uuid

# Association table for Peptide to BuildingBlock many-to-many relationship
peptide_building_block = db.Table(
    'peptide_building_block',
    db.Column('peptide_id', db.String(36), db.ForeignKey('peptides.id'), primary_key=True),
    db.Column('building_block_id', db.String(36), db.ForeignKey('building_blocks.id'), primary_key=True),
    db.Column('sequence_order', db.Integer, nullable=False, default=0)
)

class BuildingBlockModel(db.Model):
    __tablename__ = 'building_blocks'
    
    id = db.Column(db.String(36), primary_key=True, default=lambda: str(uuid.uuid4()))
    name = db.Column(db.String(255), unique=True, nullable=False)
    embeddings = db.Column(db.JSON, nullable=True)
    encodings = db.Column(db.JSON, nullable=True)
    properties = db.Column(db.JSON, nullable=True)
    metadata_json = db.Column('metadata', db.JSON, nullable=True)

class PeptideModel(db.Model):
    __tablename__ = 'peptides'
    
    id = db.Column(db.String(36), primary_key=True, default=lambda: str(uuid.uuid4()))
    name = db.Column(db.String(255), nullable=True)
    embeddings = db.Column(db.JSON, nullable=True)
    encodings = db.Column(db.JSON, nullable=True)
    properties = db.Column(db.JSON, nullable=True)
    metadata_json = db.Column('metadata', db.JSON, nullable=True)
    library_id = db.Column(db.String(36), db.ForeignKey('libraries.id'), nullable=True)
    
    building_blocks = db.relationship('BuildingBlockModel', secondary=peptide_building_block, lazy='select')

class LibraryModel(db.Model):
    __tablename__ = 'libraries'
    
    id = db.Column(db.String(36), primary_key=True, default=lambda: str(uuid.uuid4()))
    name = db.Column(db.String(255), nullable=False)
    description = db.Column(db.Text, nullable=True)
    metadata_json = db.Column('metadata', db.JSON, nullable=True)
    
    peptides = db.relationship('PeptideModel', backref='library', lazy='select')
