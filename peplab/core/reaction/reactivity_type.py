from enum import Enum

class ReactivityType(Enum):
    NUCLEOPHILE = "nucleophile"
    ELECTROPHILE = "electrophile"
    LEAVING_GROUP = "leaving_group"
    ACID = "acid"
    BASE = "base"
    RADICAL = "radical"
