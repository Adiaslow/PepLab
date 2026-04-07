import json
from typing import List, Dict, Any
from rdkit import Chem
from rdkit.Chem import Descriptors
from peplab.backend.src.infrastructure.repositories.building_block_repository import BuildingBlockRepository

class CheminformaticsService:
    def __init__(self):
        self.repo = BuildingBlockRepository()

    def analyze_sequences(self, sequences: List[List[str]]) -> List[Dict[str, Any]]:
        """
        Analyzes a list of sequences (where each sequence is a list of building block names)
        and computes RDKit properties.
        """
        results = []
        
        for seq in sequences:
            # Attempt to get the 1-letter code for each building block
            one_letter_codes = []
            is_valid = True
            
            for bb_name in seq:
                try:
                    bb = self.repo.get_building_block_by_name(bb_name)
                    # alt_name2 holds the 1 letter code
                    code = bb.properties.get("alt_name2", "") if bb.properties else ""
                    if not code or len(code) != 1:
                        is_valid = False
                        break
                    one_letter_codes.append(code)
                except Exception:
                    is_valid = False
                    break
                    
            res = {
                "sequence": seq,
                "length": len(seq)
            }
            
            if is_valid:
                fasta_str = "".join(one_letter_codes)
                # Build mol from sequence
                mol = Chem.MolFromSequence(fasta_str)
                if mol:
                    res["molecular_weight"] = round(Descriptors.MolWt(mol), 2)
                    res["exact_mass"] = round(Descriptors.ExactMolWt(mol), 4)
                    res["log_p"] = round(Descriptors.MolLogP(mol), 2)
                    res["tpsa"] = round(Descriptors.TPSA(mol), 2)
                    res["h_donors"] = Descriptors.NumHDonors(mol)
                    res["h_acceptors"] = Descriptors.NumHAcceptors(mol)
                    res["status"] = "Success"
                else:
                    res["status"] = "Failed to build mol"
                    self._fill_nulls(res)
            else:
                res["status"] = "Non-canonical blocks unsupported"
                self._fill_nulls(res)
                
            results.append(res)
            
        return results

    def _fill_nulls(self, res: dict):
        res["molecular_weight"] = None
        res["exact_mass"] = None
        res["log_p"] = None
        res["tpsa"] = None
        res["h_donors"] = None
        res["h_acceptors"] = None
