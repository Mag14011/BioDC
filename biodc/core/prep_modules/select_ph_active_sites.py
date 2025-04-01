# biodc/core/preparation/select_ph_active_sites.py

"""
Module for selecting pH-active sites for constant pH molecular dynamics,
excluding residues that are coordinated to metal centers.
"""

from pathlib import Path
from typing import Dict, List, Tuple
import logging
from dataclasses import dataclass
from Bio import PDB

from biodc.core.prep_modules.ligand_detection import LigandDetector
from biodc.utils.interaction import InteractionManager

logger = logging.getLogger(__name__)

@dataclass
class PHActiveSites:
    """Container for pH-active site selections."""
    asp_ids: str = ''
    glu_ids: str = ''
    his_ids: str = ''
    lys_ids: str = ''
    tyr_ids: str = ''
    prn_ids: str = ''

class ResidueSelector:
    def __init__(self, pdb: str, launch_dir: Path):
        self.pdb = pdb
        self.launch_dir = launch_dir
        self.pdb_path = Path(f"{pdb}.pdb")

    def _calculate_distance(self, coord1, coord2):
        return ((coord1[0] - coord2[0])**2 +
                (coord1[1] - coord2[1])**2 +
                (coord1[2] - coord2[2])**2)**0.5

    def _is_terminus(self, res_id: int, ca_data: dict, max_distance: float = 8.0) -> bool:
        current_coords = ca_data[res_id]["coords"]
        has_prev = False
        has_next = False

        for other_id, other_data in ca_data.items():
            if other_id == res_id:
                continue
            
            dist = self._calculate_distance(current_coords, other_data["coords"])
            if dist > max_distance:
                continue
            
            if other_id > res_id:
                if other_data["resname"] == "HEM":
                    return True
                has_next = True
            else:
                has_prev = True
            
        return not (has_prev and has_next)

    def find_titratable_residues(self, res_name: str) -> List[str]:
        residue_ids = set()
        terminal_ids = set()
        ca_data = {}

#       print(f"\nDEBUG: Looking for {res_name} in {self.pdb_path}")

        with open(self.pdb_path) as f:
            for line in f:
                if line.startswith("ATOM  ") and "CA " in line[12:16]:
                    res_name_pdb = line[17:20].strip()
                    res_id = int(line[22:26])
#                   print(f"DEBUG: Found CA atom {res_name_pdb} {res_id}")

                    # Store coordinates for all residues
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    ca_data[res_id] = {
                        "coords": (x, y, z),
                        "resname": res_name_pdb
                    }

                    if res_name_pdb == res_name:
                        residue_ids.add(res_id)
#                       print(f"DEBUG: Found titratable {res_name} {res_id}")

                # Check for hetatm CA atoms
                elif line.startswith("HETATM") and "CA " in line[12:16]:
                    res_name_pdb = line[17:20].strip()
                    if res_name_pdb == "HEM":
                        res_id = int(line[22:26])
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        ca_data[res_id] = {
                            "coords": (x, y, z),
                            "resname": "HEM"
                        }
#                       print(f"DEBUG: Found HEM CA at {res_id}")

                if (f"OXT {res_name}" in line) or (line.startswith("TER") and line[17:20].strip() == res_name and res_name != "PRN"):
                    term_id = int(line[22:26])
                    terminal_ids.add(term_id)
#                   print(f"DEBUG: Found explicit terminal {res_name} {term_id}")

#       print(f"DEBUG: Found residue_ids: {residue_ids}")
#       print(f"DEBUG: Found ca_data keys: {sorted(list(ca_data.keys()))}")

        # Distance and heme-based terminal detection
        for res_id in residue_ids:
            if res_id not in terminal_ids and res_name != "PRN":
                if self._is_terminus(res_id, ca_data):
                    terminal_ids.add(res_id)
#                   print(f"DEBUG: Distance/heme analysis identified {res_name} {res_id} as terminal")

        valid_residues = sorted([str(rid) for rid in residue_ids if rid not in terminal_ids], key=int)
#       print(f"DEBUG: Terminal IDs: {terminal_ids}")
#       print(f"DEBUG: Final {res_name} residues: {valid_residues}")
        return valid_residues

def get_residue_selection(
    interaction_manager: InteractionManager,
    res_type: str,
    residlist: List[str],
    select_all: bool = False,
    key_suffix: str = "1"
) -> str:
    if not residlist:
        return ''

    # Sort the residlist numerically before printing
    sorted_residlist = sorted(residlist, key=int)
    print(f"\n There are {len(sorted_residlist)} {res_type} residues that can be titrated with IDs: ", end=" ")
    print(*sorted_residlist)

    if select_all:
        selection = " ".join(residlist)
        interaction_manager._record_interaction(f"Sel{res_type}IDs_{key_suffix}", selection)
        print(f"  Automatically selected all {res_type} residues: {selection}")
        return selection

    selection = interaction_manager.prompt(
        f"Sel{res_type}IDs_{key_suffix}",
        f"  pH active {res_type} residue IDs",
        allow_empty=False
    )

    return selection

def select_ph_active_sites(
    pdb: str,
    switch: str,
    input_dict: Dict,
    launch_dir: Path
) -> Tuple[str, str, str, str, str, str]:
    interaction_manager = InteractionManager(
        launch_dir=launch_dir,
        input_dict=input_dict
    )

    select_all = input_dict.get("SelectAllTitratable", "").lower() in ['yes', 'y', 'true', '1']

    print("""
 Residues ASP, GLU, LYS, HIS, TYR, and PRN can be titrated.
 Analyzing structure for metal-coordinated residues...
    """)

    detector = LigandDetector(f"{pdb}.pdb")
    coordinated = detector.find_metal_coordinations()
    detector.print_coordination_info()

    excluded = detector.get_excluded_residues()

    if any(excluded.values()):
        print("\n* The following metal-coordinated residues will be excluded from titration:")
        for res_type, ids in excluded.items():
            if ids:
                print(f" {res_type}: {', '.join(map(str, sorted(ids)))}")

    if not select_all:
        print(f"""
 Please note: N- or C-terminal residues cannot be titrated
 because the forcefield parameters are not available.
 C-terminal residue IDs will ONLY be excluded from the
 presented lists if an OXT atom type is present in the
 PDB file. N-terminal residues will not be automatically
 excluded, so be careful not to select them.
        """)

    print(f" Reading {pdb}.pdb ...")

    selector = ResidueSelector(pdb, launch_dir)
    sites = PHActiveSites()

    residue_types = ["ASP", "GLU", "HIS", "TYR", "LYS", "PRN"]

    for res_type in residue_types:
        residlist = selector.find_titratable_residues(res_type)
        if res_type in excluded:
            residlist = [r for r in residlist if int(r) not in excluded[res_type]]

        if res_type == "ASP":
            sites.asp_ids = get_residue_selection(interaction_manager, res_type, residlist, select_all)
        elif res_type == "GLU":
            sites.glu_ids = get_residue_selection(interaction_manager, res_type, residlist, select_all)
        elif res_type == "HIS":
            sites.his_ids = get_residue_selection(interaction_manager, res_type, residlist, select_all)
        elif res_type == "LYS":
            sites.lys_ids = get_residue_selection(interaction_manager, res_type, residlist, select_all)
        elif res_type == "TYR":
            sites.tyr_ids = get_residue_selection(interaction_manager, res_type, residlist, select_all)
        elif res_type == "PRN":
            sites.prn_ids = get_residue_selection(interaction_manager, res_type, residlist, select_all)

    return (sites.asp_ids, sites.glu_ids, sites.his_ids,
            sites.lys_ids, sites.tyr_ids, sites.prn_ids)

def validate_titratable_residues(
    reordered_pdb: str,
    prep: 'PreparedStructure',
    interaction_manager: 'InteractionManager',
    launch_dir: Path
) -> 'PreparedStructure':
    """
    Validate titratable residue selections after structure reordering.
    Also handles PRN residue selection which only exists after processing.
    """
#   print(f"\nDEBUG: Using PDB file: {reordered_pdb}")
    selector = ResidueSelector(reordered_pdb.removesuffix('.pdb'), launch_dir)
    
    # Check if automatic selection is enabled
    input_dict = interaction_manager.get_input_dict()
    select_all = input_dict.get("SelectAllTitratable", "").lower() in ['yes', 'y', 'true', '1']
    
    # Dictionary mapping processed residue names to their IDs in prep
    residue_map = {
        "AS4": prep.sel_asp_ids,
        "GL4": prep.sel_glu_ids,
        "HIP": prep.sel_his_ids,
        "LYS": prep.sel_lys_ids,
        "TYR": prep.sel_tyr_ids
    }
    
#   print("\nDEBUG: Previously selected residues:")
    for res_type, selected_ids in residue_map.items():
        print(f"{res_type}: {selected_ids}")

    # Validate each residue type
    for res_type, selected_ids in residue_map.items():
        if not selected_ids:  # Skip if none were selected
            continue
            
        current_ids = selector.find_titratable_residues(res_type)
#       print(f"\nDEBUG: Found {res_type} residues: {current_ids}")
        selected_set = set(selected_ids.split())
#       print(f"DEBUG: Selected {res_type} residues: {selected_set}")
        
        # Check if any selected residues are no longer valid
        invalid_residues = [rid for rid in selected_set if rid not in current_ids]
        if invalid_residues:
            print(f"\n Warning: Some previously selected {res_type} residues are no longer valid: {' '.join(invalid_residues)}")
            print(f" Currently available {res_type} residues: {' '.join(current_ids)}")
            
            # Re-prompt for selection
            new_selection = get_residue_selection(
                interaction_manager,
                res_type,
                current_ids,
                select_all=select_all,  # Pass the select_all flag here
                key_suffix="2"
            )
            
            # Update the corresponding field in prep
            if res_type == "AS4":
                prep.sel_asp_ids = new_selection
            elif res_type == "GL4":
                prep.sel_glu_ids = new_selection
            elif res_type == "HIP":
                prep.sel_his_ids = new_selection
            elif res_type == "LYS":
                prep.sel_lys_ids = new_selection
            elif res_type == "TYR":
                prep.sel_tyr_ids = new_selection

    # Handle PRN (heme propionic acid) selection
    prn_ids = selector.find_titratable_residues("PRN")
    if prn_ids:
        print("\n PRN (heme propionic acid) residues are now available for titration.")
        prep.sel_prn_ids = get_residue_selection(
            interaction_manager,
            "PRN",
            prn_ids,
            select_all=select_all,  # Use the select_all flag here too
            key_suffix="2"
        )
    print("=" * 60)

    return prep
