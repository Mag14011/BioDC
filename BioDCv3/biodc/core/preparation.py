# biodc/core/preparation.py

"""
Structure Preparation and Relaxation module for BioDC.
Coordinates the structure preparation workflow through multiple steps.
"""

import os
import subprocess
from pathlib import Path
from typing import Dict, Tuple, List, Optional
from dataclasses import dataclass

from biodc.utils.interaction import InteractionManager
from biodc.core.prep_modules import (
    initialize,
    select_disulfides,
    select_mutate,
    ligand_detection,
    select_ph_active_sites,
    create_res_indexing,
    process_residue_indexing,
    generate_tleap,
    validate_titratable_residues,
    struct_relax
)

from biodc.utils.interaction import InteractionManager

@dataclass
class PreparedStructure:
    """Container for structure preparation data."""
    pdb: str
    disulf_list: str = ""
    disulf_array: List = None
    sel_asp_ids: str = ""
    sel_glu_ids: str = ""
    sel_his_ids: str = ""
    sel_lys_ids: str = ""
    sel_tyr_ids: str = ""
    sel_prn_ids: str = ""

def reorder_structure(out_prefix: str) -> Tuple[str, str, str]:
    """Reorder residues using CPPTRAJ."""
    reorder_input = f"""
parm {out_prefix}.prmtop
trajin {out_prefix}.rst7
fixatomorder parmout {out_prefix}_reord.prmtop
trajout {out_prefix}_orig.pdb
trajout {out_prefix}_reord.pdb topresnum
trajout {out_prefix}_reord.rst7 topresnum
run
quit
    """

    with open("ReorderRes.in", "w") as f:
        f.write(reorder_input)

    subprocess.run(
        "cpptraj -i ReorderRes.in > ReorderRes.log 2> /dev/null",
        shell=True,
        check=True
    )

    return (
        f"{out_prefix}_reord.prmtop",
        f"{out_prefix}_orig.pdb",
        f"{out_prefix}_reord.pdb"
    )

def generate_cpin(out_prefix: str, reordered_prmtop: str, prep: 'PreparedStructure') -> str:
    """Generate cpin file for constant pH dynamics."""
    res_names = []
    res_ids = []

    if prep.sel_asp_ids:
        res_names.append('AS4')
        res_ids.append(prep.sel_asp_ids)
    if prep.sel_glu_ids:
        res_names.append('GL4')
        res_ids.append(prep.sel_glu_ids)
    if prep.sel_his_ids:
        res_names.append('HIP')
        res_ids.append(prep.sel_his_ids)
    if prep.sel_lys_ids:
        res_names.append('LYS')
        res_ids.append(prep.sel_lys_ids)
    if prep.sel_tyr_ids:
        res_names.append('TYR')
        res_ids.append(prep.sel_tyr_ids)
    if prep.sel_prn_ids:
        res_names.append('PRN')
        res_ids.append(prep.sel_prn_ids)

    if not res_names or not res_ids:
        return None

    cmd = [
        "cpinutil.py",
        "-resnames", " ".join(res_names),
        "-resnums", " ".join(res_ids),
        "-p", reordered_prmtop,
        "-igb", "2",
        "-op", f"{out_prefix}_new.prmtop",
        "-o", f"{out_prefix}.cpin"
    ]

    subprocess.run(" ".join(cmd), shell=True, check=True)
    subprocess.run(f"cp {out_prefix}_reord.rst7 {out_prefix}_new.rst7", shell=True, check=True)

    return f"{out_prefix}.cpin"
def structure_preparation_and_relaxation(
    launch_dir: Path,
    forcefield_dir: Path,
    struc_dir: Path,
    input_dict: Dict
) -> Tuple[str, str]:
    
    interaction_manager = InteractionManager(
        launch_dir=launch_dir,
        input_dict=input_dict
    )

def structure_preparation_and_relaxation(
    launch_dir: Path,
    forcefield_dir: Path,
    struc_dir: Path,
    input_dict: Dict
) -> Tuple[str, str]:
    """Main function for structure preparation and relaxation."""

    interaction_manager = InteractionManager(
        launch_dir=launch_dir,
        input_dict=input_dict
    )

    print("""
 We will ask a series of questions about your structure
 in order to properly prepare it.
    """)

    # Initialize structure
    pdb = initialize(launch_dir, interaction_manager.get_input_dict())
#   print(f"Debug - pdb {pdb}")
     
    prep = PreparedStructure(pdb=pdb)
    has_titratable = False

    # Handle disulfide selection
#   print("Debug before prompt - input_dict keys:", input_dict.keys())
    if interaction_manager.yes_no_prompt(
        "SelDisulfides",
        "\nAre there disulfide linkages in your structure?"
    ):
#       print("Debug after prompt - input_dict keys:", input_dict.keys())
        prep.disulf_list, prep.disulf_array = select_disulfides(
            prep.pdb,  # Pass the PDB name directly
            interaction_manager.get_input_dict(),
            launch_dir
        )
    else:
        print(" No disulfide linkages will be present in the prepared structure.")

    # Handle mutations
    if interaction_manager.yes_no_prompt(
        "ChooseMut",
        "\nWould you like to mutate a residue?"
    ):
        prep.pdb = select_mutate(prep.pdb, launch_dir, interaction_manager.get_input_dict())
    else:
        print(" Structure will not be mutated.")

    # Handle pH active sites
    if interaction_manager.yes_no_prompt(
        "SelCpH",
        "\nDo you intend to run molecular dynamics where\n"
        " pH active residues are titrated?"
    ):
        has_titratable = True
        current_dict = interaction_manager.get_input_dict()

        # Check for automatic selection of all titratable residues
        select_all = current_dict.get("SelectAllTitratable", "").lower() in ['yes', 'y', 'true', '1']

        if select_all:
            print("\n Automatically selecting all available titratable residues...")
            # We'll pass a special flag to select_ph_active_sites
            (
                prep.sel_asp_ids, prep.sel_glu_ids, prep.sel_his_ids,
                prep.sel_lys_ids, prep.sel_tyr_ids, prep.sel_prn_ids
            ) = select_ph_active_sites(
                prep.pdb, "FirstPass", {"SelectAllTitratable": "yes"}, launch_dir
            )
        elif all(k in current_dict for k in ["SelASPIDs_1", "SelGLUIDs_1", "SelHISIDs_1"]):
            prep.sel_asp_ids = current_dict["SelASPIDs_1"]
            prep.sel_glu_ids = current_dict["SelGLUIDs_1"]
            prep.sel_his_ids = current_dict["SelHISIDs_1"]
        else:
            (
                prep.sel_asp_ids, prep.sel_glu_ids, prep.sel_his_ids,
                prep.sel_lys_ids, prep.sel_tyr_ids, prep.sel_prn_ids
            ) = select_ph_active_sites(
                prep.pdb, "FirstPass", current_dict, launch_dir
            )
    else:
        print(""" No residues in the prepared structure will be titratable
 if you decide to run molecular dynamics.""")
        prep.sel_asp_ids = ""
        prep.sel_glu_ids = ""
        prep.sel_his_ids = ""
        prep.sel_lys_ids = ""
        prep.sel_tyr_ids = ""
        prep.sel_prn_ids = ""

    # Process the structure
#   print(f"file = {prep.pdb}")
    prep.pdb = create_res_indexing(
        prep.pdb, interaction_manager.get_input_dict(), launch_dir
    )
#   print(f"file = {prep.pdb}")

    # Process the structure and get required info for tleap
#   print(f"file = {prep.pdb}")
    structure_info, indexing_data = process_residue_indexing(
        prep.pdb, prep.disulf_list,
        prep.sel_asp_ids, prep.sel_glu_ids, prep.sel_his_ids,
        prep.sel_lys_ids, prep.sel_tyr_ids, prep.sel_prn_ids,
        interaction_manager.get_input_dict(), launch_dir
    )
#   print(f"file = {prep.pdb}")
    
    # Generate tleap input and run
    out_prefix, solv_env = generate_tleap.generate_tleap_input(
        prep.pdb,
        forcefield_dir, 
        structure_info,
        indexing_data,
        interaction_manager.get_input_dict(),
        launch_dir
    )

    # Reorder the structure
    print("\n")
    print("=" * 60)
    print("Reordering residues for consistency with molecular dynamics requirements...")
    reordered_prmtop, original_pdb, reordered_pdb = reorder_structure(out_prefix)
    final_prmtop = reordered_prmtop
    print(f" Created reordered topology: {reordered_prmtop}")
    print("=" * 60)

    # Validate titratable residues after reordering
    if has_titratable:
        print("\n")
        print("=" * 60)
        print("Validating titratable residue selections after reordering...")
        prep = validate_titratable_residues(
            reordered_pdb,
            prep,
            interaction_manager,
            launch_dir
        )

    # Generate cpin file if needed
    if has_titratable:
        print("\n")
        print("=" * 60)
        print("Generating cpin file for constant pH dynamics...")
        print("=" * 60)
        cpin_file = generate_cpin(out_prefix, reordered_prmtop, prep)
        final_prmtop = f"{out_prefix}_new.prmtop"
        print(f" Created cpin file: {cpin_file}")
        print(f" Created modified topology for constant pH: {final_prmtop}")

    # Run structure relaxation
    print("\n")
    print("=" * 60)
    print("Structure Relaxation")
    print("=" * 60)
    struct_relax(
        launch_dir, struc_dir, final_prmtop.removesuffix('.prmtop'), solv_env, interaction_manager.get_input_dict()
    )

    return out_prefix, solv_env

def run(
    launch_dir: Path, 
    forcefield_dir: Path, 
    struc_dir: Path,
    input_dict: Dict
) -> Tuple[str, str]:
    """Wrapper function for structure preparation and relaxation."""
    return structure_preparation_and_relaxation(
        launch_dir,
        forcefield_dir,
        struc_dir,
        input_dict
    )
