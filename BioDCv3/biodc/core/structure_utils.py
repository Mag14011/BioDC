"""
Utility functions for structure preparation and handling in BioDC.
"""

import sys
from pathlib import Path
import subprocess

from rich.console import Console

from biodc.utils.interaction import InteractionManager

# Initialize console for rich output
console = Console()

def get_prepared_structure(launch_dir: Path, interaction_manager) -> str:
    """
    Prepare structure for energetic evaluation.
    
    Args:
        launch_dir: Project launch directory
        interaction_manager: InteractionManager for user interactions
    
    Returns:
        Path to PDB file for analysis
    
    Raises:
        SystemExit if structure preparation is incomplete
    """
#   # First, offer manual PDB selection
#   manual_select = interaction_manager.yes_no_prompt(
#       "manual_pdb_select",
#       "\nWould you like to manually select a PDB file for analysis?"
#   )
#   
#   if manual_select:
#       manual_pdb = interaction_manager.prompt(
#           "manual_pdb_path",
#           "Enter the full path to the PDB file you want to analyze:",
#           input_type=str
#       )
#       
#       pdb_path = Path(manual_pdb)
#       if not pdb_path.exists():
#           console.print("[bold red]Error:[/] Specified PDB file does not exist.")
#           sys.exit(1)
#       
#       return str(pdb_path)
    
    # Check SPR directory exists
    spr_dir = launch_dir / "SPR"
    if not spr_dir.exists():
        console.print("[bold red]Error:[/] Structure Preparation directory not found.")
        console.print("Please run the Structure Preparation & Relaxation module first.")
        sys.exit(1)
    
    # Check for minimized restart file
    min_rst7 = spr_dir / "min.rst7"
    
    # Determine topology file (prioritize *new.prmtop over *reord.prmtop)
    topology_files = list(spr_dir.glob("*new.prmtop")) or list(spr_dir.glob("*reord.prmtop"))
    
    if not topology_files:
        console.print("[bold red]Error:[/] No topology files found.")
        console.print("Please run the Structure Preparation & Relaxation module first.")
        sys.exit(1)
    
    topology = topology_files[0]
    
    # If min.rst7 exists, ask user if they want to use it
    if min_rst7.exists():
        use_existing = interaction_manager.yes_no_prompt(
            "use_existing_min",
            "\nA minimized structure (min.rst7) already exists. Would you like to use it?"
        )
        
        if not use_existing:
            # Trigger manual selection if user doesn't want to use existing
            manual_pdb = interaction_manager.prompt(
                "manual_pdb_path",
                "Enter the full path to the PDB file you want to analyze:",
                input_type=str
            )
            
            pdb_path = Path(manual_pdb)
            if not pdb_path.exists():
                console.print("[bold red]Error:[/] Specified PDB file does not exist.")
                sys.exit(1)
            
            return str(pdb_path)
    
    # Run minimization if min.rst7 doesn't exist
    if not min_rst7.exists():
        console.print("[yellow]Minimization restart file not found. Running minimization...[/]")
        
        # Determine restart file
        restart_files = list(spr_dir.glob("*new.rst7")) or list(spr_dir.glob("*reord.rst7"))
        
        if not restart_files:
            console.print("[bold red]Error:[/] No restart files found.")
            console.print("Please run the Structure Preparation & Relaxation module first.")
            sys.exit(1)
        
        restart = restart_files[0]
        
        # Run minimization
        from biodc.core.preparation.struct_relax import StructureRelaxation
        relaxation = StructureRelaxation(
            launch_dir=launch_dir,
            structure_dir=spr_dir,
            out_prefix=topology.stem.replace('_new', '').replace('_reord', ''),
            solvent_env='implicit',  # default, can be made configurable if needed
            input_dict={}  # empty dict for now
        )
        relaxation.minimize()
    
    # Create EE directory if it doesn't exist
    ee_dir = launch_dir / "EE"
    ee_dir.mkdir(exist_ok=True)
    
    # Convert restart to PDB
    output_pdb = ee_dir / "min.pdb"
    
    # Run ambpdb conversion
    print("Converting selected structure to min.pdb in the Energetic Evalulation (EE) directory.")
    conversion_cmd = (
        f"ambpdb -p {topology} -c {min_rst7} > {output_pdb}"
    )
    
    try:
        subprocess.run(conversion_cmd, shell=True, check=True)
    except subprocess.CalledProcessError:
        console.print("[bold red]Error:[/] Failed to convert restart file to PDB.")
        sys.exit(1)
    
    return str(output_pdb)
