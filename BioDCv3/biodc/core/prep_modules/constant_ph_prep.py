# biodc/core/preparation/constant_ph_prep.py

"""
Constant pH Dynamics Preparation Module for BioDC

This module handles the preparation of structures for constant pH molecular dynamics,
including residue reordering and cpin file generation.
"""

import os
import sys
import subprocess
from pathlib import Path
from typing import Dict, Tuple, List, Optional

from biodc.core.prep_modules import select_ph_active_sites
from biodc.utils.interaction import InteractionManager

class ConstantpHPreparation:
    """
    Handles preparation for constant pH molecular dynamics simulations.
    
    Provides methods for:
    - Reordering residues 
    - Generating cpin files
    - Interactive residue selection
    - Recording input choices
    """
    
    def __init__(self, launch_dir: Path, out_prefix: str):
        """
        Initialize constant pH preparation.
        
        Args:
            launch_dir: Directory where the preparation is taking place
            out_prefix: Output file prefix
        """
        self.launch_dir = Path(launch_dir)
        self.out_prefix = out_prefix
        self.interaction_manager = InteractionManager(
            launch_dir=self.launch_dir
        )
    
    def reorder_residues(self) -> Tuple[str, str, str]:
        """
        Reorder residues using CPPTRAJ to match constant pH dynamics requirements.
        
        Returns:
            Tuple of (reordered prmtop path, original pdb path, reordered pdb path)
        """
        # Prepare CPPTRAJ input for reordering
        reorder_input = f"""
parm {self.out_prefix}.prmtop
trajin {self.out_prefix}.rst7
fixatomorder parmout {self.out_prefix}_reord.prmtop
trajout {self.out_prefix}_original_resnum.pdb 
trajout {self.out_prefix}_reordered_resnum.pdb topresnum
trajout {self.out_prefix}_reord.rst7 topresnum
run
quit
        """
        
        # Write CPPTRAJ input file
        with open("ReorderRes.in", "w") as f:
            f.write(reorder_input)
        
        # Run CPPTRAJ
        subprocess.run(
            "cpptraj -i ReorderRes.in > ReorderRes.log 2> /dev/null", 
            shell=True, 
            check=True
        )
        
        return (
            f"{self.out_prefix}_reord.prmtop",
            f"{self.out_prefix}_original_resnum.pdb",
            f"{self.out_prefix}_reordered_resnum.pdb"
        )
    
    def select_titratable_residues(
        self, 
        reordered_pdb: str, 
        input_dict: Optional[Dict] = None
    ) -> Tuple[str, str, str, str, str, str]:
        """
        Interactively select titratable residues for constant pH simulation.
        
        Args:
            reordered_pdb: Path to reordered PDB
            input_dict: Optional dictionary of pre-existing inputs
        
        Returns:
            Tuple of selected residue ID strings for ASP, GLU, HIS, LYS, TYR, PRN
        """
        # Use existing select_ph_active_sites module
        return select_ph_active_sites.select_ph_active_sites(
            reordered_pdb, 
            "SecondPass", 
            input_dict or {}, 
            self.launch_dir
        )
    
    def generate_cpin_file(
        self, 
        reordered_prmtop: str, 
        selected_residues: Dict[str, str]
    ) -> str:
        """
        Generate cpin file for constant pH dynamics.
        
        Args:
            reordered_prmtop: Path to reordered topology file
            selected_residues: Dictionary of selected residue IDs
        
        Returns:
            Path to generated cpin file
        """
        # Prepare residue names
        res_names = []
        if selected_residues.get('ASP'):
            res_names.append('AS4')
        if selected_residues.get('GLU'):
            res_names.append('GL4')
        if selected_residues.get('HIS'):
            res_names.append('HIP')
        if selected_residues.get('LYS'):
            res_names.append('LYS')
        if selected_residues.get('TYR'):
            res_names.append('TYR')
        if selected_residues.get('PRN'):
            res_names.append('PRN')
        
        # Prepare residue IDs
        res_ids = []
        for key in ['ASP', 'GLU', 'HIS', 'LYS', 'TYR', 'PRN']:
            if selected_residues.get(key):
                res_ids.append(selected_residues[key])
        
        # If no residues selected, return None or raise an appropriate error
        if not res_names or not res_ids:
            print("No titratable residues selected for constant pH dynamics.")
            return None
        
        # Construct cpinutil command
        cmd = [
            "cpinutil.py",
            "-resnames", " ".join(res_names),
            "-resnums", " ".join(res_ids),
            "-p", reordered_prmtop,
            "-igb", "2",
            "-op", f"{self.out_prefix}_new.prmtop",
            "-o", f"{self.out_prefix}.cpin"
        ]
        
        # Run cpinutil
        subprocess.run(" ".join(cmd), shell=True, check=True)
        
        return f"{self.out_prefix}.cpin"
    
    def prepare_constant_ph(
        self, 
        input_dict: Optional[Dict] = None
    ) -> Dict[str, str]:
        """
        Full constant pH dynamics preparation workflow.
        
        Args:
            input_dict: Optional dictionary of pre-existing inputs
        
        Returns:
            Dictionary of generated file paths and selected residue IDs
        """
        print("""
 Now Generating the cpin file because you indicated that you want to run molecular dynamics
 with titratable residues.""")
    
        # Reorder residues
        reordered_prmtop, original_pdb, reordered_pdb = self.reorder_residues()
    
        # Check for existing titratable residue selections
        existing_selections = {
            'ASP': input_dict.get('SelASPIDs_1', ''),
            'GLU': input_dict.get('SelGLUIDs_1', ''),
            'HIS': input_dict.get('SelHISIDs_1', ''),
            'LYS': input_dict.get('SelLYSIDs_1', ''),
            'TYR': input_dict.get('SelTYRIDs_1', ''),
            'PRN': input_dict.get('SelPRNIDs_1', '')
        }
    
        # If no previous selections at all, run interactive selection
        if not any(existing_selections.values()):
            print("\n No previous titratable residue selections found. Selecting interactively.")
            selected_residues = dict(zip(
                ['ASP', 'GLU', 'HIS', 'LYS', 'TYR', 'PRN'], 
                self.select_titratable_residues(reordered_pdb, input_dict)
            ))
        else:
            # Use existing selections where available
            print("\n Using previously selected titratable residues.")
            selected_residues = {}
            
            # For each residue type, use existing selection or run interactive selection
            for res_type in ['ASP', 'GLU', 'HIS', 'LYS', 'TYR', 'PRN']:
                if existing_selections[res_type]:
                    selected_residues[res_type] = existing_selections[res_type]
                else:
                    # If no selection for this residue type, run interactive selection
                    print(f"\n No previous selection found for {res_type}.")
                    temp_selections = self.select_titratable_residues(reordered_pdb, input_dict)
                    selected_residues[res_type] = temp_selections[['ASP', 'GLU', 'HIS', 'LYS', 'TYR', 'PRN'].index(res_type)]
    
        # Generate cpin file
        cpin_file = self.generate_cpin_file(reordered_prmtop, selected_residues)
    
        # Prepare return dictionary
        result = {
            'reordered_prmtop': reordered_prmtop,
            'original_pdb': original_pdb,
            'reordered_pdb': reordered_pdb,
            'cpin_file': cpin_file
        }
    
        # Add selected residue IDs to the result
        result.update(selected_residues)
    
        return result

# Example usage function
def prepare_constant_ph_dynamics(
    out_prefix: str, 
    launch_dir: Path, 
    input_dict: Optional[Dict] = None
) -> Dict[str, str]:
    """
    Convenience function to prepare constant pH dynamics.
    
    Args:
        out_prefix: Output file prefix
        launch_dir: Working directory
        input_dict: Optional dictionary of pre-existing inputs
    
    Returns:
        Dictionary of generated file paths and selected residue IDs
    """
    prep = ConstantpHPreparation(launch_dir, out_prefix)
    return prep.prepare_constant_ph(input_dict)
