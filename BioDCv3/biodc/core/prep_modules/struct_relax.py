# biodc/core/preparation/struct_relax.py

"""
Structure Relaxation Module for BioDC

Handles energy minimization for molecular structures in explicit or implicit solvent.
"""

import os
import sys
import subprocess
from pathlib import Path
from typing import Dict, Optional, Literal

from biodc.utils.interaction import InteractionManager

class StructureRelaxation:
    """
    Manages structure relaxation (energy minimization) for molecular systems.
    
    Supports both explicit and implicit solvent minimization using SANDER or PMEMD.
    """
    
    def __init__(self, 
                 launch_dir: Path, 
                 structure_dir: Path, 
                 out_prefix: str, 
                 solvent_env: str, 
                 input_dict: Optional[Dict] = None):
        """
        Initialize structure relaxation parameters.
        
        Args:
            launch_dir: Directory where the script is launched
            structure_dir: Directory containing structure files
            out_prefix: Output file prefix
            solvent_env: Solvent environment type ('explicit' or 'implicit')
            input_dict: Optional dictionary of pre-existing inputs
        """
        self.launch_dir = Path(launch_dir)
        self.structure_dir = Path(structure_dir)
        self.out_prefix = out_prefix
        self.solvent_env = solvent_env.lower()
        self.input_dict = input_dict or {}
        
        # Create InteractionManager
        self.interaction_manager = InteractionManager(
            launch_dir=self.launch_dir, 
            input_dict=input_dict
        )
    
    def _get_topology_and_restart_files(self) -> Dict[str, Path]:
        """
        Determine available topology and restart files.
        
        Returns:
            Dictionary of available topology and restart files
        """
        possible_files = [
            (self.structure_dir / f"{self.out_prefix}_new.prmtop", 
             self.structure_dir / f"{self.out_prefix}_reord.rst7"),
            (self.structure_dir / f"{self.out_prefix}_reord.prmtop", 
             self.structure_dir / f"{self.out_prefix}_reord.rst7"),
            (self.structure_dir / f"{self.out_prefix}.prmtop", 
             self.structure_dir / f"{self.out_prefix}.rst7")
        ]
        
        for topology, restart in possible_files:
            if topology.exists() and restart.exists():
                print(f"Found {topology} and {restart}")
                return {
                    "topology": topology,
                    "restart": restart,
                    "reference": restart
                }
        
        raise FileNotFoundError(
            "No valid topology and restart files found. "
            "Something went wrong in the preceding steps!"
        )
    
    def _generate_minimization_input(self) -> str:
        """
        Generate AMBER minimization input file based on solvent environment.
        
        Returns:
            Minimization input file contents
        """
        if self.solvent_env in ["explicit", "exp"]:
            return """
&cntrl
  imin=1,            ! Perform an energy minimization
  ntb=1,             ! Constant volume
  cut=10.0,          ! Non-bonded cutoff in angstroms
  ntmin=1,           ! Steepest descent + conjugate gradient method
  ncyc=1000,         ! Number of steepest descent cycles
  maxcyc=5000,       ! Maximum number of minimization cycles
  ntwr=100,          ! Restart file written every ntwr steps
  ntwx=100,          ! Trajectory file written every ntwx steps
  ntpr=100,          ! The mdout and mdinfo files written every ntpr steps
  ntr=1,             ! Turn on positional restraints
  restraintmask='@CA,C,O,N&!:WAT|@FE,NA,NB,NC,ND,C3D,C2A,C3B,C2C,CA,CB',
  restraint_wt=10.0, ! 10 kcal/mol.A**2 restraint force constant
/
            """
        elif self.solvent_env in ["implicit", "imp"]:
            return """
&cntrl
  imin=1,            ! Perform an energy minimization
  ntb=0,             ! Non-periodic
  cut=9999,          ! Non-bonded cutoff in Å
  ntmin=1,           ! Steepest descent + conjugate gradient method 
  ncyc=200,          ! Number of steepest descent cycles
  maxcyc=500,        ! Maximum number of minimization cycles
  igb=2,             ! Generalized Born implicit solvent model
  saltcon=0.1,       ! salt concentration in M
  ntwr=100,          ! Restart file written every ntwr steps
  ntwx=100,          ! Trajectory file written every ntwx steps
  ntpr=100,          ! The mdout and mdinfo files written every ntpr steps
  ntr=1,             ! Turn on positional restraints
  restraintmask='@CA,C,O,N&!:WAT|@FE,NA,NB,NC,ND,C3D,C2A,C3B,C2C,CA,CB',
  restraint_wt=10.0, ! 10 kcal/mol.A**2 restraint force constant
/
            """
        else:
            raise ValueError(f"Invalid solvent environment: {self.solvent_env}")
    
    def _get_minimization_method(self) -> Literal["sander", "pmemd"]:
        """
        Interactively or automatically select minimization method.
        
        Returns:
            Minimization method ('sander' or 'pmemd')
        """
        method = self.interaction_manager.prompt(
            "StructRelaxCompChoice", 
            "\nRun the minimization using SANDER (S) or PMEMD (P)?", 
            choices=['S', 'P']
        )
        
        # Normalize method
        method = method.lower()
        if method in ["sander", "s"]:
            return "sander"
        elif method in ["pmemd", "p"]:
            return "pmemd"
        else:
            raise ValueError("Invalid minimization method. Choose Sander or PMEMD.")
    
    def _get_parallel_processors(self) -> Optional[int]:
        """
        Get number of processors for parallel minimization.
        
        Returns:
            Number of processors or None if using Sander
        """
        if self.minimization_method == "pmemd":
            # Prompt for number of processors
            nproc = self.interaction_manager.prompt(
                "NProc", 
                "Parallelize the minimization over how many CPUs?", 
                input_type=int
            )
            
            return nproc
        
        return None
    
    def minimize(self):
        """
        Perform energy minimization on the structure.
        """
        # Ensure we're in the correct directory
        original_cwd = Path.cwd()
        os.chdir(self.structure_dir)
        
        try:
            # Get topology and restart files
            files = self._get_topology_and_restart_files()
            
            # Write minimization input
            with open('min.in', 'w') as f:
                f.write(self._generate_minimization_input())
            
            # Determine minimization method
            self.minimization_method = self._get_minimization_method()
            
            # Prepare command
            if self.minimization_method == "sander":
                cmd = (f"sander -O -i min.in -o min.out "
                       f"-p {files['topology']} "
                       f"-c {files['restart']} "
                       f"-inf min.mdinfo -r min.rst7 "
                       f"-ref {files['reference']}")
            else:  # pmemd
                nproc = self._get_parallel_processors()
                cmd = (f"mpirun -np {nproc} pmemd.MPI -O -i min.in -o min.out "
                       f"-p {files['topology']} "
                       f"-c {files['restart']} "
                       f"-inf min.mdinfo -r min.rst7 "
                       f"-ref {files['reference']}")
            
            # Run minimization
            print("Running minimization ...")
            subprocess.run(cmd, shell=True, check=True)
            print("Minimization finished!")
        
        except Exception as e:
            print(f"Error during minimization: {e}")
            raise
        finally:
            # Return to original working directory
            os.chdir(original_cwd)

def struct_relax(launch_dir: Path, 
                 struc_dir: Path, 
                 out_prefix: str, 
                 solv_env: str, 
                 input_dict: Dict):
    """
    Convenience function to run structure relaxation.
    
    Args:
        launch_dir: Directory where the script is launched
        struc_dir: Directory containing structure files
        out_prefix: Output file prefix
        solv_env: Solvent environment type
        input_dict: Dictionary of input parameters
    """
    relaxation = StructureRelaxation(
        launch_dir, 
        struc_dir, 
        out_prefix, 
        solv_env, 
        input_dict
    )
    relaxation.minimize()
