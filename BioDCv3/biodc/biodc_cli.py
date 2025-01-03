"""
BioDC - A program that automates and accelerates the computation 
of redox currents in (polymeric) multi-heme cytochromes.

Written by Matthew J. Guberman-Pfeffer and Caleb L. Herron
"""
import os
import sys
import shutil
from datetime import datetime
from pathlib import Path
from typing import Dict, Optional
import click
from rich.console import Console
from rich.logging import RichHandler
import logging
from dataclasses import dataclass

# Ensure project root is in Python path
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, project_root)

# Explicitly import and assign to global
from biodc.utils.interaction import InteractionManager
globals()['InteractionManager'] = InteractionManager

from biodc.core import preparation 
from biodc.core.energetics import EnergeticEvaluation, EnergeticParameters
from biodc.core.kinetics import KineticEvaluation, KineticParameters
from biodc.utils import file_handling
from biodc.config import Config
from biodc.core.structure_utils import get_prepared_structure

import traceback

# Setup logging
logging.basicConfig(
    level="INFO",
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(rich_tracebacks=True)]
)
logger = logging.getLogger("biodc")
console = Console()

@dataclass
class WorkflowParameters:
    """Container for storing workflow-level parameters."""
    out_prefix: Optional[str] = None
    solv_env: Optional[str] = None
    prepared_pdb: Optional[str] = None
    energetic_params: Optional[EnergeticParameters] = None
    kinetic_params: Optional[KineticParameters] = None

class BioDC:
    """
    Main workflow orchestrator for BioDC, providing a modular and interactive approach 
    to computing redox currents in multi-heme cytochromes.
    """
    def __init__(self, config: Dict = None):
        config = config or {}
        self.launch_dir = Path.cwd()
        self.prog_dir = Path(__file__).parent
        self.forcefield_dir = self.prog_dir / "data" / "forcefield"
        self.config = config
        self._validate_installation()

        self.session = InteractionManager(
            launch_dir=self.launch_dir,
            input_dict=self.config
        )
        
        self.workflow_params = WorkflowParameters()

    def _validate_installation(self):
        if not (self.prog_dir / "data" / "forcefield").exists():
            raise FileNotFoundError("""
The forcefield directory is missing from the BioDC package.
Please ensure the forcefield files are located in:
biodc/data/forcefield/
""")

    def handle_directory(self, dir_name: str) -> Path:
        work_dir = self.launch_dir / dir_name
       
        if work_dir.exists():
            dir_fate = self.session.yes_no_prompt(
                "DirFate",
                f"\n   A directory named {dir_name} already exists.\n"
                f"   It may have been created from a prior\n"
                f"   session with BioDC. Would you like to\n"
                f"   delete it and start over?"
            )
           
            if dir_fate:
                shutil.rmtree(work_dir)
                console.print("\n   Old directory deleted.")
                work_dir.mkdir()
                console.print("   New directory Created.")
           
        else:
            console.print(f"\n Creating a directory called {dir_name} for this module.")
            work_dir.mkdir()
           
        console.print(f" Now, the current working directory is: {work_dir}")
        return work_dir

    def prepare_structure(self) -> WorkflowParameters:
        console.print("""
===================================================================  
First: Structure Preparation & Relaxation
===================================================================
""")
        work_dir = self.handle_directory("SPR")
        os.chdir(work_dir)
        
        out_prefix, solv_env = preparation.run(
            self.launch_dir,
            self.forcefield_dir,
            work_dir,
            self.session.get_input_dict()
        )
        
        self.workflow_params.out_prefix = out_prefix
        self.workflow_params.solv_env = solv_env
#       self.workflow_params.prepared_pdb = get_prepared_structure(
#           self.launch_dir, 
#           self.session
#       )
        
        return self.workflow_params

    def evaluate_energetics(self) -> WorkflowParameters:
        console.print("""
===================================================================  
Second: Energetic Evaluation
===================================================================
""")
        
        work_dir = self.handle_directory("EE")
        os.chdir(work_dir)
    
        if not self.workflow_params.prepared_pdb:
            self.workflow_params.prepared_pdb = get_prepared_structure(
                self.launch_dir, 
                self.session
            )
        
        energetic_eval = EnergeticEvaluation(
            interaction_manager=self.session,
            launch_dir=self.launch_dir,
            forcefield_dir=self.forcefield_dir,
            pdb_file=self.workflow_params.prepared_pdb
        )
        
        self.workflow_params.energetic_params = energetic_eval.run(
            self.workflow_params.prepared_pdb
        )
        
        return self.workflow_params

    def evaluate_kinetics(self) -> WorkflowParameters:
        console.print("""
===================================================================  
Third: Kinetic Evaluation
===================================================================
""")
       
        # Verify EE directory exists
        if not (self.launch_dir / "EE").exists():
            raise FileNotFoundError("EE directory not found. Please run energetic evaluation first.")
        
        work_dir = self.handle_directory("KE")
        os.chdir(work_dir)
        
        kinetic_eval = KineticEvaluation(
            interaction_manager=self.session,
            launch_dir=self.launch_dir
        )
        
        self.workflow_params.kinetic_params = kinetic_eval.run()
        
        return self.workflow_params

    def run_workflow(self) -> WorkflowParameters:
        print(f"""
================================================================== 
                       Welcome to BioDC
            A program that automates and accelerates
              the computation of redox currents in
               (polymeric) multi-heme cytochromes 

   Written by Matthew J. Guberman-Pfeffer and Caleb L. Herron
                 Last Updated: 7/15/2024

Start time: {datetime.now()}
Directory paths:
  Program:                    {self.launch_dir}
  Current working directory:  {os.getcwd()}
================================================================== 

BioDC presents a highly modular workflow that has three 
major divisions: 
   (1) Structure Preparation & Relaxation
   (2) Energetic Evaluation
   (3) Kinetic Evaluation
""")
       
        division = self.session.prompt(
            "DivSel",
            "\n Which of these divisions would you like to perform?\n"
            " (Enter zero \"0\" to be guided through the entire\n"
            " workflow.) (0/1/2/3)",
            choices=['0', '1', '2', '3']
        )
       
        try:
            if division in ('0', '1'):
                self.prepare_structure()
           
            if division in ('0', '2'):
                if division == '0':
                    os.chdir(self.launch_dir)
                self.evaluate_energetics()
           
            if division in ('0', '3'):
                if division == '0':
                    os.chdir(self.launch_dir)
                self.evaluate_kinetics()
           
            if division not in ('0', '1', '2', '3'):
                sys.exit("\n Please re-launch BioDC and select one of the available modules.")
        
        except Exception as e:
            print(f"Detailed error: {e}")
            traceback.print_exc()
            sys.exit(1)
        
        return self.workflow_params

@click.command()
@click.option('--config', '-c', type=click.Path(exists=True),
              help='Path to configuration file')
def run(config):
    """Run BioDC workflow."""
    try:
        initial_config = {}
        if config:
            initial_config = file_handling.read_input(config)
        
        biodc = BioDC(initial_config)
        biodc.run_workflow()
        
    except Exception as e:
        print(f"Detailed error: {e}")
        traceback.print_exc()
        sys.exit(1)

def main():
    """Entry point for the CLI application."""
    try:
        run(prog_name="biodc")  
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    run()
