# biodc/core/preparation/select_disulfides.py

"""
Module for selecting and defining disulfide bonds in protein structures.
Includes automatic detection of potential disulfide bonds.
"""

from pathlib import Path
from typing import Dict, Tuple, List, Optional
from Bio import PDB
from Bio.PDB.Atom import Atom
from Bio.PDB.Residue import Residue
import numpy as np

from biodc.utils.interaction import InteractionManager

# Constants
MAX_DISULFIDE_DISTANCE = 3.0  # Angstroms - slightly larger than S-S bond length to account for non-bonded Cys
OPTIMAL_SS_BOND_LENGTH = 2.05  # Optimal S-S bond length

class DisulfideFinder:
    """Handles detection and analysis of potential disulfide bonds."""
    
    def __init__(self, pdb_path: str):
        self.parser = PDB.PDBParser(QUIET=True)
        self.structure = self.parser.get_structure('protein', pdb_path)
        
    def find_potential_disulfides(self) -> List[Tuple[int, int, float]]:
        """
        Find all pairs of CYS residues that could potentially form disulfide bonds.
        Returns list of tuples: (res1_id, res2_id, distance)
        """
        potential_pairs = []
        cys_residues = []
        
        # Collect all CYS residues
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    if residue.get_resname() == "CYS" or residue.get_resname() == "CYX":
                        cys_residues.append(residue)
        
        # Check distances between all pairs
        for i, res1 in enumerate(cys_residues):
            for res2 in cys_residues[i+1:]:
                try:
                    # Get sulfur atoms
                    s1 = res1["SG"]
                    s2 = res2["SG"]
                    
                    # Calculate distance
                    distance = s1 - s2
                    
                    if distance <= MAX_DISULFIDE_DISTANCE:
                        # Get residue numbers from the PDB
                        res1_num = res1.get_id()[1]
                        res2_num = res2.get_id()[1]
                        potential_pairs.append((res1_num, res2_num, distance))
                        
                except KeyError:
                    # Skip if SG atom not found
                    continue
        
        # Sort by distance
        return sorted(potential_pairs, key=lambda x: x[2])

def select_disulfides(pdb_name: str, input_dict: Dict, launch_dir: Path) -> Tuple[str, List[List[int]]]:
    """Main function for selecting disulfide bonds.

    Args:
        pdb_path: Full path to the PDB file
        input_dict: Dictionary of input parameters
        launch_dir: Directory where script is launched
    """
    interaction_manager = InteractionManager(
        launch_dir=launch_dir,
        input_dict=input_dict
    )

    disulf_res_list = " "
    disulf_pair_id = []

    # Analyze structure for potential disulfides
    print("\n Analyzing structure for potential disulfide bonds...")
    pdb_path = f"{pdb_name}.pdb"
    finder = DisulfideFinder(pdb_path)
    potential_pairs = finder.find_potential_disulfides()

    if potential_pairs:
        print("\n Found potential disulfide bonds:")
        for i, (res1, res2, distance) in enumerate(potential_pairs, 1):
            print(f"  {i}. CYS {res1} - CYS {res2} (distance: {distance:.2f} Å)")

        use_detected = interaction_manager.yes_no_prompt(
            "UseDetectedDisulfides",
            "Would you like to review these detected pairs?"
        )

        if use_detected:
            for res1, res2, distance in potential_pairs:
                pair_str = f"{res1} {res2}"
                accept_pair = interaction_manager.yes_no_prompt(
                    f"AcceptDisulfidePair_{res1}_{res2}",
                    f"Accept disulfide pair CYS {res1} - CYS {res2} (distance: {distance:.2f} Å)?"
                )
                if accept_pair:
                    disulf_res_list += pair_str + " "
                    disulf_pair_id.append([res1, res2])

    # Allow manual entry of additional pairs
    add_manual = interaction_manager.yes_no_prompt(
        "AddManualPairs",
        "Would you like to manually add any additional disulfide pairs?"
    )

    if add_manual:
        while True:
            pair_str = interaction_manager.prompt(
                f"ManualDisulfidePair_{len(disulf_pair_id)+1}",
                f"Enter disulfide-linked Cys pair (format: res1 res2) or 'done' to finish"
            )
            if pair_str.lower() == 'done':
                break
            try:
                res1, res2 = map(int, pair_str.split())
                disulf_res_list += f"{res1} {res2} "
                disulf_pair_id.append([res1, res2])
            except ValueError:
                print("Invalid format. Please enter two numbers separated by a space.")

    if disulf_pair_id:
        with open(launch_dir / 'SPR' / 'DisulfideDefinitions.txt', 'w') as f:
            f.write(disulf_res_list)

    return disulf_res_list, disulf_pair_id
