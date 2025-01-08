# biodc/core/preparation/create_res_indexing.py

"""
Module for creating heme residue indexing files.
Handles all possible His-[distal] ligand combinations and maintains
information needed for subsequent structure preparation steps.
"""

import os
import sys
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass
import logging

from biodc.core.prep_modules.ligand_detection import LigandDetector, MetalCoordination
from biodc.utils.interaction import InteractionManager

logger = logging.getLogger(__name__)

# All possible distal ligand types and their codes
DISTAL_LIGANDS = {
    'HIS': 'H',  # His-His
    'MET': 'M',  # His-Met
    'CYS': 'C',  # His-Cys
    'TYR': 'Y',  # His-Tyr
    'ASP': 'D',  # His-Asp
    'GLU': 'E',  # His-Glu
    'ASN': 'N',  # His-Asn
    'GLN': 'Q',  # His-Gln
    'LYS': 'K'   # His-Lys
}

@dataclass
class HemeEnvironment:
    """Container for heme environment information."""
    heme_id: int            # Original heme residue ID
    cys_b: Optional[int]   # B-ring Cys for c-type (None for b-type)
    cys_c: Optional[int]   # C-ring Cys for c-type (None for b-type)
    his_p: int            # Proximal His
    distal_ligand: int    # Distal ligand residue ID
    distal_type: str      # Type of distal ligand (e.g., 'HIS', 'MET', etc.)
    is_c_type: bool       # True for c-type, False for b-type

    @property
    def ligation_code(self) -> str:
        """Get two-letter code for ligation type (e.g., 'HH', 'HM', 'HD', etc.)."""
        return f"H{DISTAL_LIGANDS[self.distal_type]}"

class ResidueIndexing:
    """Handles creation of residue indexing files for heme setup."""
    
    def __init__(self, pdb_path: Path, launch_dir: Path):
        self.pdb_path = pdb_path
        self.launch_dir = launch_dir
        self.detector = LigandDetector(str(pdb_path))
        self.environments: List[HemeEnvironment] = []

    def analyze_heme_environments(self) -> List[HemeEnvironment]:
        """Analyze all heme environments in the structure."""
        print("\n")
        print("-" * 60)
        print("Analyzing heme environments...")
        print("-" * 60)

        # First find thioether bonds to identify c-type hemes
        print("Looking for thioether bonds...")
        thioether_bonds = self.detector.find_thioether_bonds()

        # Then find metal coordinations
        print("Looking for metal coordinations...")
        coordinations = self.detector.find_metal_coordinations('FE', ['HEM', 'HEME', 'HEC'])

        # Group coordinations by heme
        heme_coords = {}
        for coord in coordinations:
            if coord.metal_id not in heme_coords:
                heme_coords[coord.metal_id] = []
            heme_coords[coord.metal_id].append(coord)

        print(f"\nFound {len(heme_coords)} hemes")

        # Process each heme
        for heme_id, coords in sorted(heme_coords.items()):
            print(f"\nAnalyzing Heme {heme_id}:")
            print(f"  Found {len(coords)} coordinating residues:")
#           for c in coords:
#               print(f"    {c.residue_name} {c.residue_id}")

            # Find proximal His and other ligands
            his_ligands = [c for c in coords if c.residue_name == 'HIS']
            print(f"  His ligands: {[c.residue_id for c in his_ligands]}")

            other_ligands = {c.residue_name: c for c in coords
                        if c.residue_name in DISTAL_LIGANDS
                        and c.residue_name != 'HIS'}
            print(f"  Other ligands: {[(k, v.residue_id) for k,v in other_ligands.items()]}")

            # Must have at least one His
            if not his_ligands:
                print(f"Error: Heme {heme_id} missing proximal histidine")
                raise ValueError(f"Heme {heme_id} missing proximal histidine")

            # Check thioether bonds to determine if c-type
            cys_ids = thioether_bonds.get(heme_id, [])
            is_c_type = len(cys_ids) == 2
            print(f"  Thioether-bonded Cys: {cys_ids}")
            print(f"  Heme type: {'c-type' if is_c_type else 'b-type'}")

            # Determine distal ligand
            if len(his_ligands) == 2:
                # His-His case
                distal_type = 'HIS'
                distal_ligand = his_ligands[1].residue_id
            else:
                # Find first non-His ligand
                for ligand_type in DISTAL_LIGANDS:
                    if ligand_type in other_ligands:
                        distal_type = ligand_type
                        distal_ligand = other_ligands[ligand_type].residue_id
                        break
                else:
                    raise ValueError(f"Heme {heme_id} missing distal ligand")

            # Create environment
            env = HemeEnvironment(
                heme_id=heme_id,
                cys_b=cys_ids[0] if is_c_type else None,
                cys_c=cys_ids[1] if is_c_type else None,
                his_p=his_ligands[0].residue_id,
                distal_ligand=distal_ligand,
                distal_type=distal_type,
                is_c_type=is_c_type
            )

            self.environments.append(env)
            print(f"  Created environment: c-type={env.is_c_type}, HisP={env.his_p}, "
                f"Distal={env.distal_ligand} ({env.distal_type})")

        return self.environments

    def validate_environments(self) -> bool:
        """Validate identified heme environments."""
        print("\n")
        print("-" * 60)
        print("Validating heme environments...")
        print("-" * 60)

        if not self.environments:
            print("Error: No heme environments were found!")
            return False

        # Check each environment in detail
        for env in self.environments:
            print(f"\nValidating Heme {env.heme_id}:")
            print(f"  Proximal His: {env.his_p}")
            print(f"  Distal ligand ({env.distal_type}): {env.distal_ligand}")

            if not env.his_p:
                print("  Error: Missing proximal His!")
                return False

            if not env.distal_ligand:
                print("  Error: Missing distal ligand!")
                return False

            if env.is_c_type:
                print(f"  C-type heme, Cys B/C: {env.cys_b}/{env.cys_c}")
                if not env.cys_b or not env.cys_c:
                    print("  Error: C-type heme missing required Cys residues!")
                    return False

            print("  Environment valid!")

        print("\nAll heme environments validated successfully!")
        return True

    def print_environment_summary(self):
        """Print summary of identified heme environments."""
        print("\n")
        print("-" * 60)
        print("Identified heme environments:")
        print("-" * 60)

        for env in self.environments:
            heme_type = "c-type" if env.is_c_type else "b-type"
            print(f"\nHeme {env.heme_id} ({heme_type}):")
            if env.is_c_type:
                print(f"  Cys B/C: {env.cys_b}/{env.cys_c}")
            print(f"  His (proximal): {env.his_p}")
            print(f"  {env.distal_type} (distal): {env.distal_ligand}")

    def write_indexing_files(self):
        """Write all required indexing files to the current working directory."""
        # Write main indexing file
        with open(Path.cwd() / "ResIndexing.txt", 'w') as f:
            for env in self.environments:
                if env.is_c_type:
                    f.write(f"{env.cys_b} {env.cys_c} {env.his_p} {env.distal_ligand} "
                           f"{env.heme_id} c {env.ligation_code}\n")
                else:
                    f.write(f"{env.his_p} {env.distal_ligand} "
                           f"{env.heme_id} b {env.ligation_code}\n")

        # Write specific ligation files
        for distal_type in DISTAL_LIGANDS:
            code = DISTAL_LIGANDS[distal_type]
            filename = f"ResIndexing_His{distal_type}Heme.txt"

            with open(Path.cwd() / filename, 'w') as f:
                for env in self.environments:
                    if env.distal_type == distal_type:
                        if env.is_c_type:
                            f.write(f"{env.cys_b} {env.cys_c} {env.his_p} {env.distal_ligand} "
                                    f"{env.heme_id} c {env.ligation_code}\n")
                        else:
                            f.write(f"{env.his_p} {env.distal_ligand} {env.heme_id} "
                                    f"b {env.ligation_code}\n")

def create_res_indexing(pdb: str, input_dict: Dict, launch_dir: Path) -> str:
    """Main function for creating residue indexing files."""
    interaction_manager = InteractionManager(
        launch_dir=launch_dir, 
        input_dict=input_dict
    )
    
    launch_dir = Path(launch_dir)
    pdb_path = Path(f"{pdb}.pdb")
    pdb_basename = pdb_path.stem

    # Check for corrected indexing
    corrected_path = Path.cwd() / "CorrectedResIndexing.txt"
    if corrected_path.exists():
        print("\n Found CorrectedResIndexing.txt! \n Copying CorrectedResIndexing.txt to ResIndexing.txt.")
        shutil.copy2(corrected_path, Path.cwd() / "ResIndexing.txt")
        return pdb_basename

    if (launch_dir / "CorrectedResIndexing.txt").exists():
        print("\n Found CorrectedResIndexing.txt in launch directory! \n Copying CorrectedResIndexing.txt to ResIndexing.txt.")
        shutil.copy2(launch_dir / "CorrectedResIndexing.txt", Path.cwd() / "ResIndexing.txt")
        return pdb_basename
    
    print("=" * 60)
    method = interaction_manager.prompt(
        "CreateResIndexingMethod", 
        "\n We need to create a file (ResIndexing.txt) that identifies\n"
        " the IDs of the Cys, His or Met residues bound to the heme group.\n"    
        " Would you like to create it automatically or manually?",
        choices=['auto', 'man']
    )
    
    if method.lower() in ('auto', 'a'):
        indexing = ResidueIndexing(pdb_path, launch_dir)
        indexing.analyze_heme_environments()
        
        if indexing.validate_environments():
            indexing.write_indexing_files()
            indexing.print_environment_summary()
            print("\n" * 60)
            print("*" * 60)
            print(""" Please verify the identified residues in the ResIndexing.txt file.
 If corrections are needed, save them in CorrectedResIndexing.txt.
 
 Note: The program identifies ligands based on distance criteria.
 Please check that the assigned ligands are correct and make any
 necessary adjustments in CorrectedResIndexing.txt.""")
            print("*" * 60)
        else:
            raise ValueError("Failed to validate heme environments")
            
    elif method.lower() in ('manual', 'm'):
        print("""
 To create ResIndexing.txt by hand:
    > Create a txt file with an editor of your choosing (e.g. 
      vi ResIndexing.txt). 
        
    > For c-type hemes, each line needs 7 space-separated fields:
      CysB CysC HisP Distal Heme c HX
        
    > For b-type hemes, each line needs 5 space-separated fields:
      HisP Distal Heme b HX
      
    Where X is the code for the distal ligand:
      H = His, M = Met, C = Cys, Y = Tyr, D = Asp,
      E = Glu, N = Asn, Q = Gln, K = Lys
        
    See documentation for detailed format explanation.""")
        
        if not (launch_dir / "ResIndexing.txt").exists():
            raise FileNotFoundError("Please create ResIndexing.txt and run the script again.")
            
    else:
        raise ValueError("Invalid method selection")
        
    return pdb_basename

if __name__ == "__main__":
    # For testing
    pdb = "example"
    input_dict = {}
    launch_dir = Path.cwd()
    create_res_indexing(pdb, input_dict, launch_dir)
