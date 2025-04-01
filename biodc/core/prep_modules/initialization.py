# biodc/core/preparation/initialization.py

"""
Initialization module for BioDC structure preparation.
Verifies required programs and handles PDB selection.
"""

import os
import sys
import shutil
import warnings
import subprocess
from pathlib import Path
from typing import Dict
from Bio import PDB
from Bio.PDB import DSSP
from collections import Counter, defaultdict
from rich.console import Console
from rich.table import Table

from biodc.utils.interaction import InteractionManager

def verify_programs() -> bool:
    """
    Verify that required programs are available in the system PATH.
    
    Returns:
        Boolean indicating if all required programs are available
    """

    program_checks = {
        'vmd': 'vmd -dispdev text -e command.vmd',  # Simple exit command file
        'tleap': ['tleap', '-h'],
        'cpptraj': ['cpptraj', '-h'],
        'sander': ['sander', '-h'],
        'pmemd': ['pmemd', '-h']
    }

    # Create temporary VMD command file
    with open('command.vmd', 'w') as f:
        f.write('quit\n')

    missing_programs = []
    for program, cmd in program_checks.items():
        try:
            if program == 'vmd':
                result = subprocess.run(
                    cmd.split(),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    timeout=10  # Longer timeout for VMD
                )
#               print(f"VMD output: {result.stdout}\nVMD stderr: {result.stderr}")
                if "VMD for" in (result.stdout + result.stderr):
                    continue
            else:
                subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=5)
        except Exception as e:
            print(f"Error running {program}: {e}")
            missing_programs.append(program)

    os.unlink('command.vmd')  # Clean up

    if missing_programs:
        print(f"\nError: The following required programs are missing: {', '.join(missing_programs)}")
        return False

    return True

def analyze_composition(structure) -> None:
    """
    Analyze residue composition of structure.
    
    Args:
        structure: Biopython Structure object
    """
    from collections import Counter, defaultdict
    
    # Define residue categories
    AMINO_ACIDS = {
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE',
        'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
    }
    HEME_TYPES = {"HEM", "HEME", "HEA", "HEC"}
    WATER_TYPES = {"HOH", "WAT", "H2O"}
    ION_TYPES = {"NA", "CL", "K", "MG", "CA", "ZN", "FE", "MN"}

    # Initialize counters
    residue_categories = defaultdict(list)
    residue_counts = Counter()
    
    # First pass: categorize and count all residues
    print("\nAnalyzing structure composition...")
    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.resname
                orig_info = (chain.id, residue.id[1], resname)
                residue_counts[resname] += 1
                
                # Categorize the residue
                if resname in AMINO_ACIDS:
                    category = "amino_acids"
                elif resname in HEME_TYPES:
                    category = "hemes"
                elif resname in WATER_TYPES:
                    category = "waters"
                elif resname in ION_TYPES:
                    category = "ions"
                else:
                    category = "unknown"
                
                residue_categories[category].append((orig_info, residue))

    # Print composition summary
    total_residues = sum(residue_counts.values())
    print("\nStructure Composition Summary:")
    print("=" * 60)
    print(f"{'Residue Type':<20}{'Count':>10}{'Percentage':>15}{'Category':>15}")
    print("-" * 60)

    # Print standard amino acids
    if residue_categories["amino_acids"]:
        for aa in sorted(AMINO_ACIDS):
            if residue_counts[aa] > 0:
                percentage = (residue_counts[aa] / total_residues) * 100
                print(f"{aa:<20}{residue_counts[aa]:>10}{percentage:>14.1f}%{'Amino Acid':>15}")

    # Print hemes
    if residue_categories["hemes"]:
        print("-" * 60)
        for heme in sorted(set(res[0][2] for res in residue_categories["hemes"])):
            percentage = (residue_counts[heme] / total_residues) * 100
            print(f"{heme:<20}{residue_counts[heme]:>10}{percentage:>14.1f}%{'Heme':>15}")

    # Print unknown residues
    if residue_categories["unknown"]:
        print("-" * 60)
        print("Non-standard residues found:")
        for resname in sorted(set(res[0][2] for res in residue_categories["unknown"])):
            percentage = (residue_counts[resname] / total_residues) * 100
            print(f"{resname:<20}{residue_counts[resname]:>10}{percentage:>14.1f}%{'Unknown':>15}")

    # Print waters and ions (stripped)
    print("-" * 60)
    water_count = sum(residue_counts[w] for w in WATER_TYPES)
    ion_count = sum(residue_counts[i] for i in ION_TYPES)
    if water_count:
        print(f"{'Waters (stripped)':<20}{water_count:>10}{(water_count/total_residues)*100:>14.1f}%{'Water':>15}")
    if ion_count:
        print(f"{'Ions (stripped)':<20}{ion_count:>10}{(ion_count/total_residues)*100:>14.1f}%{'Ion':>15}")
    print("=" * 60)

def analyze_secondary_structure(structure) -> dict:
    """
    Comprehensive secondary structure analysis using DSSP.
    Analyzes helix and beta sheet characteristics in detail.
    """

    def find_continuous_segments(structure_list):
        """Helper to find continuous segments of same secondary structure."""
        segments = []
        current_segment = []
        
        for i, (res_id, ss) in enumerate(structure_list):
            if not current_segment:
                current_segment = [res_id]
            elif structure_list[i-1][0] + 1 == res_id:
                current_segment.append(res_id)
            else:
                segments.append(current_segment)
                current_segment = [res_id]
                
        if current_segment:
            segments.append(current_segment)
            
        return segments

    dssp_executable = '/usr/local/bin/mkdssp'
    
    # Verify executable exists and is accessible
    if not os.path.exists(dssp_executable):
        print(f"\nWarning: DSSP executable not found at {dssp_executable}")
        return None
    
    # Check executable permissions
    if not os.access(dssp_executable, os.X_OK):
        print(f"\nWarning: DSSP executable at {dssp_executable} is not executable")
        return None

    try:
        DSSP.DSSP_EXECUTABLE = dssp_executable
        model = structure[0]
        
        # Write a temporary PDB file for DSSP
        io = PDB.PDBIO()
        io.set_structure(structure)
        tmp_pdb = "temp_dssp.pdb"
        io.save(tmp_pdb)
        
        # Run DSSP on the temporary file
        dssp = DSSP(model, tmp_pdb, dssp='mkdssp')
        
        # Remove temporary file
        os.remove(tmp_pdb)
        
        # Unique SS codes from DSSP
        unique_ss_codes = {
            'H': 'Alpha Helix',
            'G': '3-10 Helix',
            'I': 'Pi Helix',
            'E': 'Beta Strand',
            'B': 'Beta Bridge',
            'T': 'Turn',
            'S': 'Bend',
            '-': 'Unstructured',
            'P': 'Pro/Pro-like'  
        }

        # Calculate total residues
        total_residues = len(dssp)

        # Print detailed analysis and summary table
        print("\nDetailed Secondary Structure Analysis:")
        print("=" * 80)

        print("\nSecondary Structure Composition:")
        print("-" * 40)
        print(f"{'2° Structure':<15} {'Count':>8} {'Avg Length':>12} {'Total Residues':>16} {'% of Total':>12}")
        print("-" * 80)

        # Comprehensive tracking of each SS type
        ss_type_tracking = {code: [] for code in unique_ss_codes}

        # Collect segment lengths for each SS type
        for res_id, residue in enumerate(dssp):
            ss = residue[2]
            ss_type_tracking[ss].append(res_id)

        # Print rows for each structure type
        for ss_code, ss_name in unique_ss_codes.items():
            data = ss_type_tracking[ss_code]
            
            # Group continuous segments
            segments = []
            if data:
                current_segment = [data[0]]
                for i in range(1, len(data)):
                    if data[i] == data[i-1] + 1:
                        current_segment.append(data[i])
                    else:
                        segments.append(current_segment)
                        current_segment = [data[i]]
                segments.append(current_segment)
            
            # Calculate statistics
            count = len(segments)
            total_res = len(data)
            avg_length = sum(len(seg) for seg in segments) / count if count > 0 else 0
            
            print(f"{ss_name:<15} {count:>8} {avg_length:>12.1f} {total_res:>16} {total_res/total_residues*100:>11.1f}%")

        # Return the tracking dictionary
        return ss_type_tracking
    
    except Exception as e:
        print(f"\nWarning: Could not complete secondary structure analysis: {str(e)}")
        print("This might be due to DSSP not being installed or other structural issues.")
        return None    
    
def analyze_structure(pdb_path: str) -> None:
    """
    Analyze structure composition and secondary structure elements.
    
    Args:
        pdb_path: Path to PDB file
    """
    try:
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=PDB.PDBExceptions.PDBConstructionWarning)
            
            parser = PDB.PDBParser(QUIET=True)
            structure = parser.get_structure('protein', pdb_path)
            
            # Analyze composition
            analyze_composition(structure)
            
            # Analyze secondary structure
            analyze_secondary_structure(structure)
            
    except Exception as e:
        print(f"\nWarning: Could not complete structure analysis: {str(e)}")

def renumber_consecutive_residues(pdb_path: str) -> str:
    """
    Renumber residues in a PDB file ensuring consecutive numbering.
    Waters and ions are stripped. Chain ID set to A.
    
    Args:
        pdb_path: Path to the input PDB file
    Returns:
        Path to the renumbered PDB file
    """
    try:
        output_path = f"{os.path.splitext(pdb_path)[0]}_renumd.pdb"
        from Bio import PDB
        import warnings
        from collections import defaultdict

        # Define residue categories
        AMINO_ACIDS = {
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE',
            'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
        }
        HEME_TYPES = {"HEM", "HEME", "HEA", "HEC"}
        WATER_TYPES = {"HOH", "WAT", "H2O"}
        ION_TYPES = {"NA", "CL", "K", "MG", "CA", "ZN", "FE", "MN"}

        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=PDB.PDBExceptions.PDBConstructionWarning)
            
            parser = PDB.PDBParser(QUIET=True)
            structure = parser.get_structure('protein', pdb_path)
            
            # Create new structure
            new_structure = PDB.Structure.Structure('protein')
            new_model = PDB.Model.Model(0)
            new_structure.add(new_model)
            new_chain = PDB.Chain.Chain('A')
            new_model.add(new_chain)

            # Categorize residues
            residue_categories = defaultdict(list)
            
            for model in structure:
                for chain in model:
                    for residue in chain:
                        resname = residue.resname
                        orig_info = (chain.id, residue.id[1], resname)
                        
                        # Categorize
                        if resname in AMINO_ACIDS:
                            category = "amino_acids"
                        elif resname in HEME_TYPES:
                            category = "hemes"
                        elif resname in WATER_TYPES or resname in ION_TYPES:
                            continue  # Skip waters and ions
                        else:
                            category = "unknown"
                        
                        residue_categories[category].append((orig_info, residue))

            # Renumber and save structure
            print("\nRenumbering residues (excluding waters and ions):")
            new_residue_number = 1
            
            # Process amino acids first
            if residue_categories["amino_acids"]:
                print("\nAmino acid residues:")
                print("-" * 50)
                for (orig_chain, orig_num, resname), residue in residue_categories["amino_acids"]:
                    new_residue = residue.copy()
                    new_residue.id = (residue.id[0], new_residue_number, residue.id[2])
                    new_chain.add(new_residue)
                    print(f"{orig_chain:>3} {orig_num:>4} {resname:<4} -> {new_residue_number:>4}")
                    new_residue_number += 1

            # Process hemes next
            if residue_categories["hemes"]:
                print("\nHeme residues:")
                print("-" * 50)
                for (orig_chain, orig_num, resname), residue in residue_categories["hemes"]:
                    new_residue = residue.copy()
                    new_residue.id = (residue.id[0], new_residue_number, residue.id[2])
                    new_chain.add(new_residue)
                    print(f"{orig_chain:>3} {orig_num:>4} {resname:<4} -> {new_residue_number:>4}")
                    new_residue_number += 1

            # Process unknown residues last
            if residue_categories["unknown"]:
                print("\nNon-standard residues:")
                print("-" * 50)
                for (orig_chain, orig_num, resname), residue in residue_categories["unknown"]:
                    new_residue = residue.copy()
                    new_residue.id = (residue.id[0], new_residue_number, residue.id[2])
                    new_chain.add(new_residue)
                    print(f"{orig_chain:>3} {orig_num:>4} {resname:<4} -> {new_residue_number:>4}")
                    new_residue_number += 1

            # Save the structure
            io = PDB.PDBIO()
            io.set_structure(new_structure)
            io.save(output_path)

            print(f"\nRenumbering complete!")
            print(f"All residues have been assigned to chain A and numbered consecutively from 1 to {new_residue_number-1}")

        return output_path

    except Exception as e:
        print(f"\nError renumbering residues: {e}")
        raise

def check_consecutive_residues(pdb_path: str) -> bool:
    """
    Check if residue numbers are consecutive by reading PDB file directly,
    preserving exact file order.
    
    Args:
        pdb_path: Path to the PDB file
    
    Returns:
        Boolean indicating if residue numbers are consecutive in file order
    """
    try:
        # Read PDB file directly to preserve exact order
        residues = []
        current_res = None
        
        with open(pdb_path, 'r') as f:
            for line in f:
                if line.startswith('ATOM  ') or line.startswith('HETATM'):
                    chain_id = line[21:22]
                    res_num = int(line[22:26])
                    res_name = line[17:20].strip()
                    
                    # Only add residue if it's different from the last one we saw
                    res_id = (chain_id, res_num, res_name)
                    if res_id != current_res:
                        current_res = res_id
                        residues.append(res_id)
        
        if not residues:
            print("\nNo residues found in the PDB file.")
            return False
            
        print("\nChecking residue numbering in exact PDB file order:")
        print("-" * 50)
        
        # Print residues in file order for verification
        for chain_id, res_num, resname in residues:
            print(f"{chain_id:>3} {res_num:>4} {resname:<4}")
            
        # Check if numbers are consecutive
        for i in range(len(residues)-1):
            curr_chain, curr_num, curr_name = residues[i]
            next_chain, next_num, next_name = residues[i+1]
            if next_num != curr_num + 1:
                print(f"\nGap in numbering:")
                print(f"Current residue: Chain {curr_chain} {curr_num} {curr_name}")
                print(f"Next residue: Chain {next_chain} {next_num} {next_name}")
                return False
                
        print(f"\nAll residues are consecutively numbered in PDB file order")
        print(f"from {residues[0][1]} to {residues[-1][1]}")
        return True
        
    except Exception as e:
        print(f"\nError checking residue numbering: {e}")
        return False

def select_pdb(launch_dir: Path, input_dict: Dict) -> str:
    """Select and validate PDB file."""
    interaction_manager = InteractionManager(
        launch_dir=launch_dir,
        input_dict=input_dict
    )

    # Check if we're running in non-interactive mode with an input file
    input_file_exists = (launch_dir / "input.txt").exists()

    # If we already have a PDB in the input dictionary and input.txt exists, use it directly
    if "OriginalPDB" in input_dict and input_file_exists:
        original_pdb = input_dict["OriginalPDB"]
        print(f"\n Using PDB from input file: {original_pdb}")

        pdb_filename = f"{original_pdb}.pdb"
        pdb_path = launch_dir / pdb_filename
        current_dir = Path.cwd()

        # If SPR directory wasn't pre-created, create it
        if not current_dir.name == 'SPR':
            spr_dir = launch_dir / 'SPR'
            spr_dir.mkdir(exist_ok=True)
            current_dir = spr_dir
            os.chdir(current_dir)

        if pdb_path.is_file():
            dest_path = current_dir / pdb_filename
            if not dest_path.is_file():
                shutil.copy2(pdb_path, dest_path)
                print(f"\n PDB copied to current working directory: {dest_path}")

            if check_consecutive_residues(str(dest_path)):
                print("\n PDB has consecutive residue numbering. No renumbering needed.")

                # Skip structure analysis in non-interactive mode
                if input_file_exists:
                    print("\n Skipping structure analysis in non-interactive mode.")
                else:
                    # Structure analysis on original structure
                    print("\nAnalyzing structure composition and secondary structure...")
                    analyze_structure(str(dest_path))

                return original_pdb

            print("\n Residues are not consecutively numbered. Attempting to renumber...")
            renumbered_path = renumber_consecutive_residues(str(dest_path))
            print(f"\n PDB renumbered and saved as: {renumbered_path}")

            # Update OriginalPDB to the new renumbered name
            input_dict["OriginalPDB"] = Path(renumbered_path).stem

            # Skip structure analysis in non-interactive mode
            if input_file_exists:
                print("\n Skipping structure analysis in non-interactive mode.")
            else:
                # Structure analysis on renumbered structure
                print("\nAnalyzing structure composition and secondary structure...")
                analyze_structure(str(renumbered_path))

            return Path(renumbered_path).stem
        else:
            # If in non-interactive mode but PDB doesn't exist, we should error out
            if input_file_exists:
                raise FileNotFoundError(f"PDB file specified in input.txt ({pdb_filename}) not found.")
            print("\n That PDB does not exist unfortunately")
            # Continue to interactive selection below

    # Interactive PDB selection (only reached if no valid PDB in input_dict or PDB not found)
    print(f"\n Good! Here are the PDBs in the launch directory ({launch_dir}) :\n")
    pdb_files = [x for x in os.listdir(launch_dir) if x.endswith(".pdb")]

    if not pdb_files:
        raise FileNotFoundError("No PDB files found in the current directory.")

    for x in pdb_files:
        print(x)

    while True:
        original_pdb = interaction_manager.prompt(
            "OriginalPDB",
            "\nWhich PDB would you like to setup (omit the .pdb file extension)?",
            allow_empty=False
        )

        pdb_filename = f"{original_pdb}.pdb"
        pdb_path = launch_dir / pdb_filename
        current_dir = Path.cwd()

        # If SPR directory wasn't pre-created, create it
        if not current_dir.name == 'SPR':
            spr_dir = launch_dir / 'SPR'
            spr_dir.mkdir(exist_ok=True)
            current_dir = spr_dir
            os.chdir(current_dir)

        if pdb_path.is_file():
            dest_path = current_dir / pdb_filename
            if not dest_path.is_file():
                shutil.copy2(pdb_path, dest_path)
                print(f"\n PDB copied to current working directory: {dest_path}")

            if check_consecutive_residues(str(dest_path)):
                print("\n PDB has consecutive residue numbering. No renumbering needed.")

                # Structure analysis on original structure
                print("\nAnalyzing structure composition and secondary structure...")
                analyze_structure(str(dest_path))

                return original_pdb

            print("\n Residues are not consecutively numbered. Attempting to renumber...")
            renumbered_path = renumber_consecutive_residues(str(dest_path))
            print(f"\n PDB renumbered and saved as: {renumbered_path}")

            # Update OriginalPDB to the new renumbered name
            input_dict["OriginalPDB"] = Path(renumbered_path).stem

            # Structure analysis on renumbered structure
            print("\nAnalyzing structure composition and secondary structure...")
            analyze_structure(str(renumbered_path))

            return Path(renumbered_path).stem
        else:
            print("\n That PDB does not exist unfortunately")

def initialize(launch_dir: Path, input_dict: Dict) -> str:
    """
    Main initialization function.
    Verifies programs are in PATH and selects PDB file.

    Args:
        launch_dir: Directory where script is launched
        input_dict: Dictionary of input parameters

    Returns:
        Base name of the selected PDB file
    """
    # Verify required programs are available
    if not verify_programs():
        sys.exit("Required programs are missing. Please install and retry.")

    # Select and validate PDB, getting updated input_dict
    original_pdb = select_pdb(launch_dir, input_dict)

    return original_pdb

