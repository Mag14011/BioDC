# biodc/core/preparation/process_residues.py

"""
Module for processing and renaming residues and atoms in PDB structures.
Handles heme environments, ligands, and associated structure modifications.
Replaces the VMD/TCL-based approach with BioPython implementation.
"""

import os
import sys
import logging
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, Any, Union
from dataclasses import dataclass
from Bio import PDB
from collections import defaultdict

from biodc.utils.interaction import InteractionManager

logger = logging.getLogger(__name__)

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
    heme_id: int
    shifted_id: int
    new_id: int
    his_p: int
    distal_ligand: int
    distal_type: str
    is_c_type: bool
    cys_b: Optional[int] = None
    cys_c: Optional[int] = None
    redox_state: str = 'ox'  # 'ox' or 'red'

@dataclass
class HemeNomenclature:
    """Defines naming patterns for heme and ligand residues."""
    
    # Heme naming patterns (c-type)
    CTYPE_HEME_NAMES = {
        'HIS': {'ox': 'HCO', 'red': 'HCR'},  # His-His
        'MET': {'ox': 'MCO', 'red': 'MCR'},  # His-Met
        'CYS': {'ox': 'CCO', 'red': 'CCR'},  # His-Cys
        'TYR': {'ox': 'YCO', 'red': 'YCR'},  # His-Tyr
        'ASP': {'ox': 'DCO', 'red': 'DCR'},  # His-Asp
        'GLU': {'ox': 'ECO', 'red': 'ECR'},  # His-Glu
        'ASN': {'ox': 'NCO', 'red': 'NCR'},  # His-Asn
        'GLN': {'ox': 'QCO', 'red': 'QCR'},  # His-Gln
        'LYS': {'ox': 'KCO', 'red': 'KCR'},  # His-Lys
    }

    # Heme naming patterns (b-type)
    BTYPE_HEME_NAMES = {
        'HIS': {'ox': 'HBO', 'red': 'HBR'},  # His-His
        'MET': {'ox': 'MBO', 'red': 'MBR'},  # His-Met
        'CYS': {'ox': 'CBO', 'red': 'CBR'},  # His-Cys
        'TYR': {'ox': 'YBO', 'red': 'YBR'},  # His-Tyr
        'ASP': {'ox': 'DBO', 'red': 'DBR'},  # His-Asp
        'GLU': {'ox': 'EBO', 'red': 'EBR'},  # His-Glu
        'ASN': {'ox': 'NBO', 'red': 'NBR'},  # His-Asn
        'GLN': {'ox': 'QBO', 'red': 'QBR'},  # His-Gln
        'LYS': {'ox': 'KBO', 'red': 'KBR'},  # His-Lys
    }

    # Proximal His naming (c-type - P prefix)
    CTYPE_PROX_HIS = {
        'HIS': {'ox': 'PHO', 'red': 'PHR'},  # For His-His
        'MET': {'ox': 'PMO', 'red': 'PMR'},  # For His-Met
        'CYS': {'ox': 'PCO', 'red': 'PCR'},  # For His-Cys
        'TYR': {'ox': 'PYO', 'red': 'PYR'},  # For His-Tyr
        'ASP': {'ox': 'PDO', 'red': 'PDR'},  # For His-Asp
        'GLU': {'ox': 'PEO', 'red': 'PER'},  # For His-Glu
        'ASN': {'ox': 'PNO', 'red': 'PNR'},  # For His-Asn
        'GLN': {'ox': 'PQO', 'red': 'PQR'},  # For His-Gln
        'LYS': {'ox': 'PKO', 'red': 'PKR'},  # For His-Lys
    }
    
    # Proximal His naming (b-type - F prefix)
    BTYPE_PROX_HIS = {
        'HIS': {'ox': 'FHO', 'red': 'FHR'},  # For His-His
        'MET': {'ox': 'FMO', 'red': 'FMR'},  # For His-Met
        'CYS': {'ox': 'FCO', 'red': 'FCR'},  # For His-Cys
        'TYR': {'ox': 'FYO', 'red': 'FYR'},  # For His-Tyr
        'ASP': {'ox': 'FDO', 'red': 'FDR'},  # For His-Asp
        'GLU': {'ox': 'FEO', 'red': 'FER'},  # For His-Glu
        'ASN': {'ox': 'FNO', 'red': 'FNR'},  # For His-Asn
        'GLN': {'ox': 'FQO', 'red': 'FQR'},  # For His-Gln
        'LYS': {'ox': 'FKO', 'red': 'FKR'},  # For His-Lys
    }
    
    # Distal ligand naming (c-type - D prefix)
    CTYPE_DISTAL = {
        'HIS': {'ox': 'DHO', 'red': 'DHR'},  # For His-His
        'MET': {'ox': 'DMO', 'red': 'DMR'},  # For His-Met
        'CYS': {'ox': 'DCO', 'red': 'DCR'},  # For His-Cys
        'TYR': {'ox': 'DYO', 'red': 'DYR'},  # For His-Tyr
        'ASP': {'ox': 'DDO', 'red': 'DDR'},  # For His-Asp
        'GLU': {'ox': 'DEO', 'red': 'DER'},  # For His-Glu
        'ASN': {'ox': 'DNO', 'red': 'DNR'},  # For His-Asn
        'GLN': {'ox': 'DQO', 'red': 'DQR'},  # For His-Gln
        'LYS': {'ox': 'DKO', 'red': 'DKR'},  # For His-Lys
    }
    
    # Distal ligand naming (b-type - R prefix)
    BTYPE_DISTAL = {
        'HIS': {'ox': 'RHO', 'red': 'RHR'},  # For His-His
        'MET': {'ox': 'RMO', 'red': 'RMR'},  # For His-Met
        'CYS': {'ox': 'RCO', 'red': 'RCR'},  # For His-Cys
        'TYR': {'ox': 'RYO', 'red': 'RYR'},  # For His-Tyr
        'ASP': {'ox': 'RDO', 'red': 'RDR'},  # For His-Asp
        'GLU': {'ox': 'REO', 'red': 'RER'},  # For His-Glu
        'ASN': {'ox': 'RNO', 'red': 'RNR'},  # For His-Asn
        'GLN': {'ox': 'RQO', 'red': 'RQR'},  # For His-Gln
        'LYS': {'ox': 'RKO', 'red': 'RKR'},  # For His-Lys
    }

class PDBProcessor:
    def __init__(self, pdb_path: Path):
        """Initialize with PDB file path."""
        print(f"\nInitializing PDB Processor for: {pdb_path}")
        self.parser = PDB.PDBParser(QUIET=True)
        self.structure = self.parser.get_structure('protein', pdb_path)
        self.nomenclature = HemeNomenclature()
        self.heme_shifts = {}  # Maps original heme IDs to new IDs
        print("Structure loaded successfully")

    def calculate_heme_shifts(self, indexing_file: Path) -> Dict[int, int]:
        """First pass: Calculate new positions for all hemes."""
        print("\nCalculating heme shifts from ResIndexing.txt...")
        
        # Store both the shifts for moving hemes and original->new mapping
        self.heme_shifts = {}  # For shifting hemes
        self.id_mapping = {}   # For updating ResIndexing.txt

        heme_count = 0
        with open(indexing_file) as f:
            for line_num, line in enumerate(f, 1):
                parts = line.strip().split()
                if not parts:
                    continue
                    
                # Handle both b-type (5 fields) and c-type (7 fields)
                if len(parts) == 5:
                    original_id = int(parts[2])
                    heme_type = "b-type"
                elif len(parts) == 7:
                    original_id = int(parts[4])
                    heme_type = "c-type"
                else:
                    print(f"  Skipping invalid line {line_num}: {line.strip()}")
                    continue
                
                # Calculate new position
                new_id = original_id + (2 * heme_count)
                
                # Store both mappings
                self.heme_shifts[original_id] = new_id
                self.id_mapping[original_id] = new_id
                
                print(f"  Line {line_num}: {heme_type} heme {original_id} → {new_id} "
                    f"(shift: +{2 * heme_count})")
                heme_count += 1
        
        print(f"Found {heme_count} hemes to process")
        return self.heme_shifts

    def shift_hemes(self):
        """First pass: Move all hemes to their new positions."""
        print("\nExecuting first pass: Shifting heme positions...")
        shifted_count = 0

        # Get all hemes and sort them by ID in REVERSE order
        hemes_to_shift = []
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    original_id = residue.id[1]
                    if original_id in self.heme_shifts:
                        hemes_to_shift.append(residue)

        # Sort in reverse order to process highest IDs first
        hemes_to_shift.sort(key=lambda x: x.id[1], reverse=True)

        # Process hemes in reverse order
        for residue in hemes_to_shift:
            original_id = residue.id[1]
            new_id = self.heme_shifts[original_id]
            old_name = residue.resname
            print(f"  Moving heme {original_id} ({old_name}) → position {new_id}")
            residue.id = (residue.id[0], new_id, residue.id[2])
            shifted_count += 1

        print(f"Completed first pass: {shifted_count} hemes shifted")

    def process_disulfides(self, disulf_list: str):
        """Rename disulfide-bonded cysteines to CYX."""
        if not disulf_list:
            print("\nNo disulfide bonds to process")
            return
            
        print("\nProcessing disulfide bonds...")
        disulf_ids = [int(x) for x in disulf_list.split()]
        processed_count = 0
        
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    if (residue.get_resname() in ['CYS', 'CYX'] and 
                        residue.id[1] in disulf_ids):
                        print(f"  Renaming CYS {residue.id[1]} → CYX")
                        residue.resname = 'CYX'
                        processed_count += 1
        
        print(f"Processed {processed_count} disulfide-bonded cysteines")

    def process_titratable_residues(self, asp_ids: str, glu_ids: str, his_ids: str):
        """Rename titratable residues for constant pH simulation."""
        print("\nProcessing titratable residues...")
        
        id_mappings = {
            'ASP': (asp_ids, 'AS4'),
            'GLU': (glu_ids, 'GL4'),
            'HIS': (his_ids, 'HIP')
        }
        
        total_processed = 0
        for res_type, (id_list, new_name) in id_mappings.items():
            if not id_list:
                continue
                
            res_ids = [int(x) for x in id_list.split()]
            print(f"\n  Processing {res_type} residues...")
            processed = 0
            
            for model in self.structure:
                for chain in model:
                    for residue in chain:
                        if (residue.get_resname() == res_type and 
                            residue.id[1] in res_ids):
                            print(f"    {res_type} {residue.id[1]} → {new_name}")
                            residue.resname = new_name
                            processed += 1
            
            print(f"  Processed {processed} {res_type} residues")
            total_processed += processed
        
        print(f"Total titratable residues processed: {total_processed}")

    def rename_atoms(self, residue: PDB.Residue, old_name: str, new_name: str):
        """Rename specific atoms in a residue."""
        if old_name in residue:
            residue[old_name].id = new_name
            residue[old_name].fullname = f" {new_name} "
            print(f"    Renamed atom {old_name} → {new_name}")

    def rename_propionate_atoms(self, residue: PDB.Residue):
        """Rename propionate atoms according to schema."""
        print(f"  Renaming propionate atoms in residue {residue.id[1]}")
        atom_mappings = {
            'CAA': 'CA', 'CAD': 'CA',
            'CBA': 'CB', 'CBD': 'CB',
            'CGA': 'CG', 'CGD': 'CG',
            'O1A': 'O1', 'O1D': 'O1',
            'O2A': 'O2', 'O2D': 'O2'
        }
        
        renamed_count = 0
        for old_name, new_name in atom_mappings.items():
            if old_name in residue:
                atom = residue[old_name]
                atom.id = new_name
                atom.fullname = f" {new_name} "
                renamed_count += 1
                print(f"    {old_name} → {new_name}")
        
        print(f"  Renamed {renamed_count} propionate atoms")

    def rename_cys_atoms(self, residue: PDB.Residue, position: str, heme_name: str, original_heme_id: int):
        """
        Rename Cys atoms for heme attachment, ensuring atoms are fully detached.

        Args:
            residue (PDB.Residue): The Cysteine residue to process
            position (str): 'B' or 'C' indicating the heme binding position
            heme_name (str): Name of the associated heme (e.g., 'HCR')
            original_heme_id (int): Original heme ID from ResIndexing.txt
        """
        print(f"\n  Renaming Cys atoms for position {position} in residue {residue.id[1]}")
        print(f"  Associated Heme Name: {heme_name}")
        print(f"  Original Heme ID: {original_heme_id}")
        print(f"  ID Mapping: {self.id_mapping}")

        # Rename backbone atoms to CYO
        backbone_atoms = ['N', 'CA', 'C', 'O']
        print("  Renaming backbone atoms:")
        for atom_name in backbone_atoms:
            if atom_name in residue:
                print(f"    Preserving backbone atom {atom_name}")
                self.rename_atoms(residue, atom_name, atom_name)

        print(f"  Changing residue name to CYO")
        residue.resname = 'CYO'

        # Simplified atom naming based on heme name
        cys_atom_map = {
            'B': {'CB': 'CBB2', 'SG': 'SGB2'},
            'C': {'CB': 'CBC1', 'SG': 'SGC1'}
        }

        # Find the target heme using the original heme ID
        target_heme_id = self.id_mapping.get(original_heme_id)

        if target_heme_id is not None:
            target_heme = self.get_residue_by_id(target_heme_id)
            if target_heme:
                print(f"  Found target heme: {target_heme.get_resname()} at ID {target_heme_id}")

                # Detailed atom transfer with error handling
                if position in cys_atom_map:
                    for old_name, new_name in cys_atom_map[position].items():
                        if old_name in residue:
                            atom = residue[old_name]

                            try:
                                # Detailed logging of atom properties
#                               print(f"    Atom {old_name} details:")
#                               print(f"      Current parent: {atom.parent}")
#                               print(f"      ID: {atom.id}")
#                               print(f"      Fullname: {atom.fullname}")
#                               print(f"      Coordinates: {atom.coord}")

                                # Update atom properties
                                atom.id = new_name
                                atom.fullname = f" {new_name} "

                                # Attempt to add atom to heme
                                try:
                                    target_heme.add(atom)
                                    print(f"    Successfully moved {new_name} to {target_heme.get_resname()}")
                                except Exception as add_error:
                                    print(f"    ERROR adding atom to heme: {add_error}")
                                    # Attempt to manually add via parent's method
                                    try:
                                        parent = target_heme.parent
                                        if parent:
                                            parent.add(atom)
                                            print(f"    Manually added {new_name} via parent")
                                    except Exception as parent_error:
                                        print(f"    CRITICAL ERROR adding atom: {parent_error}")

                                # Remove from original residue
                                try:
                                    residue.detach_child(old_name)
                                    print(f"    Removed {old_name} from original residue")
                                except Exception as detach_error:
                                    print(f"    ERROR detaching atom: {detach_error}")

                            except Exception as atom_error:
                                print(f"    CRITICAL ERROR processing atom {old_name}: {atom_error}")
            else:
                print(f"  Warning: Could not find target heme with ID {target_heme_id}")
        else:
            print(f"  Warning: No mapping found for original heme ID {original_heme_id}")

    def create_propionate_residue(self, heme: PDB.Residue, prop_id: str, new_id: int) -> PDB.Residue:
        """
        Create a new PRN residue from propionate atoms.
        Moves (not copies) the atoms from the heme residue to the new PRN residue.
        
        Args:
            heme: Source heme residue containing the propionate atoms
            prop_id: Propionate identifier ('A' or 'D')
            new_id: New residue ID for the PRN residue
            
        Returns:
            New PRN residue containing the moved propionate atoms
        """
        print(f"  Creating PRN residue {new_id} from heme {heme.id[1]} propionate {prop_id}")
        new_residue = PDB.Residue.Residue((' ', new_id, ' '), 'PRN', '')
        
        # Define the atoms to move
        atom_mappings = {
            f'CA{prop_id}': 'CA',
            f'CB{prop_id}': 'CB', 
            f'CG{prop_id}': 'CG',
            f'O1{prop_id}': 'O1',
            f'O2{prop_id}': 'O2'
        }
        
        transferred_atoms = 0
        atoms_to_detach = []  # Store atoms to be removed from heme
        
        # First pass: Add atoms to new residue
        for old_name, new_name in atom_mappings.items():
            if old_name in heme:
                atom = heme[old_name]
                # Update atom properties for new residue
                atom.id = new_name
                atom.fullname = f" {new_name} "
                atom.parent = None  # Detach from old parent
                
                # Add to new residue
                new_residue.add(atom)
                atoms_to_detach.append(old_name)
                transferred_atoms += 1
                print(f"    Transferred atom {old_name} → {new_name}")
        
        # Second pass: Remove atoms from heme residue
        for atom_name in atoms_to_detach:
            # Use detach_child to properly remove the atom
            heme.detach_child(atom_name)
            print(f"    Removed {atom_name} from heme residue")
        
        print(f"  Transferred {transferred_atoms} atoms to PRN residue {new_id}")
        return new_residue

    def get_residue_by_id(self, res_id: int) -> Optional[Any]:
        """Get a residue by its ID."""
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    if residue.id[1] == res_id:
                        return residue
        return None

    def process_heme_environment(self, env: HemeEnvironment):
        """Process a complete heme environment with proper propionate handling."""
        print(f"\nProcessing heme environment:")
        print(f"  Original heme ID: {env.heme_id}")
        print(f"  Type: {'c-type' if env.is_c_type else 'b-type'}")
        print(f"  Proximal His: {env.his_p}")
        print(f"  Distal ligand: {env.distal_ligand} ({env.distal_type})")
        if env.is_c_type:
            print(f"  Cys B/C: {env.cys_b}/{env.cys_c}")
        print(f"  Redox state: {env.redox_state}")
        
        # Get the shifted heme ID
        new_heme_id = self.heme_shifts[env.heme_id]
        print(f"  New heme position: {new_heme_id}")
        
        # Find the shifted heme
        heme = self.get_residue_by_id(new_heme_id)
        if not heme:
            raise ValueError(f"Could not find shifted heme at position {new_heme_id}")
        
        # Get new names based on environment
        heme_names = (self.nomenclature.CTYPE_HEME_NAMES if env.is_c_type 
                    else self.nomenclature.BTYPE_HEME_NAMES)
        new_heme_name = heme_names[env.distal_type][env.redox_state]
        print(f"  Renaming heme to: {new_heme_name}")
        heme.resname = new_heme_name
        
        # Create PRN residues - ensure atoms are moved, not copied
        print("\n  Creating propionate residues:")
        parent_chain = heme.parent
        
        # Process propionate A
        prn_a = self.create_propionate_residue(heme, 'A', new_heme_id + 1)
        parent_chain.add(prn_a)
        print(f"  Verified propionate A transfer - atoms removed from heme")
        
        # Process propionate D
        prn_d = self.create_propionate_residue(heme, 'D', new_heme_id + 2)
        parent_chain.add(prn_d)
        print(f"  Verified propionate D transfer - atoms removed from heme")
        
        # Process ligands
        print("\n  Processing ligands:")
        if env.his_p:
            prox_his = self.get_residue_by_id(env.his_p)
            if prox_his:
                prox_names = (self.nomenclature.CTYPE_PROX_HIS if env.is_c_type 
                            else self.nomenclature.BTYPE_PROX_HIS)
                new_name = prox_names[env.distal_type][env.redox_state]
                print(f"    Proximal His {env.his_p} → {new_name}")
                prox_his.resname = new_name
        
        if env.distal_ligand:
            distal = self.get_residue_by_id(env.distal_ligand)
            if distal:
                distal_names = (self.nomenclature.CTYPE_DISTAL if env.is_c_type 
                            else self.nomenclature.BTYPE_DISTAL)
                new_name = distal_names[env.distal_type][env.redox_state]
                print(f"    Distal ligand {env.distal_ligand} → {new_name}")
                distal.resname = new_name


        # Handle c-type specific processing
        if env.is_c_type and env.cys_b and env.cys_c:
            print("\n  Processing c-type specific modifications:")
            cys_b = self.get_residue_by_id(env.cys_b)
            cys_c = self.get_residue_by_id(env.cys_c)

            # Get heme name for renaming
            heme_names = (self.nomenclature.CTYPE_HEME_NAMES if env.is_c_type
                        else self.nomenclature.BTYPE_HEME_NAMES)
            heme_name = heme_names[env.distal_type][env.redox_state]

            if cys_b:
                print(f"    Processing Cys B ({env.cys_b}):")
                self.rename_cys_atoms(cys_b, 'B', heme_name, env.heme_id)

            if cys_c:
                print(f"    Processing Cys C ({env.cys_c}):")
                self.rename_cys_atoms(cys_c, 'C', heme_name, env.heme_id)

                print("  Environment processing complete")

    def save_structure(self, output_path: Path):
        """Save the modified structure with ordered residues and TER records."""
        print("\nSaving structure with ordered residues...")

        # Ensure output is in the current working directory
        output_path = Path.cwd() / output_path.name
        print(f"Output path: {output_path}")

        # Get all residues and sort by ID
        all_residues = []
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    all_residues.append(residue)

        # Sort by residue ID
        all_residues.sort(key=lambda r: r.id[1])

        # Create new structure with sorted residues
        new_structure = PDB.Structure.Structure('protein')
        new_model = PDB.Model.Model(0)
        new_structure.add(new_model)
        new_chain = PDB.Chain.Chain('A')
        new_model.add(new_chain)

        # Add sorted residues - use deep copy to preserve all modifications
        for residue in all_residues:
            # Create a deep copy of the residue to ensure all modifications are preserved
            new_residue = residue.copy()
            new_chain.add(new_residue)

        # First save the sorted structure
        temp_output = output_path.with_suffix('.temp.pdb')
        io = PDB.PDBIO()
        io.set_structure(new_structure)
        io.save(str(temp_output))

        # Then add TER records
        self.add_ter_records_to_pdb(temp_output, output_path)
        
        # Remove temporary file
        temp_output.unlink()

        print(f"\nStructure saved successfully to {output_path}")

    def validate_structure(self, indexing_file: Path) -> bool:
        """Validate the processed structure."""
        print("\n")
        print("=" * 60)
        print("Validating processed structure...")
        print("=" * 60)

        # Track all residue IDs and their types
        used_ids = set()
        heme_prn_groups = []
        cys_issues = []
        
        # Read bound Cysteine residues from indexing file
        bound_cys = set()
        with open(indexing_file, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) == 7:  # c-type heme line
                    try:
                        cys_b = int(parts[0])
                        cys_c = int(parts[1])
                        bound_cys.update([cys_b, cys_c])
                    except (ValueError, IndexError):
                        print(f"Warning: Could not parse line: {line.strip()}")
                elif len(parts) == 5:  # b-type heme line
                    if parts[4] == 'CYS':  # Check for His-Cys ligation
                        try:
                            distal_cys = int(parts[1])  # Distal ligand ID
                            bound_cys.add(distal_cys)
                        except (ValueError, IndexError):
                            print(f"Warning: Could not parse distal Cys from line: {line.strip()}")
        
        # Collect all residue IDs and heme-PRN groups
        print("\nChecking residue organization...")
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    res_id = residue.id[1]
                    res_name = residue.resname
                    
                    # Check for duplicate IDs
                    if res_id in used_ids:
                        print(f"Error: Duplicate residue ID found: {res_id}")
                        return False
                    used_ids.add(res_id)
                    
                    # Validate Cys residues bound to hemes
                    if res_id in bound_cys:
                        if res_name != 'CYO':
                            cys_issues.append(f"Residue {res_id} not renamed to CYO (current: {res_name})")
                        
                        sidechain_atoms = [atom for atom in residue if atom.id in ['CB', 'SG', 'CBB2', 'SGB2', 'CBC1', 'SGC1']]
                        if sidechain_atoms:
                            cys_issues.append(f"Residue {res_id} still has sidechain atoms: {[atom.id for atom in sidechain_atoms]}")
                    
                    # Track heme and PRN groups
                    if any(name in res_name for name in ['HBO', 'HBR', 'HCO', 'HCR', 'MBO', 'MBR', 'MCO', 'MCR']):
                        print(f"  Found heme residue: {res_id} ({res_name})")
                        heme_prn_groups.append({'heme_id': res_id, 'prn_ids': []})
                        # Look ahead for the next two PRN residues
                        expected_prn1 = res_id + 1
                        expected_prn2 = res_id + 2
                        for prn_check in chain:
                            if prn_check.id[1] in [expected_prn1, expected_prn2] and prn_check.resname == 'PRN':
                                heme_prn_groups[-1]['prn_ids'].append(prn_check.id[1])
                                print(f"  Associated PRN {prn_check.id[1]} with heme {res_id}")

        # Print Cys validation issues
        if cys_issues:
            print("\nCysteine Residue Validation Issues:")
            for issue in cys_issues:
                print(f"  - {issue}")
            return False

        # Validate heme spacing and PRN assignments
        print("\nValidating heme-PRN relationships...")
        for i, group in enumerate(heme_prn_groups):
            # Check heme spacing
            if i > 0:
                expected_spacing = heme_prn_groups[i]['heme_id'] - heme_prn_groups[i-1]['heme_id']
                if expected_spacing != 3:
                    print(f"Error: Incorrect heme spacing between hemes {heme_prn_groups[i-1]['heme_id']} "
                        f"and {group['heme_id']} (spacing={expected_spacing}, should be 3)")
                    return False
                print(f"  Verified spacing between hemes {heme_prn_groups[i-1]['heme_id']} "
                    f"and {group['heme_id']}")
            
            # Check PRN residues
            if len(group['prn_ids']) != 2:
                print(f"Error: Heme {group['heme_id']} has incorrect number of PRN residues: "
                    f"{len(group['prn_ids'])} (should be 2)")
                return False
            
            # Check PRN spacing
            expected_prns = [group['heme_id'] + 1, group['heme_id'] + 2]
            if group['prn_ids'] != expected_prns:
                print(f"Error: Incorrect PRN residue IDs for heme {group['heme_id']}: "
                    f"{group['prn_ids']} (should be {expected_prns})")
                return False
            print(f"  Verified PRN residues {group['prn_ids']} for heme {group['heme_id']}")

        # Compare expected vs actual heme counts
        heme_counts = self.count_heme_types()
        expected_counts = defaultdict(lambda: {'ox': 0, 'red': 0})
        
        # Read expected counts from ResIndexing.txt
        with open(indexing_file) as f:
            for line in f:
                parts = line.strip().split()
                if not parts:
                    continue
                    
                if len(parts) == 5:  # b-type: HisP Distal OriginalHemeID b HX
                    key = f"{parts[4]}-{parts[3]}"  # Results in "HH-b", "HM-b", etc.
                    expected_counts[key]['red'] += 1
                elif len(parts) == 7:  # c-type
                    key = f"{parts[6]}-{parts[5]}"  # Results in "HH-c", "HM-c", etc.
                    expected_counts[key]['red'] += 1
        
        # Print comparison
        print("\nDetailed heme counts comparison:")
        print("-" * 80)
        print(f"{'Heme Type':<10} {'Expected':<8} {'Total':<8} {'Oxidized':<8} {'Reduced':<8} {'Match?':<8}")
        print("-" * 80)
        
        all_match = True
        for heme_type in sorted(set(list(heme_counts.keys()) + list(expected_counts.keys()))):
            expected_total = sum(expected_counts[heme_type].values())
            actual_total = sum(heme_counts[heme_type].values())
            actual_ox = heme_counts[heme_type]['ox']
            actual_red = heme_counts[heme_type]['red']
            matches = expected_total == actual_total
            if not matches:
                all_match = False
                
            print(f"{heme_type:<10} {expected_total:<8} {actual_total:<8} "
                f"{actual_ox:<8} {actual_red:<8} {'✓' if matches else '✗':<8}")
        
        print("-" * 80)
        if all_match:
            print("\nValidation PASSED: All heme counts match specifications.")
        else:
            print("\nValidation FAILED: Discrepancies found in heme counts.")
            print("Note: 'Total' should match 'Expected' for each type. "
                "'Oxidized' + 'Reduced' = 'Total'")

        return all_match

#   def validate_structure(self, indexing_file: Path) -> bool:
#       """Validate the processed structure."""
#       print("\nValidating processed structure...")

#       # Track all residue IDs and their types
#       used_ids = set()
#       heme_prn_groups = []
#       cys_issues = []
#       
#       # Read bound Cysteine residues from indexing file
#       bound_cys = set()
#       with open(indexing_file, 'r') as f:
#           for line in f:
#               parts = line.strip().split()
#               if len(parts) == 7:  # c-type heme line
#                   try:
#                       cys_b = int(parts[0])
#                       cys_c = int(parts[1])
#                       bound_cys.update([cys_b, cys_c])
#                   except (ValueError, IndexError):
#                       print(f"Warning: Could not parse line: {line.strip()}")
#       
#       # Collect all residue IDs and heme-PRN groups
#       print("\nChecking residue organization...")
#       for model in self.structure:
#           for chain in model:
#               for residue in chain:
#                   res_id = residue.id[1]
#                   res_name = residue.resname
#                   
#                   # Check for duplicate IDs
#                   if res_id in used_ids:
#                       print(f"Error: Duplicate residue ID found: {res_id}")
#                       return False
#                   used_ids.add(res_id)
#                   
#                   # Validate Cys residues bound to hemes
#                   if res_id in bound_cys:
#                       # Check if it's renamed to CYO
#                       if res_name != 'CYO':
#                           cys_issues.append(f"Residue {res_id} not renamed to CYO (current: {res_name})")
#                       
#                       # Check that sidechain is removed
#                       sidechain_atoms = [atom for atom in residue if atom.id in ['CB', 'SG', 'CBB2', 'SGB2', 'CBC1', 'SGC1']]
#                       if sidechain_atoms:
#                           cys_issues.append(f"Residue {res_id} still has sidechain atoms: {[atom.id for atom in sidechain_atoms]}")
#                   
#                   # Track heme and PRN groups
#                   if any(name in res_name for name in ['HBO', 'HBR', 'HCO', 'HCR', 
#                                                   'MBO', 'MBR', 'MCO', 'MCR']):
#                       print(f"  Found heme residue: {res_id} ({res_name})")
#                       heme_prn_groups.append({'heme_id': res_id, 'prn_ids': []})
#                   elif res_name == 'PRN':
#                       if heme_prn_groups:
#                           heme_prn_groups[-1]['prn_ids'].append(res_id)
#                           print(f"  Associated PRN {res_id} with heme {heme_prn_groups[-1]['heme_id']}")
#       
#       # Print Cys validation issues
#       if cys_issues:
#           print("\nCysteine Residue Validation Issues:")
#           for issue in cys_issues:
#               print(f"  - {issue}")
#           return False

#       # Validate heme spacing and PRN assignments
#       print("\nValidating heme-PRN relationships...")
#       for i, group in enumerate(heme_prn_groups):
#           # Check heme spacing
#           if i > 0:
#               expected_spacing = heme_prn_groups[i]['heme_id'] - heme_prn_groups[i-1]['heme_id']
#               if expected_spacing != 3:
#                   print(f"Error: Incorrect heme spacing between hemes {heme_prn_groups[i-1]['heme_id']} "
#                         f"and {group['heme_id']} (spacing={expected_spacing}, should be 3)")
#                   return False
#               print(f"  Verified spacing between hemes {heme_prn_groups[i-1]['heme_id']} "
#                     f"and {group['heme_id']}")
#           
#           # Check PRN residues
#           if len(group['prn_ids']) != 2:
#               print(f"Error: Heme {group['heme_id']} has incorrect number of PRN residues: "
#                     f"{len(group['prn_ids'])} (should be 2)")
#               return False
#           
#           # Check PRN spacing
#           expected_prns = [group['heme_id'] + 1, group['heme_id'] + 2]
#           if group['prn_ids'] != expected_prns:
#               print(f"Error: Incorrect PRN residue IDs for heme {group['heme_id']}: "
#                     f"{group['prn_ids']} (should be {expected_prns})")
#               return False
#           print(f"  Verified PRN residues {group['prn_ids']} for heme {group['heme_id']}")
#       
#       print("\nStructure validation complete - all checks passed!")
#       return True

    def process_structure(self, indexing_file: Path):
        """Main processing function implementing the two-pass approach."""
        print("\nStarting two-pass heme processing...")
        
        # First pass: Calculate and apply heme shifts
        self.calculate_heme_shifts(indexing_file)
        self.shift_hemes()
        
        # Second pass: Process environments
        print("\nStarting second pass: Processing heme environments...")
        
        # Validate final structure
        return self.validate_structure(indexing_file)


    def count_heme_types(self) -> Dict[str, Dict[str, int]]:
        """Count heme types in structure, with oxidation states."""
        heme_counts = defaultdict(lambda: {'ox': 0, 'red': 0})
        
        # Comprehensive mapping of heme names to their type and redox state
        heme_name_map = {
            # His-His ligated
            'HBO': ('HH-b', 'ox'), 'HBR': ('HH-b', 'red'),
            'HCO': ('HH-c', 'ox'), 'HCR': ('HH-c', 'red'),
            # His-Met ligated
            'MBO': ('HM-b', 'ox'), 'MBR': ('HM-b', 'red'),
            'MCO': ('HM-c', 'ox'), 'MCR': ('HM-c', 'red'),
            # His-Cys ligated
            'CBO': ('HC-b', 'ox'), 'CBR': ('HC-b', 'red'),
            'CCO': ('HC-c', 'ox'), 'CCR': ('HC-c', 'red'),
            # His-Tyr ligated
            'YBO': ('HY-b', 'ox'), 'YBR': ('HY-b', 'red'),
            'YCO': ('HY-c', 'ox'), 'YCR': ('HY-c', 'red'),
            # His-Asp ligated
            'DBO': ('HD-b', 'ox'), 'DBR': ('HD-b', 'red'),
            'DCO': ('HD-c', 'ox'), 'DCR': ('HD-c', 'red'),
            # His-Glu ligated
            'EBO': ('HE-b', 'ox'), 'EBR': ('HE-b', 'red'),
            'ECO': ('HE-c', 'ox'), 'ECR': ('HE-c', 'red'),
            # His-Asn ligated
            'NBO': ('HN-b', 'ox'), 'NBR': ('HN-b', 'red'),
            'NCO': ('HN-c', 'ox'), 'NCR': ('HN-c', 'red'),
            # His-Gln ligated
            'QBO': ('HQ-b', 'ox'), 'QBR': ('HQ-b', 'red'),
            'QCO': ('HQ-c', 'ox'), 'QCR': ('HQ-c', 'red'),
            # His-Lys ligated
            'KBO': ('HK-b', 'ox'), 'KBR': ('HK-b', 'red'),
            'KCO': ('HK-c', 'ox'), 'KCR': ('HK-c', 'red')
        }

        for model in self.structure:
            for chain in model:
                for residue in chain:
                    res_name = residue.get_resname()
                    if res_name in heme_name_map:
                        heme_type, redox_state = heme_name_map[res_name]
                        heme_counts[heme_type][redox_state] += 1

        return heme_counts

#   def count_heme_types(self) -> Dict[str, Dict[str, int]]:
#       """
#       Count the number of each type of heme in the processed structure,
#       breaking down by oxidized and reduced states.
#       Returns a nested dictionary with counts of each heme type and redox state.
#       """
#       heme_counts = defaultdict(lambda: {'ox': 0, 'red': 0})
#       
#       # Mapping of heme names to their type and redox state
#       heme_name_map = {
#           'HBO': ('HH-b', 'ox'), 'HBR': ('HH-b', 'red'),
#           'HCO': ('HH-c', 'ox'), 'HCR': ('HH-c', 'red'),
#           'MBO': ('HM-b', 'ox'), 'MBR': ('HM-b', 'red'),
#           'MCO': ('HM-c', 'ox'), 'MCR': ('HM-c', 'red')
#       }
#       
#       for model in self.structure:
#           for chain in model:
#               for residue in chain:
#                   res_name = residue.get_resname()
#                   
#                   if res_name in heme_name_map:
#                       heme_type, redox_state = heme_name_map[res_name]
#                       heme_counts[heme_type][redox_state] += 1
#       
#       return heme_counts

    def add_ter_records_to_pdb(self, input_pdb: Union[str, Path],
                            output_pdb: Optional[Union[str, Path]] = None) -> None:
        """
        Add TER records to a PDB file with minimal manipulation.
        
        Args:
            input_pdb (str or Path): Path to input PDB file
            output_pdb (str or Path, optional): Path to output PDB file.
                                            If None, overwrites the input file.
        """
        # Convert to Path objects
        input_pdb = Path(input_pdb)

        # If no output_pdb specified, use input_pdb
        if output_pdb is None:
            output_pdb = input_pdb
        output_pdb = Path(output_pdb)

        # Read the original PDB file
        with open(input_pdb, 'r') as f:
            lines = f.readlines()

        # Combine heme names from nomenclature object
        all_heme_names = set()
        for ligand_dict in [self.nomenclature.BTYPE_HEME_NAMES, 
                        self.nomenclature.CTYPE_HEME_NAMES]:
            for variant_dict in ligand_dict.values():
                all_heme_names.update(variant_dict.values())

        # Prepare new lines
        new_lines = []
        prev_residue_id = None
        prev_residue_name = None

        for i, line in enumerate(lines):
            # Add current line
            new_lines.append(line)

            # Only process ATOM and HETATM lines
            if not line.startswith(('ATOM', 'HETATM')):
                continue

            # Extract current residue information
            curr_residue_id = int(line[22:26])
            curr_residue_name = line[17:20].strip()

            # Determine if this is a terminus
            is_terminus = False

            # 1. OXT criteria (highest priority)
            if 'OXT' in line[12:16]:
                is_terminus = True

            # 2. Next line (if exists) is a different residue
            if i + 1 < len(lines):
                next_line = lines[i + 1]
                if next_line.startswith(('ATOM', 'HETATM')):
                    next_residue_id = int(next_line[22:26])
                    next_residue_name = next_line[17:20].strip()

                    # If next residue is different
                    if curr_residue_id != next_residue_id:
                        # Check if next residue is a heme
                        if next_residue_name in all_heme_names:
                            is_terminus = True

            # Add TER record if it's a terminus
            if is_terminus:
                ter_line = "TER\n"
                new_lines.append(ter_line)

        # Ensure a final TER record if none exists
        if not any(line.startswith('TER') for line in new_lines):
            new_lines.append("TER\n")

        # Write modified PDB
        with open(output_pdb, 'w') as f:
            f.writelines(new_lines)

@dataclass
class ForceFieldRequirement:
    """Required forcefield files for a heme type."""
    lib_file: str
    frcmod_file: str
    description: str
    implemented: bool = False  # Flag to indicate if parameters are developed

class ForceFieldValidator:
    """Validates availability of required forcefield files."""
    
    def __init__(self, forcefield_dir: Path):
        self.forcefield_dir = Path(forcefield_dir)
        # Define required files for each heme type and state
        self.requirements = {
            # b-type hemes
            # His-His (implemented)
            ('HH-b', 'ox'): ForceFieldRequirement(
                'Oxidized_HisHisLigated_b-heme_RESP.lib',
                'Oxidized_HisHisLigated_b-heme.frcmod',
                'oxidized His-His ligated b-type heme',
                implemented=True
            ),
            ('HH-b', 'red'): ForceFieldRequirement(
                'Reduced_HisHisLigated_b-heme_RESP.lib',
                'Reduced_HisHisLigated_b-heme.frcmod',
                'reduced His-His ligated b-type heme',
                implemented=True
            ),
            # His-Met (implemented)
            ('HM-b', 'ox'): ForceFieldRequirement(
                'Oxidized_HisMetLigated_b-heme_RESP.lib',
                'Oxidized_HisMetLigated_b-heme.frcmod',
                'oxidized His-Met ligated b-type heme',
                implemented=True
            ),
            ('HM-b', 'red'): ForceFieldRequirement(
                'Reduced_HisMetLigated_b-heme_RESP.lib',
                'Reduced_HisMetLigated_b-heme.frcmod',
                'reduced His-Met ligated b-type heme',
                implemented=True
            ),
            # His-Cys
            ('HC-b', 'ox'): ForceFieldRequirement(
                'Oxidized_HisCysLigated_b-heme_RESP.lib',
                'Oxidized_HisCysLigated_b-heme.frcmod',
                'oxidized His-Cys ligated b-type heme'
            ),
            ('HC-b', 'red'): ForceFieldRequirement(
                'Reduced_HisCysLigated_b-heme_RESP.lib',
                'Reduced_HisCysLigated_b-heme.frcmod',
                'reduced His-Cys ligated b-type heme'
            ),
            # His-Tyr
            ('HY-b', 'ox'): ForceFieldRequirement(
                'Oxidized_HisTyrLigated_b-heme_RESP.lib',
                'Oxidized_HisTyrLigated_b-heme.frcmod',
                'oxidized His-Tyr ligated b-type heme'
            ),
            ('HY-b', 'red'): ForceFieldRequirement(
                'Reduced_HisTyrLigated_b-heme_RESP.lib',
                'Reduced_HisTyrLigated_b-heme.frcmod',
                'reduced His-Tyr ligated b-type heme'
            ),
            # His-Asp
            ('HD-b', 'ox'): ForceFieldRequirement(
                'Oxidized_HisAspLigated_b-heme_RESP.lib',
                'Oxidized_HisAspLigated_b-heme.frcmod',
                'oxidized His-Asp ligated b-type heme'
            ),
            ('HD-b', 'red'): ForceFieldRequirement(
                'Reduced_HisAspLigated_b-heme_RESP.lib',
                'Reduced_HisAspLigated_b-heme.frcmod',
                'reduced His-Asp ligated b-type heme'
            ),
            # His-Glu
            ('HE-b', 'ox'): ForceFieldRequirement(
                'Oxidized_HisGluLigated_b-heme_RESP.lib',
                'Oxidized_HisGluLigated_b-heme.frcmod',
                'oxidized His-Glu ligated b-type heme'
            ),
            ('HE-b', 'red'): ForceFieldRequirement(
                'Reduced_HisGluLigated_b-heme_RESP.lib',
                'Reduced_HisGluLigated_b-heme.frcmod',
                'reduced His-Glu ligated b-type heme'
            ),
            # His-Asn
            ('HN-b', 'ox'): ForceFieldRequirement(
                'Oxidized_HisAsnLigated_b-heme_RESP.lib',
                'Oxidized_HisAsnLigated_b-heme.frcmod',
                'oxidized His-Asn ligated b-type heme'
            ),
            ('HN-b', 'red'): ForceFieldRequirement(
                'Reduced_HisAsnLigated_b-heme_RESP.lib',
                'Reduced_HisAsnLigated_b-heme.frcmod',
                'reduced His-Asn ligated b-type heme'
            ),
            # His-Gln
            ('HQ-b', 'ox'): ForceFieldRequirement(
                'Oxidized_HisGlnLigated_b-heme_RESP.lib',
                'Oxidized_HisGlnLigated_b-heme.frcmod',
                'oxidized His-Gln ligated b-type heme'
            ),
            ('HQ-b', 'red'): ForceFieldRequirement(
                'Reduced_HisGlnLigated_b-heme_RESP.lib',
                'Reduced_HisGlnLigated_b-heme.frcmod',
                'reduced His-Gln ligated b-type heme'
            ),
            # His-Lys
            ('HK-b', 'ox'): ForceFieldRequirement(
                'Oxidized_HisLysLigated_b-heme_RESP.lib',
                'Oxidized_HisLysLigated_b-heme.frcmod',
                'oxidized His-Lys ligated b-type heme'
            ),
            ('HK-b', 'red'): ForceFieldRequirement(
                'Reduced_HisLysLigated_b-heme_RESP.lib',
                'Reduced_HisLysLigated_b-heme.frcmod',
                'reduced His-Lys ligated b-type heme'
            ),

            # c-type hemes
            # His-His (implemented)
            ('HH-c', 'ox'): ForceFieldRequirement(
                'Oxidized_HisHisLigated_c-heme_RESP.lib',
                'Oxidized_HisHisLigated_c-heme.frcmod',
                'oxidized His-His ligated c-type heme',
                implemented=True
            ),
            ('HH-c', 'red'): ForceFieldRequirement(
                'Reduced_HisHisLigated_c-heme_RESP.lib',
                'Reduced_HisHisLigated_c-heme.frcmod',
                'reduced His-His ligated c-type heme',
                implemented=True
            ),
            # His-Met (implemented)
            ('HM-c', 'ox'): ForceFieldRequirement(
                'Oxidized_HisMetLigated_c-heme_RESP.lib',
                'Oxidized_HisMetLigated_c-heme.frcmod',
                'oxidized His-Met ligated c-type heme',
                implemented=True
            ),
            ('HM-c', 'red'): ForceFieldRequirement(
                'Reduced_HisMetLigated_c-heme_RESP.lib',
                'Reduced_HisMetLigated_c-heme.frcmod',
                'reduced His-Met ligated c-type heme',
                implemented=True
            ),
            # His-Cys
            ('HC-c', 'ox'): ForceFieldRequirement(
                'Oxidized_HisCysLigated_c-heme_RESP.lib',
                'Oxidized_HisCysLigated_c-heme.frcmod',
                'oxidized His-Cys ligated c-type heme'
            ),
            ('HC-c', 'red'): ForceFieldRequirement(
                'Reduced_HisCysLigated_c-heme_RESP.lib',
                'Reduced_HisCysLigated_c-heme.frcmod',
                'reduced His-Cys ligated c-type heme'
            ),
            # His-Tyr
            ('HY-c', 'ox'): ForceFieldRequirement(
                'Oxidized_HisTyrLigated_c-heme_RESP.lib',
                'Oxidized_HisTyrLigated_c-heme.frcmod',
                'oxidized His-Tyr ligated c-type heme'
            ),
            ('HY-c', 'red'): ForceFieldRequirement(
                'Reduced_HisTyrLigated_c-heme_RESP.lib',
                'Reduced_HisTyrLigated_c-heme.frcmod',
                'reduced His-Tyr ligated c-type heme'
            ),
            # His-Asp
            ('HD-c', 'ox'): ForceFieldRequirement(
                'Oxidized_HisAspLigated_c-heme_RESP.lib',
                'Oxidized_HisAspLigated_c-heme.frcmod',
                'oxidized His-Asp ligated c-type heme'
            ),
            ('HD-c', 'red'): ForceFieldRequirement(
                'Reduced_HisAspLigated_c-heme_RESP.lib',
                'Reduced_HisAspLigated_c-heme.frcmod',
                'reduced His-Asp ligated c-type heme'
            ),
            # His-Glu
            ('HE-c', 'ox'): ForceFieldRequirement(
                'Oxidized_HisGluLigated_c-heme_RESP.lib',
                'Oxidized_HisGluLigated_c-heme.frcmod',
                'oxidized His-Glu ligated c-type heme'
            ),
            ('HE-c', 'red'): ForceFieldRequirement(
                'Reduced_HisGluLigated_c-heme_RESP.lib',
                'Reduced_HisGluLigated_c-heme.frcmod',
                'reduced His-Glu ligated c-type heme'
            ),
            # His-Asn
            ('HN-c', 'ox'): ForceFieldRequirement(
                'Oxidized_HisAsnLigated_c-heme_RESP.lib',
                'Oxidized_HisAsnLigated_c-heme.frcmod',
                'oxidized His-Asn ligated c-type heme'
            ),
            ('HN-c', 'red'): ForceFieldRequirement(
                'Reduced_HisAsnLigated_c-heme_RESP.lib',
                'Reduced_HisAsnLigated_c-heme.frcmod',
                'reduced His-Asn ligated c-type heme'
            ),
            # His-Gln
            ('HQ-c', 'ox'): ForceFieldRequirement(
                'Oxidized_HisGlnLigated_c-heme_RESP.lib',
                'Oxidized_HisGlnLigated_c-heme.frcmod',
                'oxidized His-Gln ligated c-type heme'
            ),
            ('HQ-c', 'red'): ForceFieldRequirement(
                'Reduced_HisGlnLigated_c-heme_RESP.lib',
                'Reduced_HisGlnLigated_c-heme.frcmod',
                'reduced His-Gln ligated c-type heme'
            ),
            # His-Lys
            ('HK-c', 'ox'): ForceFieldRequirement(
                'Oxidized_HisLysLigated_c-heme_RESP.lib',
                'Oxidized_HisLysLigated_c-heme.frcmod',
                'oxidized His-Lys ligated c-type heme'
            ),
            ('HK-c', 'red'): ForceFieldRequirement(
                'Reduced_HisLysLigated_c-heme_RESP.lib',
                'Reduced_HisLysLigated_c-heme.frcmod',
                'reduced His-Lys ligated c-type heme'
            ),
        }
    
    def validate_forcefield_files(self, heme_counts: Dict[str, Dict[str, int]]) -> Tuple[bool, List[str]]:
        """
        Check if required forcefield files exist for detected heme types.
        Only enforces validation for implemented parameter sets.
        """
        missing_files = []
        print("\n")
        print("=" * 60)
        print("Checking required forcefield files...")
        print("=" * 60)
        
        for heme_type, counts in heme_counts.items():
            for redox_state, count in counts.items():
                if count > 0:
                    key = (heme_type, redox_state)
                    if key not in self.requirements:
                        print(f"\nWarning: No parameter definitions yet for {heme_type} {redox_state}")
                        continue
                        
                    req = self.requirements[key]
                    lib_path = self.forcefield_dir / req.lib_file
                    frcmod_path = self.forcefield_dir / req.frcmod_file
                    
                    if req.implemented:
                        status = "implemented"
                    else:
                        status = "not yet implemented"
                        
                    # Check files
                    if not lib_path.exists():
                        missing_files.append(req.lib_file)
                        print(f"   ! Missing the .lib file for {req.description} ({status})")
                    else:
                        print(f"   √ Found the .lib file for {req.description}")
                        
                    if not frcmod_path.exists():
                        missing_files.append(req.frcmod_file)
                        print(f"   ! Missing the .frcmod file for {req.description} ({status})")
                    else:
                        print(f"   √ Found the .frcmod file for {req.description}")
        
        # Only return failure for missing files of implemented parameter sets
        implemented_missing = [
            f for f in missing_files 
            if any(key for key, req in self.requirements.items() 
                  if req.implemented and (f == req.lib_file or f == req.frcmod_file))
        ]
        
        if implemented_missing:
            print("\nError: Missing required files for implemented parameter sets:")
            for f in implemented_missing:
                print(f"  - {f}")
        
        if missing_files != implemented_missing:
            print("\nNote: Some parameter sets are not yet implemented.")
            print("These will be available in future updates.")
        
        return len(implemented_missing) == 0, missing_files

    def get_parameter_status(self) -> Dict[str, Dict[str, Dict[str, bool]]]:
        """
        Returns dictionary showing development and availability status of all parameter sets.
        
        Returns:
            Dict with structure:
            {
                'b-type': {
                    'His-His': {'implemented': True, 'lib': True, 'frcmod': True},
                    'His-Met': {'implemented': True, 'lib': True, 'frcmod': True},
                    ...
                },
                'c-type': {
                    'His-His': {'implemented': True, 'lib': True, 'frcmod': True},
                    ...
                }
            }
        """
        status = {'b-type': {}, 'c-type': {}}
        
        for (heme_type, redox), req in self.requirements.items():
            type_key = 'c-type' if heme_type.endswith('-c') else 'b-type'
            ligand_key = heme_type[:-2]  # Remove -b or -c suffix
            
            if ligand_key not in status[type_key]:
                status[type_key][ligand_key] = {
                    'implemented': req.implemented,
                    'lib': False,
                    'frcmod': False
                }
            
            # Check file existence
            lib_exists = (self.forcefield_dir / req.lib_file).exists()
            frcmod_exists = (self.forcefield_dir / req.frcmod_file).exists()
            
            # Update status
            status[type_key][ligand_key]['lib'] |= lib_exists
            status[type_key][ligand_key]['frcmod'] |= frcmod_exists
            
        return status

    def get_alternative_charges(self, heme_type: str) -> Optional[Tuple[str, str]]:
        """
        Get alternative charge sets if available for a heme type.
        Returns (lib_file, description) if alternatives exist.
        """
        alternatives = {
            'HH-c': {
                'Henriques': (
                    'Henriques_Oxidized_HisHisLigated_c-heme_RESP.lib',
                    'Henriques et al. charge set for His-His c-type hemes',
                    """
  Henriques, J.; Costa, P. J.; Calhorda, M. J.; 
  Machuqueiro, M. Charge Parametrization of the
  DvH-c3 Heme Group: Validation Using Constant-(pH,E) 
  Molecular Dynamics Simulations. J. Phys. Chem. B 
  2013, 117 (1), 70–82."""
                )
            }
            # Add other alternative parameter sets as they become available
        }
        
        return alternatives.get(heme_type, {}).get('Henriques')

    def print_development_status(self):
        """Print a summary of parameter development status."""
        status = self.get_parameter_status()
        
        print("\nParameter Development Status:")
        print("-" * 80)
        print(f"{'Heme Type':<10} {'Ligation':<12} {'Implemented':<12} {'Files Ready':<12}")
        print("-" * 80)
        
        for heme_type in ['b-type', 'c-type']:
            for ligation, details in sorted(status[heme_type].items()):
                implemented = "✓" if details['implemented'] else "×"
                files_ready = "✓" if details['lib'] and details['frcmod'] else "×"
                print(f"{heme_type:<10} {ligation:<12} {implemented:<12} {files_ready:<12}")
        
        print("-" * 80)
        print("✓ = Ready    × = In Development")

def process_residue_indexing(pdb: str,
                           disulf_list: str,
                           sel_asp_ids: str,
                           sel_glu_ids: str,
                           sel_his_ids: str,
                           sel_lys_ids: str,
                           sel_tyr_ids: str,
                           sel_prn_ids: str,
                           input_dict: Dict,
                           launch_dir: Path) -> None:
    """
    Main function for processing residue indexing and renaming.

    Args:
        pdb: Base name of the PDB file
        disulf_list: List of disulfide-bonded residues
        sel_asp_ids: Selected Asp residue IDs
        sel_glu_ids: Selected Glu residue IDs
        sel_his_ids: Selected His residue IDs
        sel_lys_ids: Selected Lys residue IDs
        sel_tyr_ids: Selected Tyr residue IDs
        sel_prn_ids: Selected PRN residue IDs
        input_dict: Configuration dictionary
        launch_dir: Directory where the script is launched
    """
    # Create InteractionManager
    interaction_manager = InteractionManager(
        launch_dir=launch_dir,
        input_dict=input_dict
    )

    # Add ForceFieldDir if not present
    if 'ForceFieldDir' not in input_dict:
        input_dict['ForceFieldDir'] = Path(__file__).parent.parent.parent / 'data' / 'forcefield'
    
    launch_dir = Path(launch_dir)
    indexing_file = Path.cwd() / "ResIndexing.txt"

    if not indexing_file.exists():
        raise FileNotFoundError("""
ResIndexing.txt was not found!
I, unfortunately, do not know how to proceed without this file.""")

    print("\n")
    print("=" * 60)
    print("""Structure Processing:

Your ResIndexing file will now be used to re-label the residues
according to the available parameterizations of b- and c-type,
His-[ligand] ligated hemes.

Also, if you specified any disulfide residues earlier, the
participating Cys residues will be re-labeled as CYX, in
accordance with AMBER conventions.""")

    # Get user confirmation
    choice = interaction_manager.yes_no_prompt(
        "ProcessPDBChoice",
        "Would you like to proceed with the structure processing?"
    )

    if not choice:
        sys.exit("""
I'm sorry but I don't know how to proceed without processing
the structure. Hopefully this program was helpful up to this
point! Please don't hesitate to re-run this program when ready.""")

    # Determine the PDB file path
    possible_paths = [
        Path(f"{pdb}_renumd.pdb"),  # Renumbered file
        Path(f"{pdb}.pdb")          # Original file
    ]

    pdb_path = None
    for path in possible_paths:
        full_path = Path.cwd() / path
        if full_path.is_file():
            pdb_path = full_path
            break

    if not pdb_path:
        raise FileNotFoundError(f"Could not find PDB file for {pdb}")

    # Initialize PDB processor
    processor = PDBProcessor(pdb_path)

    # Process disulfides if any
    if disulf_list:
        processor.process_disulfides(disulf_list)

    # Process titratable residues
    processor.process_titratable_residues(sel_asp_ids, sel_glu_ids, sel_his_ids)

    # First pass: Calculate heme shifts and move hemes
    print("\nStarting first pass: calculating and applying heme shifts...")
    processor.heme_shifts = processor.calculate_heme_shifts(indexing_file)
    processor.shift_hemes()

    # Process heme environments
    print("\n")
    print("=" * 60)
    print("Processing heme environments...")
    print("=" * 60)
    with open(indexing_file) as f:
        for idx, line in enumerate(f, 1):
            # Skip empty or comment lines
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            # Split the line
            parts = line.split()

            # Validate line has enough parts
            if len(parts) not in [5, 7]:
                print(f"Warning: Skipping malformed line in ResIndexing.txt: {line}")
                continue

            try:
                # Determine if it's a b-type or c-type heme
                is_c_type = len(parts) == 7

                if is_c_type:
                    # C-type heme: 7 fields
                    # CysB CysC HisP Distal OriginalHemeID c HX
                    cys_b = int(parts[0])
                    cys_c = int(parts[1])
                    his_p = int(parts[2])
                    distal_ligand = int(parts[3])
                    original_heme_id = int(parts[4])
                    heme_type = parts[5]
                    ligand_type = parts[6]
                else:
                    # B-type heme: 5 fields
                    # HisP Distal OriginalHemeID b HX
                    his_p = int(parts[0])
                    distal_ligand = int(parts[1])
                    original_heme_id = int(parts[2])
                    heme_type = parts[3]
                    ligand_type = parts[4]
                    cys_b = None
                    cys_c = None

                # Map ligand type
                ligand_map = {
                    'HH': 'HIS', 'HM': 'MET', 'HC': 'CYS',
                    'HY': 'TYR', 'HD': 'ASP', 'HE': 'GLU',
                    'HN': 'ASN', 'HQ': 'GLN', 'HK': 'LYS'
                }
                distal_type = ligand_map.get(ligand_type, 'HIS')

                # Determine redox state
                redox_states = input_dict.get("RedoxState", [])
                if idx-1 < len(redox_states):
                    redox = redox_states[idx-1]
                else:
                    print("\n")
                    print("-" * 60)
                    print(f"Heme-{original_heme_id}: {heme_type}-type with {ligand_type} ligation")
                    print("-" * 60)
                    redox = interaction_manager.prompt(
                        f"RedoxState_{idx}",
                        "Model as oxidized or reduced (O/R)?",
                        choices=['O', 'R']
                    )
                    # Store the selection in input_dict
                    if "RedoxState" not in input_dict:
                        input_dict["RedoxState"] = []
                    input_dict["RedoxState"].append(redox)
                redox_state = 'ox' if redox.lower() in ('o', 'ox') else 'red'

                # Create environment object
                env = HemeEnvironment(
                    heme_id=original_heme_id,
                    shifted_id=processor.heme_shifts[original_heme_id],
                    new_id=processor.heme_shifts[original_heme_id],
                    his_p=his_p,
                    distal_ligand=distal_ligand,
                    distal_type=distal_type,
                    is_c_type=is_c_type,
                    cys_b=cys_b,
                    cys_c=cys_c,
                    redox_state=redox_state
                )

                # Process the environment
                processor.process_heme_environment(env)

            except Exception as e:
                print(f"Error processing line {line}: {e}")
                continue

    # Save final structure with proper ordering
    output_path = launch_dir / "processed.pdb"
    processor.save_structure(output_path)
    print("\nStructure processing completed.")

    # Get the residue mapping before validation
    id_mapping = processor.id_mapping
    print(f"Obtained ID mapping with {len(id_mapping)} entries")

    # Now validate
    validation_passed = processor.validate_structure(indexing_file)

    if validation_passed:
        print(f"\nFinal structure saved to: {output_path}")

        # Update ResIndexing.txt with new heme IDs...
        print("\n")
        print("=" * 60)
        print("Updating ResIndexing.txt with new residue IDs...")
        print("=" * 60)
        with open(indexing_file, 'r') as f:
            lines = f.readlines()

        updated_lines = []
        for line in lines:
            parts = line.strip().split()
            if not parts:  # Keep empty lines
                updated_lines.append(line)
                continue

            try:
                if len(parts) == 5:  # b-type: HisP Distal Heme b HX
                    old_heme_id = int(parts[2])
                    if old_heme_id in processor.id_mapping:
                        new_heme_id = processor.id_mapping[old_heme_id]
                        updated_line = f"{parts[0]} {parts[1]} {old_heme_id} {new_heme_id} {parts[3]} {parts[4]}\n"
                        print(f"  Updated b-type entry: heme {old_heme_id} → {new_heme_id}")
                    else:
                        print(f"  Warning: No mapping found for b-type heme {old_heme_id}")
                        updated_line = line

                elif len(parts) == 7:  # c-type: CysB CysC HisP Distal Heme c HX
                    old_heme_id = int(parts[4])
                    if old_heme_id in processor.id_mapping:
                        new_heme_id = processor.id_mapping[old_heme_id]
                        updated_line = f"{parts[0]} {parts[1]} {parts[2]} {parts[3]} {old_heme_id} {new_heme_id} {parts[5]} {parts[6]}\n"
                        print(f"  Updated c-type entry: heme {old_heme_id} → {new_heme_id}")
                    else:
                        print(f"  Warning: No mapping found for c-type heme {old_heme_id}")
                        updated_line = line
                else:
                    updated_line = line

                updated_lines.append(updated_line)

            except (ValueError, IndexError) as e:
                print(f"Warning: Error processing line: {line.strip()}")
                updated_lines.append(line)
                continue

        # Write updated content back to file
        with open(indexing_file, 'w') as f:
            f.writelines(updated_lines)

        print(f"\nResIndexing.txt updated successfully")

    # Get structure info
    structure_info = processor.count_heme_types()  # Original format

    # Transform structure_info into format needed by generate_tleap
    tleap_structure_info = {}
    for heme_type, states in structure_info.items():
        for redox_state, count in states.items():
            if count > 0:  # Only include types that exist
                tleap_structure_info[(heme_type, redox_state)] = count

    # Collect indexing data from updated ResIndexing.txt
    indexing_data = []

    # Add disulfide information if present
    if disulf_list and disulf_list.strip():
        # Convert space-separated string to pairs of integers
        nums = [int(x) for x in disulf_list.strip().split()]
        disulf_pairs = list(zip(nums[::2], nums[1::2]))  # Create pairs
        indexing_data.append({
            'type': 'disulfide',
            'pairs': disulf_pairs
    })

    with open(indexing_file) as f:
        for line in f:
            parts = line.strip().split()
            if not parts:
                continue

            try:
                # Only need to handle updated format:
                # C-type: CysB CysC HisP Distal OldHemeID NewHemeID c HX (8 fields)
                # B-type: HisP Distal OldHemeID NewHemeID b HX (6 fields)
                is_c_type = len(parts) == 8

                if is_c_type:
                    data = {
                        'heme_id': int(parts[4]),      # OldHemeID
                        'shifted_id': int(parts[5]),   # NewHemeID
                        'new_id': int(parts[5]),       # NewHemeID
                        'his_p': int(parts[2]),        # HisP
                        'distal_ligand': int(parts[3]), # Distal
                        'distal_type': parts[7],       # HX
                        'is_c_type': True,
                        'cys_b': int(parts[0]),        # CysB
                        'cys_c': int(parts[1])         # CysC
                    }
                else:
                    data = {
                        'heme_id': int(parts[2]),      # OldHemeID
                        'shifted_id': int(parts[3]),   # NewHemeID
                        'new_id': int(parts[3]),       # NewHemeID
                        'his_p': int(parts[0]),        # HisP
                        'distal_ligand': int(parts[1]), # Distal
                        'distal_type': parts[5],       # HX
                        'is_c_type': False
                    }

                # Get redox state
                heme_count = sum(1 for entry in indexing_data if isinstance(entry, dict) and 'type' not in entry)
                idx = heme_count + 1
                redox_states = input_dict.get("RedoxState", [])
                if idx-1 < len(redox_states):
                    redox = redox_states[idx-1]
                else:
                    redox = interaction_manager.prompt(
                        f"RedoxState_{idx}",
                        "Model as oxidized or reduced (O/R)?",
                        choices=['O', 'R']
                    )
                data['redox_state'] = 'ox' if redox.lower() in ('o', 'ox') else 'red'

                indexing_data.append(data)
            except Exception as e:
                print(f"Warning: Error processing line: {line.strip()}: {e}")
                continue

    # Validate forcefield files
    validator = ForceFieldValidator(input_dict['ForceFieldDir'])
    heme_counts = processor.count_heme_types()

    success, missing_files = validator.validate_forcefield_files(heme_counts)

    if not success:
        print("\nWarning: Missing required forcefield files:")
        for file in missing_files:
            print(f"  - {file}")
        print("\nPlease ensure all required forcefield files are present.")
        raise FileNotFoundError("Missing required forcefield files")

    # Check for alternative charge sets
    for heme_type in heme_counts:
        if heme_counts[heme_type]['ox'] > 0:  # If we have oxidized hemes
            alt_charges = validator.get_alternative_charges(heme_type)
            if alt_charges and "FFchoice" not in input_dict:
                lib_file, description, reference = alt_charges
                print(f"\nAlternative charges available for {heme_type}:")
                print(f"Description: {description}")
                print(f"Reference: {reference}")
                choice = interaction_manager.yes_no_prompt(
                    "FFchoice",
                    f"Would you like to use these alternative charges?"
                )

    # Return the transformed structure info
    return tleap_structure_info, indexing_data

