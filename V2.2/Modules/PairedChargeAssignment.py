################################################################################
# Refactored PairedChargeAssignment Module
# This module handles the generation of topology files and coordinates for pairs 
# of heme groups in different redox states.
################################################################################

from dataclasses import dataclass
from enum import Enum
from typing import List, Dict, Tuple, Optional, Set
import os
import sys
import itertools
import subprocess
from pathlib import Path

class HemeType(Enum):
    B = "b"
    C = "c"

class LigationType(Enum):
    HIS_HIS = "HH"
    HIS_MET = "HM"

class RedoxState(Enum):
    OXIDIZED = "o"
    REDUCED = "r"

@dataclass
class HemeDefinition:
    """Represents a single heme definition with its properties"""
    heme_id: int
    heme_type: HemeType
    ligation_type: LigationType
    cys_b: Optional[int] = None  # Only for c-type hemes
    cys_c: Optional[int] = None  # Only for c-type hemes
    ligand_proximal: int = 0
    ligand_distal: int = 0

    def is_c_type(self) -> bool:
        return self.heme_type == HemeType.C

    def is_his_his(self) -> bool:
        return self.ligation_type == LigationType.HIS_HIS

class PDBProcessor:
    """Handles PDB file processing and TER record insertion"""
    
    @staticmethod
    def insert_ter_records(pdb_file: str, output_file: str):
        """Insert TER records at appropriate positions in PDB file"""
        try:
            with open(pdb_file, 'r') as file:
                lines = file.readlines()

            modified_lines = []
            
            i = 0
            while i < len(lines):
                line = lines[i]
                
                # Skip non-ATOM/HETATM lines but preserve them
                if not line.startswith(('ATOM', 'HETATM')):
                    if not line.startswith('TER'):  # Don't duplicate TER records
                        modified_lines.append(line)
                    i += 1
                    continue

                # Case 1: After OXT atom
                if line[12:16].strip() == "OXT":
                    modified_lines.append(line)
                    if i + 1 < len(lines) and not lines[i+1].startswith('TER'):
                        chain = line[21]
                        resnum = int(line[22:26])
                        ter_record = f"TER   {str(len(modified_lines)+1).rjust(5)}      {resnum:4}    {chain}\n"
                        modified_lines.append(ter_record)
                    i += 1
                    continue

                # Case 2: Before N-terminal residue
                if line[12:16].strip() == "N":
                    # Check for H1, H2, H3 pattern
                    if (i + 3 < len(lines) and
                        all(lines[i+j].startswith('ATOM') and 
                            lines[i+j][12:16].strip() == f"H{j}" 
                            for j in range(1, 4))):
                        if not (modified_lines and modified_lines[-1].startswith('TER')):
                            chain = line[21]
                            resnum = int(line[22:26])
                            ter_record = f"TER   {str(len(modified_lines)+1).rjust(5)}      {resnum:4}    {chain}\n"
                            modified_lines.append(ter_record)
                    # Check for H2, H3 pattern
                    elif (i + 2 < len(lines) and
                          lines[i+1].startswith('ATOM') and lines[i+1][12:16].strip() == "H2" and
                          lines[i+2].startswith('ATOM') and lines[i+2][12:16].strip() == "H3"):
                        if not (modified_lines and modified_lines[-1].startswith('TER')):
                            chain = line[21]
                            resnum = int(line[22:26])
                            ter_record = f"TER   {str(len(modified_lines)+1).rjust(5)}      {resnum:4}    {chain}\n"
                            modified_lines.append(ter_record)

                modified_lines.append(line)
                i += 1

            # Add final TER if the last line was an ATOM/HETATM record
            if modified_lines and modified_lines[-1].startswith(("ATOM", "HETATM")):
                last_line = modified_lines[-1]
                chain = last_line[21]
                resnum = int(last_line[22:26])
                ter_record = f"TER   {str(len(modified_lines)+1).rjust(5)}      {resnum:4}    {chain}\n"
                modified_lines.append(ter_record)

            # Write the modified file
            with open(output_file, 'w') as file:
                file.writelines(modified_lines)

        except Exception as e:
            print(f"Error processing PDB file: {str(e)}")
            raise

class ForceFieldParameters:
    """Manages force field parameters and their loading"""
    def __init__(self, force_field_dir: str, ff_choice: str):
        self.dir = force_field_dir
        self.ff_choice = ff_choice.lower()

    def get_lib_name(self, heme_type: HemeType, ligation: LigationType, 
                     redox: RedoxState) -> str:
        """Get library file name based on heme properties"""
        prefix = "Henriques_" if (self.ff_choice == "henriques" and 
                                 heme_type == HemeType.C and 
                                 ligation == LigationType.HIS_HIS) else ""
        redox_state = "Oxidized" if redox == RedoxState.OXIDIZED else "Reduced"
        ligation_text = "HisHisLigated" if ligation == LigationType.HIS_HIS else "HisMetLigated"
        heme_text = f"{heme_type.value}-heme"
        
        return f"{prefix}{redox_state}_{ligation_text}_{heme_text}_RESP.lib"

    def get_frcmod_name(self, heme_type: HemeType, ligation: LigationType, 
                        redox: RedoxState) -> str:
        """Get force field modification file name"""
        redox_state = "Oxidized" if redox == RedoxState.OXIDIZED else "Reduced"
        ligation_text = "HisHisLigated" if ligation == LigationType.HIS_HIS else "HisMetLigated"
        heme_text = f"{heme_type.value}-heme"
        
        return f"{redox_state}_{ligation_text}_{heme_text}.frcmod"

class VMDScriptGenerator:
    """Generates VMD scripts for redox state modifications"""
    def __init__(self, output_file: str):
        print(f"Opening {output_file} for writing...")
        self.file = open(output_file, 'w')
        self._write_header()
        
    def _write_header(self):
        print("\nmol new RefState.pdb", file=self.file)
        
    def add_heme_modification(self, heme: HemeDefinition, redox: RedoxState):
        """Add commands for modifying a specific heme"""
        state = "Oxidized" if redox == RedoxState.OXIDIZED else "Reduced"
        residue_prefixes = self._get_residue_prefixes(heme, redox)
        selection_patterns = self._get_selection_patterns(heme.heme_type, heme.ligation_type)
        
        print(f"\n #-------------------------------------------------------------------------", file=self.file)
        print(f" #{state} Heme-{heme.heme_id}: {heme.ligation_type.value} ligated {heme.heme_type.value}-type heme", file=self.file)
        
        self._write_atom_selections(heme, selection_patterns)
        self._write_residue_modifications(heme, residue_prefixes, state)

    def _get_selection_patterns(self, heme_type: HemeType, ligation: LigationType) -> Dict[str, str]:
        """Get patterns that match both oxidized and reduced states for selections"""
        patterns = {
            (HemeType.C, LigationType.HIS_HIS): {
                'his_p': 'PHO PHR',  # Proximal His (both oxidized and reduced)
                'his_d': 'DHO DHR',  # Distal His
                'hem': 'HCO HCR'     # c-type heme
            },
            (HemeType.C, LigationType.HIS_MET): {
                'his_p': 'PMO PMR',  # Proximal His
                'met_d': 'DMO DMR',  # Distal Met
                'hem': 'MCO MCR'     # c-type heme
            },
            (HemeType.B, LigationType.HIS_HIS): {
                'his_p': 'FHO FHR',  # Proximal His
                'his_d': 'RHO RHR',  # Distal His
                'hem': 'HBO HBR'     # b-type heme
            },
            (HemeType.B, LigationType.HIS_MET): {
                'his_p': 'FMO FMR',  # Proximal His
                'met_d': 'RMO RMR',  # Distal Met
                'hem': 'MBO MBR'     # b-type heme
            }
        }
        return patterns[(heme_type, ligation)]

    def _get_residue_prefixes(self, heme: HemeDefinition, redox: RedoxState) -> Dict[str, str]:
        """Get the target residue names based on desired redox state"""
        is_oxidized = redox == RedoxState.OXIDIZED
        
        prefix_mapping = {
            (HemeType.C, LigationType.HIS_HIS): {
                'his_p': 'PHO' if is_oxidized else 'PHR',
                'his_d': 'DHO' if is_oxidized else 'DHR',
                'hem': 'HCO' if is_oxidized else 'HCR'
            },
            (HemeType.C, LigationType.HIS_MET): {
                'his_p': 'PMO' if is_oxidized else 'PMR',
                'met_d': 'DMO' if is_oxidized else 'DMR',
                'hem': 'MCO' if is_oxidized else 'MCR'
            },
            (HemeType.B, LigationType.HIS_HIS): {
                'his_p': 'FHO' if is_oxidized else 'FHR',
                'his_d': 'RHO' if is_oxidized else 'RHR',
                'hem': 'HBO' if is_oxidized else 'HBR'
            },
            (HemeType.B, LigationType.HIS_MET): {
                'his_p': 'FMO' if is_oxidized else 'FMR',
                'met_d': 'RMO' if is_oxidized else 'RMR',
                'hem': 'MBO' if is_oxidized else 'MBR'
            }
        }
        return prefix_mapping[(heme.heme_type, heme.ligation_type)]

    def _write_atom_selections(self, heme: HemeDefinition, selection_patterns: Dict[str, str]):
        """Write VMD atom selection commands using patterns that match both states"""
        is_his_met = heme.ligation_type == LigationType.HIS_MET
        ligand_key = 'met_d' if is_his_met else 'his_d'
        
        print(f"""
 #Define atom groups
   set HISp  [atomselect top "resname {selection_patterns['his_p']} and resid {heme.ligand_proximal}"]
   set {'METd' if is_his_met else 'HISd'}  [atomselect top "resname {selection_patterns[ligand_key]} and resid {heme.ligand_distal}"]
   set HEM   [atomselect top "resname {selection_patterns['hem']} and resid {heme.heme_id}"]""", file=self.file)

    def _write_residue_modifications(self, heme: HemeDefinition, prefixes: Dict[str, str], state: str):
        """Write residue name modification commands"""
        is_his_met = heme.ligation_type == LigationType.HIS_MET
        ligand_key = 'met_d' if is_his_met else 'his_d'
        ligand_type = "Met" if is_his_met else "His"
        
        print(f"""
 #Change residue names for {state.lower()} state:
   #The proximal His residue 
     $HISp  set resname {prefixes['his_p']}; #Proximal His for {state.lower()} {ligand_type}-{ligand_type} ligated heme.
   #The distal {ligand_type} residue
     ${'METd' if is_his_met else 'HISd'}  set resname {prefixes[ligand_key]}; #Distal   {ligand_type} for {state.lower()} {ligand_type}-{ligand_type} ligated heme.
   #The heme group
     $HEM   set resname {prefixes['hem']}; #{heme.heme_type.value}-type {ligand_type}-{ligand_type} ligated {state.lower()} heme""", file=self.file)

    def add_final_commands(self, output_prefix: str):
        """Write commands to save the PDB file"""
        print(f"""
 set sel [atomselect top "all and not resname WAT 'Na+' 'Cl-'"]
 $sel writepdb {output_prefix}.pdb
 #-------------------------------------------------------------------------""", file=self.file)

    def finalize(self):
        """Write final exit command and close file"""
        print("\nexit", file=self.file)
        self.file.close()

class TLeapScriptGenerator:
    """Generates TLeap input scripts for topology generation"""
    def __init__(self, ff_params: ForceFieldParameters):
        self.ff_params = ff_params
        self.defined_bonds: Set[Tuple[int, int]] = set()
        # Track what we've written to avoid duplicates
        self.written_atom_types: Set[Tuple[HemeType, LigationType, RedoxState]] = set()
        self.written_ff_params: Set[Tuple[HemeType, LigationType, RedoxState]] = set()
        
    def generate_script(self, heme_pair: Tuple[HemeDefinition, HemeDefinition], output_file: str):
        """Generate the complete TLeap script"""
        try:
            print(f"Opening {output_file} for writing...")
            with open(output_file, 'w') as f:
                print("Writing header...")
                self._write_header(f)
                print("Writing atom types...")
                self._write_atom_types(f, heme_pair)
                print("Writing forcefield loads...")
                self._write_forcefield_loads(f, heme_pair)
                print("Writing structure loads...")
                self._write_structure_loads(f, heme_pair[0].heme_id, heme_pair[1].heme_id)
                print("Writing disulfide bonds...")
                self._write_disulfide_bonds(f)
                print("Writing heme bonds...")
                self._write_heme_bonds(f, heme_pair)
                print("Writing footer...")
                self._write_footer(f, heme_pair[0].heme_id, heme_pair[1].heme_id)
        except Exception as e:
            print(f"Failed to generate TLeap script: {str(e)}")
            raise

    def _write_header(self, f):
        """Write the initial TLeap commands"""
        print("""
# Load parameters
 source leaprc.constph
 source leaprc.conste
 source leaprc.gaff
 source leaprc.water.tip3p

 addAtomTypes {""", file=f)

    def _write_atom_types(self, f, heme_pair: Tuple[HemeDefinition, HemeDefinition]):
        """Write atom type definitions based on heme types and their states"""
        try:
            print("Reading all heme types from SelResIndexing.txt...")
            # First, get the unique heme types from SelResIndexing.txt
            heme_types = set()
            with open("SelResIndexing.txt") as fp:
                for line in fp:
                    parts = line.strip().split()
                    heme_type = HemeType(parts[-2])
                    ligation_type = LigationType(parts[-1])
                    heme_types.add((heme_type, ligation_type))
        
            print(f"Found {len(heme_types)} unique heme type combinations")
        
            # Write atom types for all heme types in the structure
            for heme_type, ligation_type in heme_types:
                print(f"Processing {heme_type.value}-type {ligation_type.value} heme")
                for redox_state in [RedoxState.OXIDIZED, RedoxState.REDUCED]:
                    key = (heme_type, ligation_type, redox_state)
                    if key not in self.written_atom_types:
                        print(f"  Writing {redox_state.value} state atom types")
                        if heme_type == HemeType.C:
                            if ligation_type == LigationType.HIS_HIS:
                                self._write_c_type_his_his_atoms(f, redox_state)
                            else:
                                self._write_c_type_his_met_atoms(f, redox_state)
                        else:  # b-type
                            if ligation_type == LigationType.HIS_HIS:
                                self._write_b_type_his_his_atoms(f, redox_state)
                            else:
                                self._write_b_type_his_met_atoms(f, redox_state)
                        self.written_atom_types.add(key)
                    else:
                        print(f"  Skipping {redox_state.value} state (already written)")
        
            print("\n }", file=f)
        except Exception as e:
            print(f"Error in _write_atom_types: {str(e)}")
            raise

    def _write_b_type_his_his_atoms(self, f, redox_state: RedoxState):
        """Write atom types for b-type His-His ligated heme"""
        if redox_state == RedoxState.OXIDIZED:
            print("""
        { "M1"  "Fe" "sp3" } #M1&Y1-Y6:
        { "Y1"  "N" "sp3" }  #Oxidized
        { "Y2"  "N" "sp3" }  #His-His
        { "Y3"  "N" "sp3" }  #Ligated
        { "Y4"  "N" "sp3" }  #b-Heme
        { "Y5"  "N" "sp3" }
        { "Y6"  "N" "sp3" }""", end=" ", file=f)
        else:  # REDUCED
            print("""
        { "M2"  "Fe" "sp3" } #M2&Z1-Z6:
        { "Z1"  "N" "sp3" }  #Reduced
        { "Z2"  "N" "sp3" }  #His-His
        { "Z3"  "N" "sp3" }  #Ligated
        { "Z4"  "N" "sp3" }  #b-Heme
        { "Z5"  "N" "sp3" }
        { "Z6"  "N" "sp3" }""", end=" ", file=f)

    def _write_b_type_his_met_atoms(self, f, redox_state: RedoxState):
        """Write atom types for b-type His-Met ligated heme"""
        if redox_state == RedoxState.OXIDIZED:
            print("""
        { "M3"  "Fe" "sp3" } #M3&W1-W6:
        { "W1"  "S" "sp3" }  #Oxidized
        { "W2"  "N" "sp3" }  #His-Met
        { "W3"  "N" "sp3" }  #Ligated
        { "W4"  "N" "sp3" }  #b-Heme
        { "W5"  "N" "sp3" }
        { "W6"  "N" "sp3" }""", end=" ", file=f)
        else:  # REDUCED
            print("""
        { "M4"  "Fe" "sp3" } #M4&X1-X6:
        { "X1"  "S" "sp3" }  #Reduced
        { "X2"  "N" "sp3" }  #His-Met
        { "X3"  "N" "sp3" }  #Ligated
        { "X4"  "N" "sp3" }  #b-Heme
        { "X5"  "N" "sp3" }
        { "X6"  "N" "sp3" }""", end=" ", file=f)

    def _write_c_type_his_his_atoms(self, f, redox_state: RedoxState):
        """Write atom types for c-type His-His ligated heme"""
        if redox_state == RedoxState.OXIDIZED:
            print("""
        { "M7"  "Fe" "sp3" } #M7&S1-S6:
        { "S1"  "N" "sp3" }  #Oxidized
        { "S2"  "N" "sp3" }  #His-His
        { "S3"  "N" "sp3" }  #Ligated
        { "S4"  "N" "sp3" }  #c-Heme
        { "S5"  "N" "sp3" }
        { "S6"  "N" "sp3" }""", end=" ", file=f)
        else:  # REDUCED
            print("""
        { "M8"  "Fe" "sp3" } #M8&T1-T6:
        { "T1"  "N" "sp3" }  #Reduced
        { "T2"  "N" "sp3" }  #His-His
        { "T3"  "N" "sp3" }  #Ligated
        { "T4"  "N" "sp3" }  #c-Heme
        { "T5"  "N" "sp3" }
        { "T6"  "N" "sp3" }""", end=" ", file=f)

    def _write_c_type_his_met_atoms(self, f, redox_state: RedoxState):
        """Write atom types for c-type His-Met ligated heme"""
        if redox_state == RedoxState.OXIDIZED:
            print("""
        { "M5"  "Fe" "sp3" } #M5&U1-U6:
        { "U1"  "N" "sp3" }  #Oxidized
        { "U2"  "S" "sp3" }  #His-Met
        { "U3"  "N" "sp3" }  #Ligated
        { "U4"  "N" "sp3" }  #c-Heme
        { "U5"  "N" "sp3" }
        { "U6"  "N" "sp3" }""", end=" ", file=f)
        else:  # REDUCED
            print("""
        { "M6"  "Fe" "sp3" } #M6&V1-V6:
        { "V1"  "N" "sp3" }  #Reduced
        { "V2"  "S" "sp3" }  #His-Met
        { "V3"  "N" "sp3" }  #Ligated
        { "V4"  "N" "sp3" }  #c-Heme
        { "V5"  "N" "sp3" }
        { "V6"  "N" "sp3" }""", end=" ", file=f)

    def _write_forcefield_loads(self, f, heme_pair: Tuple[HemeDefinition, HemeDefinition]):
        """Write force field loading commands for all heme types"""
        # Get unique heme types from SelResIndexing.txt
        heme_types = set()
        with open("SelResIndexing.txt") as fp:
            for line in fp:
                parts = line.strip().split()
                heme_type = HemeType(parts[-2])
                ligation_type = LigationType(parts[-1])
                heme_types.add((heme_type, ligation_type))

        # Load force field parameters for all heme types
        for heme_type, ligation_type in heme_types:
            for redox_state in [RedoxState.OXIDIZED, RedoxState.REDUCED]:
                key = (heme_type, ligation_type, redox_state)
                if key not in self.written_ff_params:
                    lib_name = self.ff_params.get_lib_name(heme_type, ligation_type, redox_state)
                    frcmod_name = self.ff_params.get_frcmod_name(heme_type, ligation_type, redox_state)
                
                    print(f"""
 loadamberparams {self.ff_params.dir}/{frcmod_name}
 loadoff {self.ff_params.dir}/{lib_name}""", file=f)
                    self.written_ff_params.add(key)

    def _write_structure_loads(self, f, heme1_id: int, heme2_id: int):
        """Write structure loading commands"""
        print(f"""
# Load PDB files
 oo = loadpdb o{heme1_id}-o{heme2_id}.pdb
 or = loadpdb o{heme1_id}-r{heme2_id}.pdb
 ro = loadpdb r{heme1_id}-o{heme2_id}.pdb
 rr = loadpdb r{heme1_id}-r{heme2_id}.pdb""", file=f)

    def _write_disulfide_bonds(self, f):
        """Write disulfide bond definitions if present"""
        if not os.path.isfile("DisulfideDefinitions.txt"):
            return
            
        with open("DisulfideDefinitions.txt") as dsl:
            disulfide_pairs = [list(map(int, line.split())) for line in dsl]
            
        if disulfide_pairs:
            print("\n# Define Disulfide linkages", file=f)
            for pair in disulfide_pairs:
                for prefix in ['oo', 'or', 'ro', 'rr']:
                    print(f" bond {prefix}.{pair[0]}.SG {prefix}.{pair[1]}.SG", file=f)

    def _write_heme_bonds(self, f, heme_pair: Tuple[HemeDefinition, HemeDefinition]):
        """Write all heme-related bond definitions"""
        for heme in heme_pair:
            self._write_single_heme_bonds(f, heme)

    def _write_single_heme_bonds(self, f, heme: HemeDefinition):
        """Write bonds for a single heme"""
        print(f"""
#------------------------------------------------------------
#For heme {heme.heme_id}

#Bond ligating atoms to Fe center""", file=f)

        # Write Fe-ligand bonds
        for prefix in ['oo', 'or', 'ro', 'rr']:
            print(f" bond {prefix}.{heme.ligand_proximal}.NE2 {prefix}.{heme.heme_id}.FE", file=f)
            atom_type = "NE2" if heme.is_his_his() else "SD"
            print(f" bond {prefix}.{heme.ligand_distal}.{atom_type} {prefix}.{heme.heme_id}.FE", file=f)

        # Write backbone bonds
        print(f"\n#Bond axially coordinated residues to preceding and proceeding residues", file=f)
        for residue in [heme.ligand_proximal, heme.ligand_distal]:
            for bond in [(residue-1, residue), (residue, residue+1)]:
                if bond not in self.defined_bonds:
                    for prefix in ['oo', 'or', 'ro', 'rr']:
                        print(f" bond {prefix}.{bond[0]}.C {prefix}.{bond[1]}.N", file=f)
                    self.defined_bonds.add(bond)

        # Write thioether bonds for c-type hemes
        if heme.is_c_type():
            print(f"\n#Bond heme thioethers to protein backbone", file=f)
            for prefix in ['oo', 'or', 'ro', 'rr']:
                print(f" bond {prefix}.{heme.cys_b}.CA {prefix}.{heme.heme_id}.CBB2", file=f)
                print(f" bond {prefix}.{heme.cys_c}.CA {prefix}.{heme.heme_id}.CBC1", file=f)

        # Write propionic acid bonds
        print(f"\n#Bond propionic acids to heme", file=f)
        for prefix in ['oo', 'or', 'ro', 'rr']:
            print(f" bond {prefix}.{heme.heme_id}.C2A {prefix}.{heme.heme_id+1}.CA", file=f)
            print(f" bond {prefix}.{heme.heme_id}.C3D {prefix}.{heme.heme_id+2}.CA", file=f)

    def _write_footer(self, f, heme1_id: int, heme2_id: int):
        """Write final commands to save topology and coordinate files"""
        print(f"""
# Save topology and coordinate files
 saveamberparm oo o{heme1_id}-o{heme2_id}.prmtop o{heme1_id}-o{heme2_id}.rst7
 saveamberparm or o{heme1_id}-r{heme2_id}.prmtop o{heme1_id}-r{heme2_id}.rst7
 saveamberparm ro r{heme1_id}-o{heme2_id}.prmtop r{heme1_id}-o{heme2_id}.rst7
 saveamberparm rr r{heme1_id}-r{heme2_id}.prmtop r{heme1_id}-r{heme2_id}.rst7

quit""", file=f)

class PairedChargeAssignment:
    """Main class for handling paired charge assignments"""
    def __init__(self, force_field_dir: str, ff_choice: str, ref_redox_state: str, input_dict: dict):
        try:
            print(f"Initializing PairedChargeAssignment...")
            print(f"Force field directory: {force_field_dir}")
            print(f"Force field choice: {ff_choice}")
            print(f"Reference redox state: {ref_redox_state}")
            self.ff_params = ForceFieldParameters(force_field_dir, ff_choice)
            self.ref_redox_state = RedoxState.OXIDIZED if ref_redox_state == "O" else RedoxState.REDUCED
            self.input_dict = input_dict
            print("Successfully initialized PairedChargeAssignment")
            self.process_heme_pairs()
        except Exception as e:
            print(f"\nError during PairedChargeAssignment initialization:")
            print(f"Error type: {type(e).__name__}")
            print(f"Error message: {str(e)}")
            print(f"Force field dir: {force_field_dir}")
            print(f"FF choice: {ff_choice}")
            print(f"Redox state: {ref_redox_state}")
            raise

    def process_heme_pairs(self):
        """Main processing method"""
        try:
            print("Starting to load heme definitions...")
            heme_definitions = self._load_heme_definitions()
            print("Successfully loaded heme definitions")
        
            print("Validating heme counts...")
            self._validate_heme_counts(heme_definitions)
            print("Successfully validated heme counts")
        
            print("Generating heme pairs...")
            heme_pairs = list(itertools.combinations(heme_definitions, r=2))
            print(f"   There are {len(heme_definitions)} hemes")
            print(f"   There is/are {len(heme_pairs)} unique heme-heme pairwise interactions")
        
            for pair in heme_pairs:
                print(f"\nProcessing heme pair {pair[0].heme_id}-{pair[1].heme_id}...")
                try:
                    self._process_single_pair(pair)
                    print(f"Successfully completed processing pair {pair[0].heme_id}-{pair[1].heme_id}")
                except Exception as e:
                    print(f"Error processing heme pair {pair[0].heme_id}-{pair[1].heme_id}:")
                    print(f"Error type: {type(e).__name__}")
                    print(f"Error message: {str(e)}")
                    raise
                
        except Exception as e:
            print(f"\nError in process_heme_pairs:")
            print(f"Error type: {type(e).__name__}")
            print(f"Error message: {str(e)}")
            print("Full heme definitions and state:")
            if 'heme_definitions' in locals():
                print(f"Number of heme definitions: {len(heme_definitions)}")
                for heme in heme_definitions:
                    print(f"Heme: {heme}")
            else:
                print("Heme definitions not yet loaded")
            raise

    def _load_heme_definitions(self) -> List[HemeDefinition]:
        """Load heme definitions from ResIndexing.txt"""
        if not os.path.isfile("SelResIndexing.txt"):
            raise FileNotFoundError("""
SelResIndexing.txt is missing.
Something went wrong when you defined
the linear sequence of hemes.

This problem must be resolved before 
proceeding.""")
            
        definitions = []
        with open("SelResIndexing.txt") as f:
            for line in f:
                heme_def = self._parse_heme_line(line)
                if heme_def:
                    definitions.append(heme_def)
        return definitions
        
    def _parse_heme_line(self, line: str) -> Optional[HemeDefinition]:
        """Parse a single line from ResIndexing.txt into a HemeDefinition"""
        parts = line.strip().split()
        if len(parts) < 6:
            return None
            
        heme_type = HemeType(parts[-2])
        ligation_type = LigationType(parts[-1])
        
        if heme_type == HemeType.C:
            return HemeDefinition(
                heme_id=int(parts[5]),
                heme_type=heme_type,
                ligation_type=ligation_type,
                cys_b=int(parts[0]),
                cys_c=int(parts[1]),
                ligand_proximal=int(parts[2]),
                ligand_distal=int(parts[3])
            )
        else:
            return HemeDefinition(
                heme_id=int(parts[3]),
                heme_type=heme_type,
                ligation_type=ligation_type,
                ligand_proximal=int(parts[0]),
                ligand_distal=int(parts[1])
            )

    def _validate_heme_counts(self, heme_definitions: List[HemeDefinition]):
        """Validate the counts of different heme types"""
        counts = {
            'c_HH': sum(1 for h in heme_definitions if h.is_c_type() and h.is_his_his()),
            'c_HM': sum(1 for h in heme_definitions if h.is_c_type() and not h.is_his_his()),
            'b_HH': sum(1 for h in heme_definitions if not h.is_c_type() and h.is_his_his()),
            'b_HM': sum(1 for h in heme_definitions if not h.is_c_type() and not h.is_his_his())
        }
        
        total = sum(counts.values())
        if total != len(heme_definitions):
            raise ValueError(f"""
 The total number of hemes ({len(heme_definitions)}) in ResIndexing.txt 
 does NOT equal the sum of c-type His-His ({counts['c_HH']}), 
 b-type His-His ({counts['b_HH']}), c-type His-Met ({counts['c_HM']}),
 and b-type His-Met ({counts['b_HM']}) hemes. These are the only
 types of hemes that can currently be analyzed with BioDC.
 Please revise ResIndexing.txt before re-running BioDC.""")

    def _process_single_pair(self, pair: Tuple[HemeDefinition, HemeDefinition]):
        """Process a single pair of hemes"""
        try:
            h1_id, h2_id = pair[0].heme_id, pair[1].heme_id
            print(f"""
 ------------------------------------------------
 Topologies for Heme Pair {h1_id} - {h2_id} 
 ------------------------------------------------""")
        
            # Generate VMD script
            print(f"Generating VMD script PairInt_{h1_id}-{h2_id}.tcl...")
            vmd_file = f"PairInt_{h1_id}-{h2_id}.tcl"
            vmd_gen = VMDScriptGenerator(vmd_file)
        
            # Generate all redox state combinations
            print("Writing VMD commands for redox state combinations...")
            for state1 in RedoxState:
                for state2 in RedoxState:
                    print(f"  Processing {state1.value}{h1_id}-{state2.value}{h2_id} state...")
                    vmd_gen.add_heme_modification(pair[0], state1)
                    vmd_gen.add_heme_modification(pair[1], state2)
                    vmd_gen.add_final_commands(f"{state1.value}{h1_id}-{state2.value}{h2_id}")
        
            # Close the VMD script
            vmd_gen.finalize()
        
            # Generate TLeap script
            print(f"Generating TLeap input file GeneratePairIntTopologiesForHems{h1_id}-{h2_id}.in...")
            tleap_file = f"GeneratePairIntTopologiesForHems{h1_id}-{h2_id}.in"
            tleap_gen = TLeapScriptGenerator(self.ff_params)
            tleap_gen.generate_script(pair, tleap_file)
        
            # Run external tools
            self._run_external_tools(h1_id, h2_id)
        except Exception as e:
            print(f"\nError in _process_single_pair: {str(e)}")
            raise

    def _run_external_tools(self, h1_id: int, h2_id: int):
            """Run VMD and TLeap with error checking"""
            # Run VMD
            print(f"\nUsing VMD to generate redox microstate PDBs for Heme-{h1_id} and Heme-{h2_id}...")
            vmd_command = f"vmd -e PairInt_{h1_id}-{h2_id}.tcl > PairInt_{h1_id}-{h2_id}.log"
            subprocess.run(vmd_command, shell=True)

            # Process the generated PDB files to add TER records
            print("\nProcessing PDB files to add TER records...")
            for state1 in RedoxState:
                for state2 in RedoxState:
                    pdb_file = f"{state1.value}{h1_id}-{state2.value}{h2_id}.pdb"
                    if os.path.exists(pdb_file):
                        print(f"  Adding TER records to {pdb_file}...")
                        temp_pdb = f"{state1.value}{h1_id}-{state2.value}{h2_id}_temp.pdb"
                        PDBProcessor.insert_ter_records(pdb_file, temp_pdb)
                        os.replace(temp_pdb, pdb_file)  # Replace original with processed version

            # Check VMD results
            errors = []
            for state1 in RedoxState:
                for state2 in RedoxState:
                    final_pdb = f"{state1.value}{h1_id}-{state2.value}{h2_id}.pdb"
                
                    if os.path.isfile(final_pdb):
                        print(f"  √ {state1.value}{h1_id} -- {state2.value}{h2_id}: VMD finished successfully!")
                    else:
                        errors.append(f"Failed to generate {final_pdb}")

            if errors:
                error_msg = "\n".join(errors)
                raise RuntimeError(f"""
    VMD failed to generate the PDB for one or more redox microstates:
    {error_msg}
    Please inspect PairInt_{h1_id}-{h2_id}.log before re-running this module""")

            # Run TLeap
            print(f"\nUsing TLEaP to build the redox microstate topologies for Heme-{h1_id} and Heme-{h2_id}...")
            tleap_command = f"tleap -s -f GeneratePairIntTopologiesForHems{h1_id}-{h2_id}.in > GeneratePairIntTopologiesForHems{h1_id}-{h2_id}.log"
            subprocess.run(tleap_command, shell=True)

            # Check TLeap results
            errors = []
            for state1 in RedoxState:
                for state2 in RedoxState:
                    prmtop = f"{state1.value}{h1_id}-{state2.value}{h2_id}.prmtop"
                    rst7 = f"{state1.value}{h1_id}-{state2.value}{h2_id}.rst7"
                
                    if os.path.isfile(prmtop) and os.path.isfile(rst7):
                        print(f"  √ {state1.value}{h1_id} -- {state2.value}{h2_id}: TLEaP finished successfully!")
                    else:
                        errors.append(f"Missing {prmtop} or {rst7}")

            if errors:
                error_msg = "\n".join(errors)
                raise RuntimeError(f"""
    TLEaP failed to build the topologies for one or more redox microstates:
    {error_msg}
    Please inspect GeneratePairIntTopologiesForHems{h1_id}-{h2_id}.log before re-running this module""")
        
def main(force_field_dir: str, ff_choice: str, ref_redox_state: str, input_dict: dict):
    """Main entry point for the paired charge assignment process"""
    try:
        processor = PairedChargeAssignment(force_field_dir, ff_choice, ref_redox_state, input_dict)
        processor.process_heme_pairs()
    except Exception as e:
        print(f"\nError: {str(e)}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python PairedChargeAssignment.py <force_field_dir> <ff_choice> <ref_redox_state>")
        sys.exit(1)

    # When called from command line, use empty dict
    main(sys.argv[1], sys.argv[2], sys.argv[3], {})
