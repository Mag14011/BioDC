"""
Module for generating and managing redox states of heme-containing proteins.
Handles both single-heme state changes and heme-pair microstates.
Supports all heme types (b and c) and various ligations.
"""

import copy
import logging
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set, Any, Union
from enum import Enum
from dataclasses import dataclass

from Bio.PDB import PDBIO, PDBParser, Structure, Model, Chain, Residue
from Bio.PDB.PDBIO import Select

logger = logging.getLogger(__name__)

class HemeType(Enum):
    B = "b"  # b-type hemes
    C = "c"  # c-type hemes

class LigationType(Enum):
    HIS_HIS = "HH"  # His-His
    HIS_MET = "HM"  # His-Met
    HIS_CYS = "HC"  # His-Cys
    HIS_TYR = "HY"  # His-Tyr 
    HIS_ASP = "HD"  # His-Asp
    HIS_GLU = "HE"  # His-Glu
    HIS_ASN = "HN"  # His-Asn
    HIS_GLN = "HQ"  # His-Gln
    HIS_LYS = "HK"  # His-Lys

class RedoxState(Enum):
    OXIDIZED = "ox"
    REDUCED = "red"

class ForceFieldParameters:
    """Manages force field parameters and their loading for heme calculations."""
    
    def __init__(self, forcefield_dir: Path):
        """
        Initialize ForceFieldParameters with directory containing parameter files.
        
        Args:
            forcefield_dir: Path to directory containing forcefield files
        """
        self.dir = Path(forcefield_dir)
        self._available_params: Set[Tuple[HemeType, LigationType, RedoxState]] = set()
        self._scan_available_parameters()
#       print("\nInitializing ForceFieldParameters...")
#       print("Forcefield directory:", forcefield_dir)

    def _scan_available_parameters(self):
        """Scan forcefield directory to determine which parameter sets exist."""
#       print("\nScanning available parameters...")

        patterns = {
            (HemeType.B, LigationType.HIS_HIS): "HisHisLigated_b-heme",
            (HemeType.B, LigationType.HIS_MET): "HisMetLigated_b-heme",
            (HemeType.B, LigationType.HIS_CYS): "HisCysLigated_b-heme",
            (HemeType.B, LigationType.HIS_TYR): "HisTyrLigated_b-heme",
            (HemeType.B, LigationType.HIS_ASP): "HisAspLigated_b-heme",
            (HemeType.B, LigationType.HIS_GLU): "HisGluLigated_b-heme",
            (HemeType.B, LigationType.HIS_ASN): "HisAsnLigated_b-heme",
            (HemeType.B, LigationType.HIS_GLN): "HisGlnLigated_b-heme",
            (HemeType.B, LigationType.HIS_LYS): "HisLysLigated_b-heme",
            (HemeType.C, LigationType.HIS_HIS): "HisHisLigated_c-heme",
            (HemeType.C, LigationType.HIS_MET): "HisMetLigated_c-heme",
            (HemeType.C, LigationType.HIS_CYS): "HisCysLigated_c-heme",
            (HemeType.C, LigationType.HIS_TYR): "HisTyrLigated_c-heme",
            (HemeType.C, LigationType.HIS_ASP): "HisAspLigated_c-heme",
            (HemeType.C, LigationType.HIS_GLU): "HisGluLigated_c-heme",
            (HemeType.C, LigationType.HIS_ASN): "HisAsnLigated_c-heme",
            (HemeType.C, LigationType.HIS_GLN): "HisGlnLigated_c-heme",
            (HemeType.C, LigationType.HIS_LYS): "HisLysLigated_c-heme",
        }

        for (heme_type, ligation), pattern in patterns.items():
            for redox in RedoxState:
                lib_name = f"{'Oxidized' if redox == RedoxState.OXIDIZED else 'Reduced'}_{pattern}_RESP.lib"
                if (self.dir / lib_name).exists():
                    self._available_params.add((heme_type, ligation, redox))
#       print("Available parameters:", self._available_params)

    def parameter_exists(self, heme_type: HemeType, ligation: LigationType, redox: RedoxState) -> bool:
        """
        Check if parameters exist for given heme configuration.
        
        Args:
            heme_type: Type of heme (b or c)
            ligation: Type of ligation
            redox: Redox state
            
        Returns:
            Whether parameters exist for this configuration
        """
        return (heme_type, ligation, redox) in self._available_params

    def get_lib_name(self, heme_type: HemeType, ligation: Union[str, LigationType], 
                    redox: RedoxState) -> str:
        """Get library file name for given heme configuration."""
        # Convert ligation to enum if needed
        ligation = self._convert_ligand_type(ligation)
        
        pattern = next(
            pattern for (ht, lig), pattern in {
                (HemeType.B, LigationType.HIS_HIS): "HisHisLigated_b-heme",
                (HemeType.B, LigationType.HIS_MET): "HisMetLigated_b-heme",
                (HemeType.B, LigationType.HIS_CYS): "HisCysLigated_b-heme",
                (HemeType.B, LigationType.HIS_TYR): "HisTyrLigated_b-heme",
                (HemeType.B, LigationType.HIS_ASP): "HisAspLigated_b-heme",
                (HemeType.B, LigationType.HIS_GLU): "HisGluLigated_b-heme",
                (HemeType.B, LigationType.HIS_ASN): "HisAsnLigated_b-heme",
                (HemeType.B, LigationType.HIS_GLN): "HisGlnLigated_b-heme",
                (HemeType.B, LigationType.HIS_LYS): "HisLysLigated_b-heme",
                (HemeType.C, LigationType.HIS_HIS): "HisHisLigated_c-heme",
                (HemeType.C, LigationType.HIS_MET): "HisMetLigated_c-heme",
                (HemeType.C, LigationType.HIS_CYS): "HisCysLigated_c-heme",
                (HemeType.C, LigationType.HIS_TYR): "HisTyrLigated_c-heme",
                (HemeType.C, LigationType.HIS_ASP): "HisAspLigated_c-heme",
                (HemeType.C, LigationType.HIS_GLU): "HisGluLigated_c-heme",
                (HemeType.C, LigationType.HIS_ASN): "HisAsnLigated_c-heme",
                (HemeType.C, LigationType.HIS_GLN): "HisGlnLigated_c-heme",
                (HemeType.C, LigationType.HIS_LYS): "HisLysLigated_c-heme",
            }.items()
            if ht == heme_type and lig == ligation
        )
        return f"{'Oxidized' if redox == RedoxState.OXIDIZED else 'Reduced'}_{pattern}_RESP.lib"

    def get_frcmod_name(self, heme_type: HemeType, ligation: Union[str, LigationType], 
                       redox: RedoxState) -> str:
        """Get force field modification file name for given heme configuration."""
        # Convert ligation to enum if needed
        ligation = self._convert_ligand_type(ligation)
        return self.get_lib_name(heme_type, ligation, redox).replace('_RESP.lib', '.frcmod')

    def get_atom_types(self, heme_type: HemeType, ligation: Union[str, LigationType], 
                      redox: RedoxState) -> Optional[Dict[str, str]]:
        """Get atom type definitions for the heme configuration."""
        # Convert ligation to enum if needed
        ligation = self._convert_ligand_type(ligation)
        
        atom_types = {
            # B-type hemes
            (HemeType.B, LigationType.HIS_HIS, RedoxState.OXIDIZED): {
                'metal': 'M1', 'nitrogens': ['Y1', 'Y2', 'Y3', 'Y4', 'Y5', 'Y6']
            },
            (HemeType.B, LigationType.HIS_HIS, RedoxState.REDUCED): {
                'metal': 'M2', 'nitrogens': ['Z1', 'Z2', 'Z3', 'Z4', 'Z5', 'Z6']
            },
            (HemeType.B, LigationType.HIS_MET, RedoxState.OXIDIZED): {
                'metal': 'M3', 'nitrogens': ['Y7', 'Y8', 'Y9', 'Y10', 'Y11'],
                'sulfurs': ['Y12']
            },
            (HemeType.B, LigationType.HIS_MET, RedoxState.REDUCED): {
                'metal': 'M4', 'nitrogens': ['Z7', 'Z8', 'Z9', 'Z10', 'Z11'],
                'sulfurs': ['Z12']
            },
            
            # C-type hemes
            (HemeType.C, LigationType.HIS_HIS, RedoxState.OXIDIZED): {
                'metal': 'M5', 'nitrogens': ['Y13', 'Y14', 'Y15', 'Y16', 'Y17', 'Y18']
            },
            (HemeType.C, LigationType.HIS_HIS, RedoxState.REDUCED): {
                'metal': 'M6', 'nitrogens': ['Z13', 'Z14', 'Z15', 'Z16', 'Z17', 'Z18']
            },
            (HemeType.C, LigationType.HIS_MET, RedoxState.OXIDIZED): {
                'metal': 'M7', 'nitrogens': ['Y19', 'Y20', 'Y21', 'Y22', 'Y23'],
                'sulfurs': ['Y24']
            },
            (HemeType.C, LigationType.HIS_MET, RedoxState.REDUCED): {
                'metal': 'M8', 'nitrogens': ['Z19', 'Z20', 'Z21', 'Z22', 'Z23'],
                'sulfurs': ['Z24']
            }
        }
        return atom_types.get((heme_type, ligation, redox))

    def get_ligand_atom_type(self, ligation: LigationType) -> str:
        """
        Get the ligating atom type for a given ligation.
        
        Args:
            ligation: Type of ligation
            
        Returns:
            Atom type string for the ligating atom
        """
        return {
            LigationType.HIS_HIS: "NE2",  # His-His
            LigationType.HIS_MET: "SD",   # His-Met
            LigationType.HIS_CYS: "SG",   # His-Cys
            LigationType.HIS_TYR: "OH",   # His-Tyr
            LigationType.HIS_ASP: "OD1",  # His-Asp
            LigationType.HIS_GLU: "OE1",  # His-Glu
            LigationType.HIS_ASN: "ND2",  # His-Asn
            LigationType.HIS_GLN: "NE2",  # His-Gln
            LigationType.HIS_LYS: "NZ",   # His-Lys
        }[ligation]

    def _convert_ligand_type(self, ligand_type: Union[str, LigationType]) -> LigationType:
        """Convert ligand type string to LigationType enum if needed."""
        if isinstance(ligand_type, LigationType):
            return ligand_type
            
        # If it's a string, ensure it's in the correct format
        if ligand_type in [e.value for e in LigationType]:
            return LigationType(ligand_type)
        
        # Try converting from full name to short code
        full_to_short = {
            'HIS': 'HH',
            'MET': 'HM',
            'CYS': 'HC',
            'TYR': 'HY',
            'ASP': 'HD',
            'GLU': 'HE',
            'ASN': 'HN',
            'GLN': 'HQ',
            'LYS': 'HK'
        }
        
        if ligand_type in full_to_short:
            return LigationType(full_to_short[ligand_type])
            
        raise ValueError(f"Invalid ligand type: {ligand_type}")

class TLeapScriptGenerator:
    """Generates TLeap input scripts for topology generation."""
    def __init__(self, ff_params: ForceFieldParameters, launch_dir: Path):
        self.ff_params = ff_params
        self.launch_dir = launch_dir
        self.structure_info = {}  # Initialize empty dict
        self.defined_bonds = set()
        self.written_atom_types = set()
        self.written_ff_params = set()

    def generate_single_heme_script(self,
                                  reference_hemes: Dict[int, Dict],
                                  selected_heme: Dict,
                                  output_file: Path,
                                  reference_state: RedoxState) -> bool:
        """Generate TLeap script for single selected heme in both states."""
        try:
            # Create structure info dict for all hemes
            self.structure_info = {}
        
            # Add selected heme
            heme_key = f"{selected_heme['distal_type']}-{selected_heme['type']}"
            self.structure_info[heme_key] = {
                'ox': 1 if reference_state == RedoxState.OXIDIZED else 0,
                'red': 1 if reference_state == RedoxState.REDUCED else 0
            }
        
            # Add reference hemes
            for ref_heme in reference_hemes.values():
                ref_key = f"{ref_heme['distal_type']}-{ref_heme['type']}"
                if ref_key not in self.structure_info:
                    self.structure_info[ref_key] = {
                        'ox': 0,
                        'red': 0
                    }
                # Set count for reference state
                if reference_state == RedoxState.OXIDIZED:
                    self.structure_info[ref_key]['ox'] += 1
                else:
                    self.structure_info[ref_key]['red'] += 1

            self.defined_bonds = set()
            self.written_atom_types = set()
            self.written_ff_params = set()

            with open(output_file, 'w') as f:
                self._write_header(f)
                
                # Write atom types for all unique hemes
                self._write_atom_types(f, reference_hemes, selected_heme, reference_state)
                
                # Write forcefield loads
                self._write_forcefield_loads(f, reference_hemes, selected_heme, reference_state)
                
                # Load both PDB structures
                self._write_structure_loads_single(f, selected_heme['heme_id'])
                
                # Write bonds for all hemes in both structures
                self._write_bonds_single(f, reference_hemes, selected_heme)
                
                self._write_footer_single(f, selected_heme['heme_id'])
                
            return True
            
        except Exception as e:
            logging.error(f"Failed to generate TLeap script: {str(e)}")
            raise
            
    def generate_pair_script(self,
                           reference_hemes: Dict[int, Dict],
                           heme1_config: Dict,
                           heme2_config: Dict,
                           output_file: Path,
                           reference_state: RedoxState) -> bool:
        """Generate TLeap script for heme pair in all four states."""
        try:
            self.defined_bonds = set()
            self.written_atom_types = set()
            self.written_ff_params = set()

            with open(output_file, 'w') as f:
                self._write_header(f)
                
                # Write atom types for all unique hemes
                self._write_atom_types_pair(f, reference_hemes, 
                                          heme1_config, heme2_config, reference_state)
                
                # Write forcefield loads
                self._write_forcefield_loads_pair(f, reference_hemes, 
                                                heme1_config, heme2_config, reference_state)
                
                # Load all four PDB structures
                self._write_structure_loads_pair(f, heme1_config['heme_id'], 
                                               heme2_config['heme_id'])
                
                # Write bonds for all hemes in all structures
                self._write_bonds_pair(f, reference_hemes, heme1_config, heme2_config)
                
                self._write_footer_pair(f, heme1_config['heme_id'], 
                                      heme2_config['heme_id'])
                
            return True
            
        except Exception as e:
            logging.error(f"Failed to generate TLeap script: {str(e)}")
            raise
            
    def _write_header(self, f):
        """Write initial TLeap commands."""
        print("""
# Load base force field parameters
source leaprc.constph
source leaprc.conste
source leaprc.gaff
source leaprc.water.tip3p

addAtomTypes {""", file=f)

    def _write_atom_types(self, f, reference_hemes: Dict[int, Dict],
                         selected_heme: Dict, reference_state: RedoxState):
        """Write all necessary atom type definitions."""
        # Write reference heme atom types (in reference state only)
#       print("Structure info in _write_atom_types:", self.structure_info)

        for heme_id, config in reference_hemes.items():
            heme_type = HemeType(config['type'])
            ligation = LigationType(config['distal_type'])
            key = (heme_type, ligation, reference_state)
            
            if key not in self.written_atom_types:
                atom_types = self.ff_params.get_atom_types(heme_type, ligation, reference_state)
                if atom_types:
                    self._write_heme_atom_types(f, atom_types)
                    self.written_atom_types.add(key)
        
        # Write selected heme atom types (both states)
        heme_type = HemeType(selected_heme['type'])
        ligation = LigationType(selected_heme['distal_type'])
        
        for state in [RedoxState.OXIDIZED, RedoxState.REDUCED]:
            key = (heme_type, ligation, state)
            if key not in self.written_atom_types:
                atom_types = self.ff_params.get_atom_types(heme_type, ligation, state)
                if atom_types:
                    self._write_heme_atom_types(f, atom_types)
                    self.written_atom_types.add(key)
        
        print("}", file=f)
        print("", file=f)

    def _write_heme_atom_types(self, f, atom_types: Dict[str, Any]):
        """Write atom types for a specific heme configuration."""
        # Write metal atom type
        print(f'    {{ "{atom_types["metal"]}"  "Fe" "sp3" }}', file=f)
        
        # Write nitrogen atom types
        for n_type in atom_types['nitrogens']:
            print(f'    {{ "{n_type}"  "N" "sp3" }}', file=f)
            
        # Write oxygen atom types if present
        if 'oxygens' in atom_types:
            for o_type in atom_types['oxygens']:
                print(f'    {{ "{o_type}"  "O" "sp3" }}', file=f)
                
        # Write sulfur atom types if present
        if 'sulfurs' in atom_types:
            for s_type in atom_types['sulfurs']:
                print(f'    {{ "{s_type}"  "S" "sp3" }}', file=f)

    def _write_forcefield_loads(self, f, reference_hemes: Dict[int, Dict],
                              selected_heme: Dict, reference_state: RedoxState):
        """Write commands to load all necessary forcefield parameters."""
        print("\n# Load forcefield parameters", file=f)

#       print("\nWriting forcefield loads:")
#       print("Selected heme:", selected_heme)
#       print("Reference hemes:", reference_hemes)
#       print("Reference state:", reference_state)

#       for state in [RedoxState.OXIDIZED, RedoxState.REDUCED]:
#           print(f"Checking state: {state}")

        # Load parameters for reference hemes (reference state only)
        for config in reference_hemes.values():
            heme_type = HemeType(config['type'])
            ligation = LigationType(config['distal_type'])
            key = (heme_type, ligation, reference_state)
            
            if key not in self.written_ff_params:
                lib_name = self.ff_params.get_lib_name(heme_type, ligation, reference_state)
                frcmod_name = self.ff_params.get_frcmod_name(heme_type, ligation, reference_state)
#               print(f"loadamberparams {self.ff_params.dir.absolute()}/{frcmod_name}")
#               print(f"loadoff {self.ff_params.dir.absolute()}/{lib_name}")
                print(f"loadamberparams {self.ff_params.dir.absolute()}/{frcmod_name}", file=f)
                print(f"loadoff {self.ff_params.dir.absolute()}/{lib_name}", file=f)
                self.written_ff_params.add(key)
        
        # Load parameters for selected heme (both states)
        heme_type = HemeType(selected_heme['type'])
        ligation = LigationType(selected_heme['distal_type'])
        
        for state in [RedoxState.OXIDIZED, RedoxState.REDUCED]:
            key = (heme_type, ligation, state)
            if key not in self.written_ff_params:
                lib_name = self.ff_params.get_lib_name(heme_type, ligation, state)
                frcmod_name = self.ff_params.get_frcmod_name(heme_type, ligation, state)
                print(f"loadamberparams {self.ff_params.dir.absolute()}/{frcmod_name}", file=f)
                print(f"loadoff {self.ff_params.dir.absolute()}/{lib_name}", file=f)
                self.written_ff_params.add(key)

    def _write_structure_loads_single(self, f, heme_id: int):
        """Write commands to load PDB files for single heme case."""
        print(f"""
# Load PDB files for oxidized and reduced states
ox = loadpdb {self.launch_dir.absolute()}/EE/o_{heme_id}_{heme_id}.pdb
red = loadpdb {self.launch_dir.absolute()}/EE/r_{heme_id}_{heme_id}.pdb""", file=f)

    def _write_structure_loads_pair(self, f, heme1_id: int, heme2_id: int):
        """Write commands to load PDB files for heme pair case."""
        print(f"""
# Load PDB files for all microstates
oo = loadpdb {self.launch_dir.absolute()}/EE/oo_{heme1_id}_{heme2_id}.pdb
or = loadpdb {self.launch_dir.absolute()}/EE/or_{heme1_id}_{heme2_id}.pdb
ro = loadpdb {self.launch_dir.absolute()}/EE/ro_{heme1_id}_{heme2_id}.pdb
rr = loadpdb {self.launch_dir.absolute()}/EE/rr_{heme1_id}_{heme2_id}.pdb""", file=f)

    def _write_bonds_single(self, f, reference_hemes: Dict[int, Dict], selected_heme: Dict):
        """Write all bond definitions for single heme case."""
        print("\n# Bond definitions", file=f)
        
        # Write bonds for reference hemes in both structures
        for config in reference_hemes.values():
            self._write_heme_bonds(f, config, structures=['ox', 'red'])
            
        # Write bonds for selected heme in both structures
        self._write_heme_bonds(f, selected_heme, structures=['ox', 'red'])

    def _write_bonds_pair(self, f, reference_hemes: Dict[int, Dict], 
                         heme1_config: Dict, heme2_config: Dict):
        """Write all bond definitions for heme pair case."""
        print("\n# Bond definitions", file=f)
        
        # Write bonds for reference hemes in all structures
        for config in reference_hemes.values():
            self._write_heme_bonds(f, config, structures=['oo', 'or', 'ro', 'rr'])
            
        # Write bonds for selected hemes in all structures
        self._write_heme_bonds(f, heme1_config, structures=['oo', 'or', 'ro', 'rr'])
        self._write_heme_bonds(f, heme2_config, structures=['oo', 'or', 'ro', 'rr'])

    def _make_canonical_bond(self, res1_id: int, res2_id: int, atom1: str, atom2: str) -> Tuple[int, int, str, str]:
        """
        Create a canonical representation of a bond that is direction-independent.
        Always puts the smaller residue ID first, or if residue IDs are equal,
        puts the alphabetically first atom type first.
        
        Args:
            res1_id: First residue ID
            res2_id: Second residue ID
            atom1: First atom type
            atom2: Second atom type
            
        Returns:
            Tuple representing the bond in canonical form
        """
        if res1_id < res2_id:
            return (res1_id, res2_id, atom1, atom2)
        elif res2_id < res1_id:
            return (res2_id, res1_id, atom2, atom1)
        else:  # res1_id == res2_id
            if atom1 <= atom2:
                return (res1_id, res2_id, atom1, atom2)
            else:
                return (res1_id, res2_id, atom2, atom1)

    def _write_heme_bonds(self, f, heme_config: Dict, structures: List[str]):
        """
        Write all bonds for a specific heme, avoiding duplicates within each structure.

        Args:
            f: File handle
            heme_config: Heme configuration dictionary
            structures: List of structure prefixes ('ox'/'red' or 'oo'/'or'/'ro'/'rr')
        """
        heme_id = heme_config['heme_id']
        is_c_type = heme_config['type'] == 'c'
        ligation = LigationType(heme_config['distal_type'])
        print(f"\n# Bond definitions for heme {heme_id}", file=f)

        # Reset defined_bonds for each structure set
        structure_defined_bonds = {struct: set() for struct in structures}

        for struct in structures:
            # Create tuples to represent each bond
            bonds_to_write = []

            # Proximal histidine bond
            prox_bond = self._make_canonical_bond(
                heme_config['his_p'], heme_id, 'NE2', 'FE'
            )
            if prox_bond not in structure_defined_bonds[struct]:
                bonds_to_write.append(
                    (f"{struct}.{heme_config['his_p']}.NE2", f"{struct}.{heme_id}.FE")
                )
                structure_defined_bonds[struct].add(prox_bond)

            # Distal ligand bond
            atom_type = self.ff_params.get_ligand_atom_type(ligation)
            dist_bond = self._make_canonical_bond(
                heme_config['distal'], heme_id, atom_type, 'FE'
            )
            if dist_bond not in structure_defined_bonds[struct]:
                bonds_to_write.append(
                    (f"{struct}.{heme_config['distal']}.{atom_type}", f"{struct}.{heme_id}.FE")
                )
                structure_defined_bonds[struct].add(dist_bond)

            # C-type specific bonds
            if is_c_type and 'cys_b' in heme_config and 'cys_c' in heme_config:
                cys_b_bond = self._make_canonical_bond(
                    heme_config['cys_b'], heme_id, 'CA', 'CBB2'
                )
                if cys_b_bond not in structure_defined_bonds[struct]:
                    bonds_to_write.append(
                        (f"{struct}.{heme_config['cys_b']}.CA", f"{struct}.{heme_id}.CBB2")
                    )
                    structure_defined_bonds[struct].add(cys_b_bond)

                cys_c_bond = self._make_canonical_bond(
                    heme_config['cys_c'], heme_id, 'CA', 'CBC1'
                )
                if cys_c_bond not in structure_defined_bonds[struct]:
                    bonds_to_write.append(
                        (f"{struct}.{heme_config['cys_c']}.CA", f"{struct}.{heme_id}.CBC1")
                    )
                    structure_defined_bonds[struct].add(cys_c_bond)

            # Propionic acid bonds
            prop_a_bond = self._make_canonical_bond(
                heme_id, heme_id+1, 'C2A', 'CA'
            )
            if prop_a_bond not in structure_defined_bonds[struct]:
                bonds_to_write.append(
                    (f"{struct}.{heme_id}.C2A", f"{struct}.{heme_id+1}.CA")
                )
                structure_defined_bonds[struct].add(prop_a_bond)

            prop_d_bond = self._make_canonical_bond(
                heme_id, heme_id+2, 'C3D', 'CA'
            )
            if prop_d_bond not in structure_defined_bonds[struct]:
                bonds_to_write.append(
                    (f"{struct}.{heme_id}.C3D", f"{struct}.{heme_id+2}.CA")
                )
                structure_defined_bonds[struct].add(prop_d_bond)

            # Write all new bonds for this structure
            for atom1, atom2 in bonds_to_write:
                print(f"bond {atom1} {atom2}", file=f)

    def _write_footer_single(self, f, heme_id: int):
        """Write commands to save topology and coordinate files for single heme case."""
        print(f"""
# Save topology and coordinate files
saveamberparm ox {self.launch_dir.absolute()}/EE/o_{heme_id}_{heme_id}.prmtop {self.launch_dir.absolute()}/EE/o_{heme_id}_{heme_id}.rst7
saveamberparm red {self.launch_dir.absolute()}/EE/r_{heme_id}_{heme_id}.prmtop {self.launch_dir.absolute()}/EE/r_{heme_id}_{heme_id}.rst7

quit""", file=f)

    def _write_footer_pair(self, f, heme1_id: int, heme2_id: int):
        """Write commands to save topology and coordinate files for heme pair case."""
        print(f"""
# Save topology and coordinate files
saveamberparm oo {self.launch_dir.absolute()}/EE/oo_{heme1_id}_{heme2_id}.prmtop {self.launch_dir.absolute()}/EE/oo_{heme1_id}_{heme2_id}.rst7
saveamberparm or {self.launch_dir.absolute()}/EE/or_{heme1_id}_{heme2_id}.prmtop {self.launch_dir.absolute()}/EE/or_{heme1_id}_{heme2_id}.rst7
saveamberparm ro {self.launch_dir.absolute()}/EE/ro_{heme1_id}_{heme2_id}.prmtop {self.launch_dir.absolute()}/EE/ro_{heme1_id}_{heme2_id}.rst7
saveamberparm rr {self.launch_dir.absolute()}/EE/rr_{heme1_id}_{heme2_id}.prmtop {self.launch_dir.absolute()}/EE/rr_{heme1_id}_{heme2_id}.rst7

quit""", file=f)

    def _write_atom_types_pair(self, f, reference_hemes: Dict[int, Dict],
                              heme1_config: Dict, heme2_config: Dict,
                              reference_state: RedoxState):
        """Write all necessary atom type definitions for pair case."""
#       print("Structure info in _write_atom_types:", self.structure_info) 

        # Write reference heme atom types (in reference state only)
        for heme_id, config in reference_hemes.items():
            heme_type = HemeType(config['type'])
            ligation = LigationType(config['distal_type'])
            key = (heme_type, ligation, reference_state)
            
            if key not in self.written_atom_types:
                atom_types = self.ff_params.get_atom_types(heme_type, ligation, reference_state)
                if atom_types:
                    self._write_heme_atom_types(f, atom_types)
                    self.written_atom_types.add(key)
        
        # Write selected hemes atom types (both states for each)
        for config in [heme1_config, heme2_config]:
            heme_type = HemeType(config['type'])
            ligation = LigationType(config['distal_type'])
            
            for state in [RedoxState.OXIDIZED, RedoxState.REDUCED]:
                key = (heme_type, ligation, state)
                if key not in self.written_atom_types:
                    atom_types = self.ff_params.get_atom_types(heme_type, ligation, state)
                    if atom_types:
                        self._write_heme_atom_types(f, atom_types)
                        self.written_atom_types.add(key)
        
        print("}", file=f)
        print("", file=f)

    def _write_forcefield_loads_pair(self, f, reference_hemes: Dict[int, Dict],
                                   heme1_config: Dict, heme2_config: Dict,
                                   reference_state: RedoxState):
        """Write commands to load all necessary forcefield parameters for pair case."""
        print("\n# Load forcefield parameters", file=f)
        
        # Load parameters for reference hemes (reference state only)
        for config in reference_hemes.values():
            heme_type = HemeType(config['type'])
            ligation = LigationType(config['distal_type'])
            key = (heme_type, ligation, reference_state)
            
            if key not in self.written_ff_params:
                lib_name = self.ff_params.get_lib_name(heme_type, ligation, reference_state)
                frcmod_name = self.ff_params.get_frcmod_name(heme_type, ligation, reference_state)
                print(f"loadamberparams {self.ff_params.dir.absolute()}/{frcmod_name}", file=f)
                print(f"loadoff {self.ff_params.dir.absolute()}/{lib_name}", file=f)
                self.written_ff_params.add(key)
        
        # Load parameters for selected hemes (both states for each)
        for config in [heme1_config, heme2_config]:
            heme_type = HemeType(config['type'])
            ligation = LigationType(config['distal_type'])
            
            for state in [RedoxState.OXIDIZED, RedoxState.REDUCED]:
                key = (heme_type, ligation, state)
                if key not in self.written_ff_params:
                    lib_name = self.ff_params.get_lib_name(heme_type, ligation, state)
                    frcmod_name = self.ff_params.get_frcmod_name(heme_type, ligation, state)
                    print(f"loadamberparams {self.ff_params.dir.absolute()}/{frcmod_name}", file=f)
                    print(f"loadoff {self.ff_params.dir.absolute()}/{lib_name}", file=f)
                    self.written_ff_params.add(key)

class HemeStateGenerator:
    """Generates PDB files for different redox states of heme-containing proteins."""

    class RefStateSelect(Select):
        """Select class for filtering water and ions from PDB files."""
        def accept_residue(self, residue):
            return residue.get_resname() not in ['WAT', 'Na+', 'Cl-']
    
    def __init__(self, input_pdb_path: str, launch_dir: Path):
        """
        Initialize HemeStateGenerator with an input PDB file and launch directory.
        
        Args:
            input_pdb_path (str): Path to the input PDB file
            launch_dir (Path): Path to the launch directory containing SPR subdirectory
        """
        self.parser = PDBParser(QUIET=True)
        self.structure = self.parser.get_structure('input', input_pdb_path)
        self.launch_dir = Path(launch_dir)  # Ensure it's a Path object
        
        # Nomenclature mappings directly matching the HemeNomenclature class
        self.DISTAL_LIGANDS = {
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

        # Heme naming patterns (c-type)
        self.CTYPE_HEME_NAMES = {
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
        self.BTYPE_HEME_NAMES = {
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
        
        # Proximal His naming (c-type)
        self.CTYPE_PROX_HIS = {
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
        
        # Proximal His naming (b-type)
        self.BTYPE_PROX_HIS = {
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
        
        # Distal ligand naming (c-type)
        self.CTYPE_DISTAL = {
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
        
        # Distal ligand naming (b-type)
        self.BTYPE_DISTAL = {
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

        self.SHORT_TO_FULL = {
            'HH': 'HIS',  # His-His
            'HM': 'MET',  # His-Met
            'HC': 'CYS',  # His-Cys
            'HY': 'TYR',  # His-Tyr
            'HD': 'ASP',  # His-Asp
            'HE': 'GLU',  # His-Glu
            'HN': 'ASN',  # His-Asn
            'HQ': 'GLN',  # His-Gln
            'HK': 'LYS'   # His-Lys
        }
        
        self.FULL_TO_SHORT = {v: k for k, v in self.SHORT_TO_FULL.items()}

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
        
        # Combine heme names from both b-type and c-type dictionaries
        all_heme_names = set()
        for ligand_dict in [self.BTYPE_HEME_NAMES, self.CTYPE_HEME_NAMES]:
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
                # Use the current line's atom serial number to generate TER record
                try:
                    last_serial = int(line[6:11])
#                   ter_line = f"TER{' ' * 6}{last_serial + 1:6d}\n"
                    ter_line = f"TER{' ' * 6}                    \n"
                except ValueError:
                    ter_line = "TER\n"
                
                new_lines.append(ter_line)
        
        # Ensure a final TER record if none exists
        if not any(line.startswith('TER') for line in new_lines):
            new_lines.append("TER\n")
        
        # Write modified PDB
        with open(output_pdb, 'w') as f:
            f.writelines(new_lines)

    def parse_heme_environments(self, heme_ids: List[int]) -> List[Dict]:
        """
        Parse ResIndexing.txt to extract configuration for specified heme IDs.
        
        Format:
        b-type: his_p distal original_heme_id new_heme_id type ligand_type
        c-type: CysB CysC his_p distal original_heme_id new_heme_id type ligand_type
    
        Args:
            heme_ids (List[int]): List of heme residue IDs to extract
    
        Returns:
            List of configuration dictionaries for the specified hemes
    
        Raises:
            ValueError if any specified heme ID is not found in ResIndexing.txt
        """
        spr_dir = self.launch_dir / "SPR"
        indexing_file = spr_dir / "ResIndexing.txt"
    
        if not indexing_file.exists():
            raise FileNotFoundError(f"ResIndexing.txt not found in {spr_dir}")
    
        heme_configs = []
        found_heme_ids = set()
    
        with open(indexing_file, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if not parts:
                    continue
            
                # Determine heme type and parse accordingly based on number of fields
                if len(parts) == 6:  # b-type heme
                    his_p, distal = int(parts[0]), int(parts[1])
                    orig_heme_id = int(parts[2])
                    heme_id = int(parts[3])  # new heme ID
                    heme_type = parts[4]  # 'b'
                    ligand_type = parts[5]  # 'HX'
                    is_c_type = False
                    cys_b = cys_c = None
                elif len(parts) == 8:  # c-type heme
                    cys_b, cys_c = int(parts[0]), int(parts[1])
                    his_p, distal = int(parts[2]), int(parts[3])
                    orig_heme_id = int(parts[4])
                    heme_id = int(parts[5])  # new heme ID
                    heme_type = parts[6]  # 'c'
                    ligand_type = parts[7]  # 'HX'
                    is_c_type = True
                else:
                    continue
            
                # Check if this heme is in the requested list (using new heme ID)
                if heme_id in heme_ids:
                    # Map ligand type
#                   ligand_map = {
#                       'HH': 'HIS', 'HM': 'MET', 'HC': 'CYS',
#                       'HY': 'TYR', 'HD': 'ASP', 'HE': 'GLU',
#                       'HN': 'ASN', 'HQ': 'GLN', 'HK': 'LYS'
#                   }
#                   distal_type = ligand_map.get(ligand_type, 'HIS')

                    ligand_map = {
                        'HH': 'HH',  # Keep short code instead of converting to 'HIS'
                        'HM': 'HM',  # Keep 'HM' instead of 'MET'
                        'HC': 'HC',  # Keep 'HC' instead of 'CYS'
                        'HY': 'HY',  # Keep 'HY' instead of 'TYR'
                        'HD': 'HD',  # Keep 'HD' instead of 'ASP'
                        'HE': 'HE',  # Keep 'HE' instead of 'GLU'
                        'HN': 'HN',  # Keep 'HN' instead of 'ASN'
                        'HQ': 'HQ',  # Keep 'HQ' instead of 'GLN'
                        'HK': 'HK'   # Keep 'HK' instead of 'LYS'
                    }
                    distal_type = ligand_map.get(ligand_type, 'HH')  # Default to HH if unknown

                    heme_config = {
                        'heme_id': heme_id,
                        'original_heme_id': orig_heme_id,
                        'type': 'c' if is_c_type else 'b',
                        'his_p': his_p,
                        'distal': distal,
                        'distal_type': distal_type
                    }
                
                    # Add c-type specific info
                    if is_c_type:
                        heme_config['cys_b'] = cys_b
                        heme_config['cys_c'] = cys_c
                
                    heme_configs.append(heme_config)
                    found_heme_ids.add(heme_id)
    
        # Check if all requested heme IDs were found
        missing_hemes = set(heme_ids) - found_heme_ids
        if missing_hemes:
            raise ValueError(f"Could not find heme IDs in ResIndexing.txt: {missing_hemes}")
    
        return heme_configs

    def analyze_heme_types(self, structure) -> Dict[str, List[int]]:
        """
        Analyze a structure and return information about unique heme types present.
        
        Args:
            structure: BioPython Structure object to analyze
            
        Returns:
            Dict mapping heme residue names to lists of residue IDs
            Example: {
                'HCO': [42, 55],  # Oxidized His-His c-type hemes
                'MCR': [68],      # Reduced His-Met c-type heme
                'HBO': [73]       # Oxidized His-His b-type heme
            }
        """
        heme_types = {}
        
        for model in structure:
            for chain in model:
                for residue in chain:
                    res_name = residue.resname
                    res_id = residue.get_id()[1]  # Get residue number
                    
                    # Check if it's any type of heme
                    is_heme = False
                    for ligand_types in [self.BTYPE_HEME_NAMES, self.CTYPE_HEME_NAMES]:
                        for ligand_variants in ligand_types.values():
                            if res_name in ligand_variants.values():
                                is_heme = True
                                if res_name not in heme_types:
                                    heme_types[res_name] = []
                                heme_types[res_name].append(res_id)
                                break
                        if is_heme:
                            break
                        
        return heme_types

    def generate_reference_state(self, state: str = 'ox'):
        """Generate a reference state PDB with all hemes in specified redox state."""
        working_structure = copy.deepcopy(self.structure)

        for model in working_structure:
            for chain in model:
                for residue in list(chain):
                    res_name = residue.resname

                    # Check and rename b-type hemes and ligands
                    for ligand_type in self.BTYPE_HEME_NAMES:
                        if res_name in self.BTYPE_HEME_NAMES[ligand_type].values():
                            residue.resname = self.BTYPE_HEME_NAMES[ligand_type][state]
                            break
                        elif res_name in self.BTYPE_PROX_HIS[ligand_type].values():
                            residue.resname = self.BTYPE_PROX_HIS[ligand_type][state]
                            break
                        elif res_name in self.BTYPE_DISTAL[ligand_type].values():
                            residue.resname = self.BTYPE_DISTAL[ligand_type][state]
                            break

                    # Check and rename c-type hemes and ligands
                    for ligand_type in self.CTYPE_HEME_NAMES:
                        if res_name in self.CTYPE_HEME_NAMES[ligand_type].values():
                            residue.resname = self.CTYPE_HEME_NAMES[ligand_type][state]
                            break
                        elif res_name in self.CTYPE_PROX_HIS[ligand_type].values():
                            residue.resname = self.CTYPE_PROX_HIS[ligand_type][state]
                            break
                        elif res_name in self.CTYPE_DISTAL[ligand_type].values():
                            residue.resname = self.CTYPE_DISTAL[ligand_type][state]
                            break

        # Write output PDB
        output_path = self.launch_dir / "EE" / f"RefState_{state}.pdb"
        io = PDBIO()
        io.set_structure(working_structure)
        io.save(str(output_path), self.RefStateSelect())

        # Add TER records
        self.add_ter_records_to_pdb(output_path)

        return output_path

    def process_specific_heme_states(
        self, 
        heme_id: int, 
        is_c_type: bool, 
        his_p: int, 
        distal_ligand: int, 
        distal_type: str,  # Will be in short form (HH, HM, etc.)
        cys_b: Optional[int] = None, 
        cys_c: Optional[int] = None
    ) -> Tuple[str, str]:
        """
        Process both oxidized and reduced states for a specific heme.
        All other hemes remain in their original states.
        
        Args:
            heme_id: Heme residue ID
            is_c_type: Whether this is a c-type or b-type heme
            his_p: Proximal histidine residue ID
            distal_ligand: Distal ligand residue ID
            distal_type: Type of distal ligand in short form (e.g., 'HH', 'HM')
            cys_b: Optional cysteine B for c-type hemes
            cys_c: Optional cysteine C for c-type hemes
        
        Returns:
            Tuple of (oxidized_pdb_path, reduced_pdb_path)
        """
        # Convert short code to full name for dictionary lookups
        full_distal_type = self.SHORT_TO_FULL[distal_type]
        
        # Create deep copies to avoid modifying original structure
        working_structure_ox = copy.deepcopy(self.structure)
        working_structure_red = copy.deepcopy(self.structure)
        
        # Determine which nomenclature to use
        if is_c_type:
            heme_names = self.CTYPE_HEME_NAMES
            prox_his_names = self.CTYPE_PROX_HIS
            distal_names = self.CTYPE_DISTAL
        else:
            heme_names = self.BTYPE_HEME_NAMES
            prox_his_names = self.BTYPE_PROX_HIS
            distal_names = self.BTYPE_DISTAL
        
        # Prepare rename mappings for oxidized and reduced states
        ox_renames = {
            heme_id: heme_names[full_distal_type]['ox'],
            his_p: prox_his_names[full_distal_type]['ox'],
            distal_ligand: distal_names[full_distal_type]['ox']
        }
        
        red_renames = {
            heme_id: heme_names[full_distal_type]['red'],
            his_p: prox_his_names[full_distal_type]['red'],
            distal_ligand: distal_names[full_distal_type]['red']
        }
        
        # Handle c-type heme cysteines
        if is_c_type and cys_b and cys_c:
            ox_renames[cys_b] = 'CYO'
            ox_renames[cys_c] = 'CYO'
            red_renames[cys_b] = 'CYO'
            red_renames[cys_c] = 'CYO'
        
        # Apply renames to oxidized state
        for model in working_structure_ox:
            for chain in model:
                for residue in list(chain):
                    if residue.get_id()[1] in ox_renames:
                        residue.resname = ox_renames[residue.get_id()[1]]
        
        # Apply renames to reduced state
        for model in working_structure_red:
            for chain in model:
                for residue in list(chain):
                    if residue.get_id()[1] in red_renames:
                        residue.resname = red_renames[residue.get_id()[1]]
        
        # Write state PDBs
        ox_output_path = self.launch_dir / "EE" / f"o_{heme_id}_{heme_id}.pdb"
        red_output_path = self.launch_dir / "EE" / f"r_{heme_id}_{heme_id}.pdb"
        
        io_ox = PDBIO()
        io_ox.set_structure(working_structure_ox)
        io_ox.save(str(ox_output_path), self.RefStateSelect())
        
        io_red = PDBIO()
        io_red.set_structure(working_structure_red)
        io_red.save(str(red_output_path), self.RefStateSelect())
        
        # Add TER records
        self.add_ter_records_to_pdb(ox_output_path)
        self.add_ter_records_to_pdb(red_output_path)
    
        return ox_output_path, red_output_path

    def generate_heme_pair_microstates(
        self,
        heme1_id: int,
        is_heme1_c_type: bool,
        heme1_his_p: int,
        heme1_distal_ligand: int,
        heme1_distal_type: str,  # Will be in short form (HH, HM, etc.)
        heme2_id: int,
        is_heme2_c_type: bool,
        heme2_his_p: int,
        heme2_distal_ligand: int,
        heme2_distal_type: str,  # Will be in short form (HH, HM, etc.)
        heme1_cys_b: Optional[int] = None,
        heme1_cys_c: Optional[int] = None,
        heme2_cys_b: Optional[int] = None,
        heme2_cys_c: Optional[int] = None
    ) -> Dict[str, str]:
        """
        Generate all four redox microstates for a pair of hemes.

        Args:
            heme1_* and heme2_*: Configuration details for each heme
            distal_type args should be in short form (HH, HM, etc.)

        Returns:
            Dictionary of microstate names to their corresponding PDB paths
            Microstate names: 'RR', 'RO', 'OR', 'OO'
            Where R = Reduced, O = Oxidized
        """
        def generate_microstate(
            heme1_redox: str,
            heme2_redox: str
        ) -> str:
            # Convert short codes to full names for dictionary lookups
            heme1_full_type = self.SHORT_TO_FULL[heme1_distal_type]
            heme2_full_type = self.SHORT_TO_FULL[heme2_distal_type]
            
            # Create a deep copy of the structure
            working_structure = copy.deepcopy(self.structure)

            # Process each heme
            for heme_config in [
                {
                    'heme_id': heme1_id,
                    'type': 'c' if is_heme1_c_type else 'b',
                    'his_p': heme1_his_p,
                    'distal': heme1_distal_ligand,
                    'distal_type': heme1_full_type,
                    'cys_b': heme1_cys_b,
                    'cys_c': heme1_cys_c,
                    'redox': heme1_redox
                },
                {
                    'heme_id': heme2_id,
                    'type': 'c' if is_heme2_c_type else 'b',
                    'his_p': heme2_his_p,
                    'distal': heme2_distal_ligand,
                    'distal_type': heme2_full_type,
                    'cys_b': heme2_cys_b,
                    'cys_c': heme2_cys_c,
                    'redox': heme2_redox
                }
            ]:
                # Determine which nomenclature to use
                if heme_config['type'] == 'c':
                    heme_names = self.CTYPE_HEME_NAMES
                    prox_his_names = self.CTYPE_PROX_HIS
                    distal_names = self.CTYPE_DISTAL
                else:
                    heme_names = self.BTYPE_HEME_NAMES
                    prox_his_names = self.BTYPE_PROX_HIS
                    distal_names = self.BTYPE_DISTAL

                # Prepare rename mappings
                renames = {
                    heme_config['heme_id']: heme_names[heme_config['distal_type']][heme_config['redox']],
                    heme_config['his_p']: prox_his_names[heme_config['distal_type']][heme_config['redox']],
                    heme_config['distal']: distal_names[heme_config['distal_type']][heme_config['redox']]
                }

                # Handle c-type heme cysteines
                if heme_config['type'] == 'c':
                    if heme_config['cys_b']:
                        renames[heme_config['cys_b']] = 'CYO'
                    if heme_config['cys_c']:
                        renames[heme_config['cys_c']] = 'CYO'

                # Apply renames
                for model in working_structure:
                    for chain in model:
                        for residue in list(chain):
                            if residue.get_id()[1] in renames:
                                residue.resname = renames[residue.get_id()[1]]

            # Generate output filename
            output_filename = self.launch_dir / "EE" / f"{heme1_redox[0]}{heme2_redox[0]}_{heme1_id}_{heme2_id}.pdb"
            
            # Write PDB
            io = PDBIO()
            io.set_structure(working_structure)
            io.save(str(output_filename), self.RefStateSelect())
            
            # Add TER records
            self.add_ter_records_to_pdb(output_filename)

            return output_filename

        # Generate all four microstates
        return {
            'RR': generate_microstate('red', 'red'),  # Both reduced
            'RO': generate_microstate('red', 'ox'),   # Heme1 reduced, Heme2 oxidized
            'OR': generate_microstate('ox', 'red'),   # Heme1 oxidized, Heme2 reduced
            'OO': generate_microstate('ox', 'ox')     # Both oxidized
        }    

class RedoxStateManager:
    """Manages the generation and tracking of all redox states for a protein."""
    
    def __init__(self, 
                 input_pdb: str, 
                 heme_ids: List[int], 
                 launch_dir: Path,
                 forcefield_dir: Path,
                 reference_state: RedoxState):
        """
        Initialize RedoxStateManager.
        
        Args:
            input_pdb: Path to input PDB file
            heme_ids: List of heme residue IDs to process
            launch_dir: Path to launch directory
            forcefield_dir: Path to forcefield parameters directory
            reference_state: Redox state for reference PDB (OXIDIZED or REDUCED)
        """
        self.launch_dir = Path(launch_dir)
        self.ee_dir = self.launch_dir / "EE"
        self.ee_dir.mkdir(exist_ok=True)
        
        self.heme_generator = HemeStateGenerator(input_pdb, launch_dir)
        self.heme_ids = heme_ids
        self.reference_state = reference_state
        
        # Initialize forcefield handlers
        self.ff_params = ForceFieldParameters(forcefield_dir)
        self.tleap_gen = TLeapScriptGenerator(self.ff_params, launch_dir)
        
        # Generate reference state first
        self.reference_pdb = self._generate_reference_state()
        
        # After reference state is generated, reinitialize generator with it
        self.heme_generator = HemeStateGenerator(self.reference_pdb, launch_dir)
        
        # Track generated files
        self.generated_pdbs = {}
        self.generated_files = {}


    def _get_structure_info(self, pdb_path) -> Dict[str, Dict[str, int]]:
        """Analyze PDB to get heme types and counts"""
        structure = self.heme_generator.parser.get_structure('ref', pdb_path)
        return self.heme_generator.analyze_heme_types(structure)

    def generate_all_states(self):
        """
        Generate all required states:
        1. Reference state with all hemes in specified redox state
        2. Single-heme states (oxidized and reduced) for each heme
        3. Pair microstates for all possible heme pairs
        """

        print("\n=== Redox State Generation Workflow ===")
        print(f" Reference State: {self.reference_state.value.upper()}")
        print(f" Total Hemes to Process: {len(self.heme_ids)}")
        print(f" Heme IDs: {self.heme_ids}\n")
    
        # Reference state was already generated in __init__
        print("1. Generating Reference State PDB...")
        print(f"   Reference PDB: {self.reference_pdb}")

        
        # Process single heme states
        print("\n2. Generating Single Heme States...")
        self._generate_single_heme_states()
        
        # Process heme pairs
        print("\n3. Generating Heme Pair Microstates...")
        unique_pairs = [(i, j) for i in self.heme_ids for j in self.heme_ids if i < j]
        total_pairs = len(unique_pairs)
        print(f"   Total Number of Unique Heme Pair Combinations: {total_pairs}: {unique_pairs}")
        self._generate_pair_states()
        
        # Print summary
        print("\n4. Generating Summary...")
        self._print_generation_summary()

        print("\n=== Redox State Generation Complete ===")

    def _generate_reference_state(self) -> str:
        """Generate reference state PDB with all hemes in specified state."""
        ref_pdb = self.heme_generator.generate_reference_state(
            state=self.reference_state.value
        )
        return ref_pdb

    def _generate_single_heme_states(self):
        """Generate oxidized and reduced states for each heme individually."""
        logger.info("Generating single heme states...")
        
        for heme_id in self.heme_ids:
            logger.info(f"Processing heme {heme_id}")
            
            # Get heme configuration
            heme_config = self.heme_generator.parse_heme_environments([heme_id])[0]
            
            # Generate both oxidized and reduced PDBs
            ox_pdb, red_pdb = self.heme_generator.process_specific_heme_states(
                heme_id=heme_id,
                is_c_type=(heme_config['type'] == 'c'),
                his_p=heme_config['his_p'],
                distal_ligand=heme_config['distal'],
                distal_type=heme_config['distal_type'],
                cys_b=heme_config.get('cys_b'),
                cys_c=heme_config.get('cys_c')
            )
            
            # Track generated PDBs
            self.generated_pdbs[f"ox_{heme_id}"] = ox_pdb
            self.generated_pdbs[f"red_{heme_id}"] = red_pdb
            
            # Generate and run TLeaP for this heme
            try:
                self._process_single_heme_tleap(heme_id, heme_config)
            except Exception as e:
                logger.error(f"Failed to process TLeaP for heme {heme_id}: {str(e)}")
                raise

    def _generate_pair_states(self):
        """Generate all microstate combinations for each possible heme pair."""
        logger.info("Generating heme pair microstates...")
        
        # Generate all possible pairs
        pairs = [(i, j) for i in self.heme_ids for j in self.heme_ids if i < j]
        
        for heme1_id, heme2_id in pairs:
            logger.info(f"Processing heme pair {heme1_id}-{heme2_id}")
            
            # Get configurations for both hemes
            heme1_config = self.heme_generator.parse_heme_environments([heme1_id])[0]
            heme2_config = self.heme_generator.parse_heme_environments([heme2_id])[0]
            
            # Generate all microstates
            microstates = self.heme_generator.generate_heme_pair_microstates(
                heme1_id=heme1_id,
                is_heme1_c_type=(heme1_config['type'] == 'c'),
                heme1_his_p=heme1_config['his_p'],
                heme1_distal_ligand=heme1_config['distal'],
                heme1_distal_type=heme1_config['distal_type'],
                heme1_cys_b=heme1_config.get('cys_b'),
                heme1_cys_c=heme1_config.get('cys_c'),
                
                heme2_id=heme2_id,
                is_heme2_c_type=(heme2_config['type'] == 'c'),
                heme2_his_p=heme2_config['his_p'],
                heme2_distal_ligand=heme2_config['distal'],
                heme2_distal_type=heme2_config['distal_type'],
                heme2_cys_b=heme2_config.get('cys_b'),
                heme2_cys_c=heme2_config.get('cys_c')
            )
            
            # Track generated PDBs
            for state, pdb in microstates.items():
                self.generated_pdbs[f"{state.lower()}_{heme1_id}_{heme2_id}"] = pdb
            
            # Generate and run TLeaP for this pair
            try:
                self._process_pair_tleap(heme1_id, heme2_id, heme1_config, heme2_config)
            except Exception as e:
                logger.error(f"Failed to process TLeaP for heme pair {heme1_id}-{heme2_id}: {str(e)}")
                raise

    def _process_single_heme_tleap(self, heme_id: int, heme_config: Dict):
        """
        Generate and run TLeaP for a single selected heme.
        
        Args:
            heme_id: ID of selected heme
            heme_config: Configuration dictionary for selected heme
        """

#       print("\nProcessing TLeap for heme:", heme_id)
#       print("All heme IDs:", self.heme_ids)

        reference_hemes = {}
        ref_structure_info = self._get_structure_info(self.reference_pdb)
        
        # Load the structure first instead of passing the path
        ref_structure = self.heme_generator.parser.get_structure('ref', self.reference_pdb)
        all_heme_types = self.heme_generator.analyze_heme_types(ref_structure)
        
        for hid in self.heme_ids:
            if hid != heme_id:  # Skip the selected heme
                config = self.heme_generator.parse_heme_environments([hid])[0]
                reference_hemes[hid] = config
        
        # Generate TLeaP input file
        tleap_file = self.ee_dir / f"tleap_{heme_id}.in"
        success = self.tleap_gen.generate_single_heme_script(
            reference_hemes=reference_hemes,
            selected_heme=heme_config,
            output_file=tleap_file,
            reference_state=self.reference_state
        )
        
        if not success:
            raise RuntimeError(f"Failed to generate TLeaP script for heme {heme_id}")
        
        # Run TLeaP
        self._run_tleap(tleap_file)
        
        # Track generated files
        self.generated_files[f"ox_{heme_id}"] = {
            'prmtop': self.ee_dir / f"o_{heme_id}_{heme_id}.prmtop",
            'rst7': self.ee_dir / f"o_{heme_id}_{heme_id}.rst7"
        }
        self.generated_files[f"red_{heme_id}"] = {
            'prmtop': self.ee_dir / f"r_{heme_id}_{heme_id}.prmtop",
            'rst7': self.ee_dir / f"r_{heme_id}_{heme_id}.rst7"
        }

    def _process_pair_tleap(self, heme1_id: int, heme2_id: int, 
                           heme1_config: Dict, heme2_config: Dict):
        """
        Generate and run TLeaP for a pair of selected hemes.
        
        Args:
            heme1_id: ID of first selected heme
            heme2_id: ID of second selected heme
            heme1_config: Configuration dictionary for first heme
            heme2_config: Configuration dictionary for second heme
        """
        # Get reference heme configurations (all hemes except selected pair)
        reference_hemes = {}
        # Parse the reference PDB path into a structure
        ref_structure = self.heme_generator.parser.get_structure('ref', str(self.reference_pdb))
        all_heme_types = self.heme_generator.analyze_heme_types(ref_structure)
        
        for hid in self.heme_ids:
            if hid not in (heme1_id, heme2_id):  # Skip the selected hemes
                config = self.heme_generator.parse_heme_environments([hid])[0]
                reference_hemes[hid] = config
        
        # Generate TLeaP input file
        tleap_file = self.ee_dir / f"tleap_{heme1_id}_{heme2_id}.in"
        success = self.tleap_gen.generate_pair_script(
            reference_hemes=reference_hemes,
            heme1_config=heme1_config,
            heme2_config=heme2_config,
            output_file=tleap_file,
            reference_state=self.reference_state
        )
        
        if not success:
            raise RuntimeError(f"Failed to generate TLeaP script for heme pair {heme1_id}-{heme2_id}")
        
        # Run TLeaP
        self._run_tleap(tleap_file)
        
        # Track generated files
        states = ['oo', 'or', 'ro', 'rr']
        for state in states:
            state_key = f"{state.lower()}_{heme1_id}_{heme2_id}"
            self.generated_files[state_key] = {
                'prmtop': self.ee_dir / f"{state}_{heme1_id}_{heme2_id}.prmtop",
                'rst7': self.ee_dir / f"{state}_{heme1_id}_{heme2_id}.rst7"
            }

    def _run_tleap(self, input_file: Path):
        """
        Run TLeaP with given input file.
        
        Args:
            input_file: Path to TLeaP input file
        """
        log_file = input_file.with_suffix('.log')
        try:
            print(f"tleap -s -f {input_file} > {log_file}")
            subprocess.run(
                f"tleap -s -f {input_file} > {log_file}",
                shell=True,
                check=True
            )
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"TLeaP failed for {input_file}. Check {log_file} for details.") from e

        # Check for successful completion by looking at the last line
        try:
            with open(log_file, 'r') as f:
                # Read the last line
                last_line = f.readlines()[-1].strip()

                # Check if the last line indicates zero errors
                if not last_line.startswith('Exiting LEaP:') or 'Errors = 0' not in last_line:
                    raise RuntimeError(f"TLEaP did not complete successfully. Check {log_file} for details.")
                else:
                    print("TLEaP successfully completed.")
        except (IndexError, IOError) as e:
            raise RuntimeError(f"Could not read log file {log_file}: {str(e)}")
       
    def _print_generation_summary(self):
        """
        Print a summary of generated files and PDBs.
        """
        # Categorize PDBs
        pdb_categories = {
            'Single Heme States': [name for name in self.generated_pdbs.keys() if '_' in name],
            'Microstate Pairs': [name for name in self.generated_pdbs.keys() if len(name.split('_')) == 1]
        }

        print("\n--- Detailed Generation Summary ---")
    
        # Summary statistics
        print(f"Total Generated PDBs: {len(self.generated_pdbs)}")
        print(f"Generated PDBs:")
        for category, pdbs in pdb_categories.items():
            print(f"\n{category} PDBs:")
            for name in pdbs:
                print(f"  {name}: {self.generated_pdbs[name]}")
    
        # Print generated topology and restart files
        print(f"Total Generated Topology/Restart Files: {len(self.generated_files)}")
        print("\nGenerated Topology and Restart Files:")
        for name, files in self.generated_files.items():
            print(f"  {name}:")
            for file_type, path in files.items():
                print(f"    {file_type}: {path}")

# Example usage
if __name__ == "__main__":
    # Setup logging
    logging.basicConfig(level=logging.INFO,
                       format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    
    # Example usage
    input_pdb = "input.pdb"
    heme_ids = [1280, 1274, 1277, 1271, 1268, 1283]
    launch_dir = Path("launch_dir")
    forcefield_dir = Path("forcefield_dir")
    
    try:
        manager = RedoxStateManager(
            input_pdb=input_pdb,
            heme_ids=heme_ids,
            launch_dir=launch_dir,
            forcefield_dir=forcefield_dir,
            reference_state=RedoxState.REDUCED
        )
        
        # Generate all states
        manager.generate_all_states()
        
    except Exception as e:
        logger.error(f"Error during state generation: {str(e)}")
        raise
