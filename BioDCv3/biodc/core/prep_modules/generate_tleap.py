# biodc/core/preparation/generate_tleap.py

"""
Module for generating AMBER/tleap input files.
Handles forcefield setup, parameter loading, and structure preparation.
Maintains full interactivity for structure preparation options while
supporting input file driven execution.
"""

import os
import sys
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass

from biodc.core.prep_modules import constant_ph_prep
from biodc.utils.interaction import InteractionManager

@dataclass
class TleapSettings:
    """Settings for tleap script generation."""
    out_prefix: str
    solvent_type: str  # 'explicit' or 'implicit'
    box_type: Optional[str] = None  # 'rectangular' or 'octahedral'
    buffer_size: Optional[float] = None
    na_count: int = 0
    cl_count: int = 0
    ff_choice: Optional[str] = None  # For alternative charge sets

@dataclass
class HemeAtomTypes:
    """Container for heme-specific atom types."""
    metal: str        # Fe atom type
    nitrogens: List[str]  # Porphyrin nitrogen types
    extra: Optional[str] = None  # Additional atom type (e.g., Met sulfur)
    description: str = ""

@dataclass
class IndexingData:
    """Container for data from ResIndexing.txt."""
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

class TleapGenerator:
    """Generates tleap input script based on structure analysis."""

    # Add this new mapping
    LIGAND_TYPE_MAP = {
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

    # Defined atom types for implemented heme types
    ATOMTYPES = {
        # b-type His-His oxidized
        ('HH-b', 'ox'): HemeAtomTypes(
            metal="M1",
            nitrogens=["Y1", "Y2", "Y3", "Y4", "Y5", "Y6"],
            description="Oxidized His-His b-type"
        ),
        # b-type His-His reduced
        ('HH-b', 'red'): HemeAtomTypes(
            metal="M2",
            nitrogens=["Z1", "Z2", "Z3", "Z4", "Z5", "Z6"],
            description="Reduced His-His b-type"
        ),
        # b-type His-Met oxidized
        ('HM-b', 'ox'): HemeAtomTypes(
            metal="M3",
            nitrogens=["W1", "W2", "W3", "W4", "W5", "W6"],
            extra="S1",  # Met sulfur
            description="Oxidized His-Met b-type"
        ),
        # b-type His-Met reduced
        ('HM-b', 'red'): HemeAtomTypes(
            metal="M4",
            nitrogens=["X1", "X2", "X3", "X4", "X5", "X6"],
            extra="S2",  # Met sulfur
            description="Reduced His-Met b-type"
        ),
        # c-type His-His oxidized
        ('HH-c', 'ox'): HemeAtomTypes(
            metal="M7",
            nitrogens=["S1", "S2", "S3", "S4", "S5", "S6"],
            description="Oxidized His-His c-type"
        ),
        # c-type His-His reduced
        ('HH-c', 'red'): HemeAtomTypes(
            metal="M8",
            nitrogens=["T1", "T2", "T3", "T4", "T5", "T6"],
            description="Reduced His-His c-type"
        ),
        # c-type His-Met oxidized
        ('HM-c', 'ox'): HemeAtomTypes(
            metal="M5",
            nitrogens=["U1", "U2", "U3", "U4", "U5", "U6"],
            extra="S3",  # Met sulfur
            description="Oxidized His-Met c-type"
        ),
        # c-type His-Met reduced
        ('HM-c', 'red'): HemeAtomTypes(
            metal="M6",
            nitrogens=["V1", "V2", "V3", "V4", "V5", "V6"],
            extra="S4",  # Met sulfur
            description="Reduced His-Met c-type"
        )
    }
    
    # Define bonding patterns for ligands
    LIGAND_BONDS = {
        'HIS': 'NE2',  # His coordinates through NE2
        'MET': 'SD',   # Met coordinates through sulfur
        'CYS': 'SG',   # Cys coordinates through sulfur
        'TYR': 'OH',   # Tyr coordinates through oxygen
        'ASP': 'OD1',  # Asp coordinates through carboxylate O
        'GLU': 'OE1',  # Glu coordinates through carboxylate O
        'ASN': 'OD1',  # Asn coordinates through O
        'GLN': 'OE1',  # Gln coordinates through O
        'LYS': 'NZ'    # Lys coordinates through terminal N
    }

    def __init__(self, 
                 forcefield_dir: Path,
                 structure_info: Dict[str, Dict[str, int]],
                 indexing_data: List[IndexingData],
                 launch_dir: Path):
        """
        Initialize tleap generator.
        
        Args:
            forcefield_dir: Directory containing forcefield files
            structure_info: Output from PDBProcessor.count_heme_types()
            indexing_data: Parsed data from ResIndexing.txt
            launch_dir: Directory for output files
        """
        self.forcefield_dir = Path(forcefield_dir)
        self.structure_info = structure_info
        self.indexing_data = indexing_data
        self.launch_dir = launch_dir
        self.settings = None  # Will be set by get_settings

    def get_settings(self, input_dict: Dict) -> TleapSettings:
        interaction_manager = InteractionManager(
            launch_dir=self.launch_dir, 
            input_dict=input_dict
        )

        print("\n")
        print("=" * 60)
        print("""TLEaP Setup:

Now, we submit the processed PDB with appropriately renamed atoms/residues
and renumbered residues to TLEaP of the AmberTools package to generate 
topology and coordinate files.""")

        out_prefix = interaction_manager.prompt(
            "OutPrefix", 
            "\nPrefix for output parm/rst7"
        )

        solvent_type = interaction_manager.prompt(
            "SolvEnv", 
            "\nShould the structure be prepared with an explicit or implicit solventV?", 
            choices=['exp', 'imp']
        )

        settings = TleapSettings(
            out_prefix=out_prefix,
            solvent_type=solvent_type
        )

        if solvent_type == 'exp':
            settings.box_type = interaction_manager.prompt(
                "BoxShape", 
                "Using a rectangular or an octahedral box?", 
                choices=['rec', 'octahed']
            )

            settings.buffer_size = interaction_manager.prompt(
                "BufferSize", 
                "With how much of a solvent buffer (in angstroms)?", 
                input_type=float
            )

            settings.na_count = interaction_manager.prompt(
                "NaCount", 
                "How many Na+ ions? (0 = enough for charge neutrality)", 
                input_type=int
            )

            settings.cl_count = interaction_manager.prompt(
                "ClCount", 
                "How many Cl- ions? (0 = enough for charge neutrality)", 
                input_type=int
            )

        if any(htype == 'HH-c' and counts['ox'] > 0 
            for htype, counts in self.structure_info.items()):
            settings.ff_choice = interaction_manager.yes_no_prompt(
                "FFchoice", 
                """
    Alternative charges are available for oxidized His-His c-type hemes:
    [1] Henriques, J.; Costa, P. J.; Calhorda, M. J.; 
        Machuqueiro, M. Charge Parametrization of the
        DvH-c3 Heme Group: Validation Using Constant-(pH,E) 
        Molecular Dynamics Simulations. J. Phys. Chem. B 
        2013, 117 (1), 70–82.

    Would you like to use these alternative charges?"""
            )

        self.settings = settings
        return settings

    def generate_script(self) -> str:
        """Generate complete tleap input script."""
        if not self.settings:
            raise ValueError("Must call get_settings before generating script")
            
        script = []
        
        # Add header and basic forcefield loading
        script.extend([
            "# Load parameters",
            "source leaprc.constph",
            "source leaprc.conste",
            "source leaprc.gaff",
            "source leaprc.water.tip3p",
            ""
        ])
        
        # Add atom types
        script.extend(self._generate_atomtypes())
        
        # Add references for implemented forcefields
        script.extend(self._generate_references())
        
        # Load heme parameters
        script.extend(self._generate_parameter_loading())
        
        # Load structure
        script.extend([
            f"\n# Load PDB",
            f"{self.settings.out_prefix} = loadpdb processed.pdb",
            ""
        ])
        
        # Add all bond definitions
        script.extend(self._generate_bonds())
        
        # Add solvent if requested
        if self.settings.solvent_type == 'exp':
            script.extend(self._generate_solvation())
        
        # Save outputs
        script.extend([
            "\n# Save topology and coordinate files",
            f"saveamberparm {self.settings.out_prefix} "
            f"{self.settings.out_prefix}.prmtop {self.settings.out_prefix}.rst7",
            "\nquit"
        ])
        
        return "\n".join(script)
        
    def _generate_atomtypes(self) -> List[str]:
        """Generate atom type definitions based on structure content."""
        print("Structure info in _generate_atomtypes:", self.structure_info) # Before first line

        script = ["\naddAtomTypes {"]
        
        for (heme_type, redox), count in self.structure_info.items():
            if count > 0 and (heme_type, redox) in self.ATOMTYPES:
                types = self.ATOMTYPES[(heme_type, redox)]
                # Add metal
                script.append(f'    {{ "{types.metal}"  "Fe" "sp3" }}  # {types.description}')
                # Add nitrogens
                for n_type in types.nitrogens:
                    script.append(f'    {{ "{n_type}"  "N" "sp3" }}')
                # Add extra types if any
                if types.extra:
                    script.append(f'    {{ "{types.extra}"  "S" "sp3" }}')
                
        script.append("}")
        return script
        
    def _generate_references(self) -> List[str]:
        """Generate reference documentation in script."""
        script = []
        
        # Add b-type references if needed
        if any(htype.endswith('-b') for htype, _ in self.structure_info):
            script.extend([
                "\n# References for b-type heme forcefield parameters:",
                "#    Bonded parameters for the macrocycle come from:",
                "#      Yang, Longhua, Åge A. Skjevik, Wen-Ge Han Du, Louis Noodleman,",
                "#      Ross C. Walker, and Andreas W. Götz. Data for molecular",
                "#      dynamics simulations of B-type cytochrome c oxidase with",
                "#      the Amber force field. Data in brief 8 (2016): 1209-1214.",
                "#",
                "#    Bonded parameters for the Fe center and atomic partial charges were derived",
                "#    by Guberman-Pfeffer using the Metal Center Parameter Builder.",
                "#    The B3LYP approximate density functional was used with the mixed basis set",
                "#    (LANL2TZ(f) for Fe and 6-31G(d) for 2nd row elements."
            ])

        # Add c-type references if needed
        if any(htype.endswith('-c') for htype, _ in self.structure_info):
            script.extend([
                "\n# References for c-type heme forcefield parameters:",
                "#    Bonded parameters for the macrocycle come from:",
                "#      Crespo, A.; Martí, M. A.; Kalko, S. G.; Morreale, A.; Orozco, M.;",
                "#      Gelpi, J. L.; Luque, F. J.; Estrin, D. A. Theoretical Study of the",
                "#      Truncated Hemoglobin HbN: Exploring the Molecular Basis of the NO",
                "#      Detoxification Mechanism. J. Am. Chem. Soc. 2005, 127 (12), 4433–4444."
            ])

            if any(htype == 'HH-c' for htype, _ in self.structure_info):
                script.extend([
                    "#",
                    "#    Alternative charges from:",
                    "#      Henriques, J.; Costa, P. J.; Calhorda, M. J.; Machuqueiro, M.",
                    "#      Charge Parametrization of the DvH-c3 Heme Group: Validation Using",
                    "#      Constant-(pH,E) Molecular Dynamics Simulations.",
                    "#      J. Phys. Chem. B 2013, 117 (1), 70–82."
                ])

        return script


    def _generate_parameter_loading(self) -> List[str]:
        """Generate parameter loading commands."""
        script = []
        
        for (heme_type, redox), count in self.structure_info.items():
            if count > 0:
                redox_prefix = "Oxidized" if redox == "ox" else "Reduced"
                
                # Handle His-His c-type with alternative charges
                if (heme_type == 'HH-c' and redox == 'ox' and 
                    self.settings.ff_choice and 
                    self.settings.ff_choice.lower() in ('yes', 'y')):
                    script.extend([
                        f"\n# Load {redox_prefix} His-His c-type parameters with Henriques charges",
                        f"loadamberparams {self.forcefield_dir}/Oxidized_HisHisLigated_c-heme.frcmod",
                        f"loadoff {self.forcefield_dir}/Henriques_Oxidized_HisHisLigated_c-heme_RESP.lib"
                    ])
                else:
                    # Standard parameter loading
                    ligation = "HisHis" if heme_type.startswith("HH") else "HisMet"
                    heme_class = "b" if heme_type.endswith('-b') else "c"
                    script.extend([
                        f"\n# Load {redox_prefix} {ligation} {heme_class}-type parameters",
                        f"loadamberparams {self.forcefield_dir}/{redox_prefix}_{ligation}Ligated_{heme_class}-heme.frcmod",
                        f"loadoff {self.forcefield_dir}/{redox_prefix}_{ligation}Ligated_{heme_class}-heme_RESP.lib"
                    ])
        
        return script

    def _generate_bonds(self) -> List[str]:
        """Generate all bond definitions for the structure."""
        script = []
        defined_bonds = set()  # Track defined bonds to avoid duplicates

        # Handle disulfide bonds
        for data in self.indexing_data:
            if isinstance(data, dict) and data.get('type') == 'disulfide':
                script.append("\n#------------------------------------------------------------")
                script.append("#Disulfide bonds:")
                for cys1, cys2 in data['pairs']:
                    script.append(
                        f"bond {self.settings.out_prefix}.{cys1}.SG "
                        f"{self.settings.out_prefix}.{cys2}.SG"
                    )
                script.append("#------------------------------------------------------------\n")

        # Handle heme bonds
        for heme in self.indexing_data:
            if not isinstance(heme, dict) or heme.get('type') != 'disulfide':
                script.append(f"\n#------------------------------------------------------------")
                script.append(f"#For heme {heme.new_id} (originally {heme.heme_id}):")

                # Bond Fe to ligands
                script.append(f"\n#Bond ligating atoms to Fe center")
                script.append(
                    f"bond {self.settings.out_prefix}.{heme.his_p}.NE2 "
                    f"{self.settings.out_prefix}.{heme.new_id}.FE"
                )

                # Convert abbreviated code to full residue name
                full_type = self.LIGAND_TYPE_MAP.get(heme.distal_type)
                if not full_type:
                    raise ValueError(f"Unknown ligand type: {heme.distal_type}")

                bonding_atom = self.LIGAND_BONDS[full_type]
                script.append(
                    f"bond {self.settings.out_prefix}.{heme.distal_ligand}.{bonding_atom} "
                    f"{self.settings.out_prefix}.{heme.new_id}.FE"
                )
                
                # Add backbone connections
                script.append(f"\n#Bond axially coordinated residues to preceeding and proceeding residues")
                for ligand_id in [heme.his_p, heme.distal_ligand]:
                    for bond in [(ligand_id-1, ligand_id), (ligand_id, ligand_id+1)]:
                        if bond not in defined_bonds:
                            script.append(
                                f"bond {self.settings.out_prefix}.{bond[0]}.C "
                                f"{self.settings.out_prefix}.{bond[1]}.N"
                            )
                            defined_bonds.add(bond)
                
                # For c-type hemes, add Cys thioether bonds
                if heme.is_c_type and heme.cys_b and heme.cys_c:
                    script.append(f"\n#Bond heme thioethers to protein backbone")
                    script.append(
                        f"bond {self.settings.out_prefix}.{heme.cys_b}.CA "
                        f"{self.settings.out_prefix}.{heme.new_id}.CBB2"
                    )
                    script.append(
                        f"bond {self.settings.out_prefix}.{heme.cys_c}.CA "
                        f"{self.settings.out_prefix}.{heme.new_id}.CBC1"
                    )
                
                # Add propionate bonds
                script.append(f"\n#Bond propionic acids to heme")
                script.append(
                    f"bond {self.settings.out_prefix}.{heme.new_id}.C2A "
                    f"{self.settings.out_prefix}.{heme.new_id+1}.CA"
                )
                script.append(
                    f"bond {self.settings.out_prefix}.{heme.new_id}.C3D "
                    f"{self.settings.out_prefix}.{heme.new_id+2}.CA"
                )
                
                script.append(f"#------------------------------------------------------------")
        
        return script

    def _generate_solvation(self) -> List[str]:
        """Generate solvation and ion commands."""
        script = []
        
        # Add solvent box
        if self.settings.box_type == 'rec':
            script.append(f"\n#Solvate with rectangular box")
            script.append(f"solvateBox {self.settings.out_prefix} "
                        f"TIP3PBOX {self.settings.buffer_size}")
        else:  # octahedral
            script.append(f"\n#Solvate with octahedral box")
            script.append(f"solvateOct {self.settings.out_prefix} "
                        f"TIP3PBOX {self.settings.buffer_size}")
            
        # Add ions if requested
        if self.settings.na_count > 0 or self.settings.cl_count > 0:
            script.append("\n#Add ions")
            if self.settings.na_count > 0:
                script.append(f"addions {self.settings.out_prefix} Na+ {self.settings.na_count}")
            if self.settings.cl_count > 0:
                script.append(f"addions {self.settings.out_prefix} Cl- {self.settings.cl_count}")
            
        return script

def generate_tleap_input(pdb: str,
                        forcefield_dir: Path,
                        structure_info: Dict[str, Dict[str, int]],
                        indexing_data: List[Dict],
                        input_dict: Dict,
                        launch_dir: Path) -> Tuple[str, str]:
    """
    Main function for generating tleap input script.
    
    Args:
        pdb: PDB file prefix
        forcefield_dir: Directory containing forcefield files
        structure_info: Output from PDBProcessor.count_heme_types()
        indexing_data: Parsed data from ResIndexing.txt
        input_dict: User input dictionary
        launch_dir: Directory for output files
    
    Returns:
        Tuple of (output_prefix, solvent_type)
    """
    # Convert indexing data to proper format
    parsed_indexing = []
    for data in indexing_data:
        if isinstance(data, dict) and data.get('type') == 'disulfide':
            # Pass through disulfide data as-is
            parsed_indexing.append(data)
        else:
            # Convert heme data to IndexingData
            parsed_indexing.append(IndexingData(**data))

    # Initialize generator
    generator = TleapGenerator(forcefield_dir, structure_info, parsed_indexing, launch_dir)
    
    # Get settings interactively or from input dict
    settings = generator.get_settings(input_dict)
    
    # Generate script
    script = generator.generate_script()
    
    # Write script
    with open(launch_dir / "SPR" / "tleap.in", "w") as f:
        f.write(script)
    
    print("\nGenerated tleap.in script for structure preparation.")
    
    # Run tleap
    print("Running tleap...")
    subprocess.run(["tleap", "-s", "-f", "tleap.in"], 
                  stdout=open("tleap.log", "w"),
                  stderr=subprocess.STDOUT,
                  check=True)
    print("Tleap completed! Please check tleap.log for any errors.")
    
    return settings.out_prefix, settings.solvent_type
