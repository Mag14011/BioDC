"""
Module for calculating reorganization energy using SASA-based approach.
Supports multiple SASA calculation backends and handles complete heme environments.
"""

import os
import logging
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Protocol, Tuple, Set, Union
from dataclasses import dataclass
from abc import ABC, abstractmethod

from biodc.utils.interaction import InteractionManager

logger = logging.getLogger(__name__)

@dataclass
class SASAResult:
    """Container for SASA calculation results."""
    donor_sasa: float
    acceptor_sasa: float
    distance: float

@dataclass
class HemeEnvironment:
    """Container for heme and its coordinating residues."""
    heme_id: int
    his_p: int  # Proximal His
    distal_ligand: int
    distal_type: str  # H, M, C, etc.
    is_c_type: bool
    cys_b: Optional[int] = None  # Only for c-type
    cys_c: Optional[int] = None  # Only for c-type

    def get_all_residues(self) -> List[int]:
        """Get list of all residues in this heme's environment."""
        residues = [self.heme_id, self.his_p, self.distal_ligand]
        if self.is_c_type:
            residues.extend([self.cys_b, self.cys_c])
        return [r for r in residues if r is not None]

class ResidueIndexParser:
    """Parser for ResIndexing.txt file."""

    def __init__(self, launch_dir: Path):
        self.launch_dir = Path(launch_dir)
        self.spr_dir = self.launch_dir / "SPR"

    def find_indexing_file(self) -> Path:
        """Find ResIndexing.txt in either SPR or current directory."""
        possible_locations = [
            self.spr_dir / "ResIndexing.txt",
            Path.cwd() / "ResIndexing.txt"
        ]

        for loc in possible_locations:
            if loc.exists():
                logger.info(f"Found ResIndexing.txt at {loc}")
                return loc
        raise FileNotFoundError(
            "Could not find ResIndexing.txt in SPR or current directory"
        )

    def parse_file(self) -> Dict[int, HemeEnvironment]:
        """Parse ResIndexing.txt and return mapping of heme ID to environment."""
        index_file = self.find_indexing_file()
        environments = {}

        with open(index_file) as f:
            for line_num, line in enumerate(f, 1):
                parts = line.strip().split()
                if not parts:
                    continue

                try:
                    if len(parts) == 8:  # c-type
                        env = HemeEnvironment(
                            heme_id=int(parts[5]),  # Use new_HemeID
                            his_p=int(parts[2]),
                            distal_ligand=int(parts[3]),
                            distal_type=parts[7][1],  # Second char of HX
                            is_c_type=True,
                            cys_b=int(parts[0]),
                            cys_c=int(parts[1])
                        )
                    elif len(parts) == 6:  # b-type
                        env = HemeEnvironment(
                            heme_id=int(parts[3]),  # Use new_HemeID
                            his_p=int(parts[0]),
                            distal_ligand=int(parts[1]),
                            distal_type=parts[5][1],  # Second char of HX
                            is_c_type=False
                        )
                    else:
                        logger.warning(f"Skipping line {line_num}: Invalid format")
                        continue

                    environments[env.heme_id] = env
                    logger.debug(f"Parsed environment for heme {env.heme_id}")

                except (ValueError, IndexError) as e:
                    logger.warning(f"Error parsing line {line_num}: {line.strip()}: {e}")

        return environments

class SASACalculator(ABC):
    """Abstract base class for SASA calculation engines."""
    
    @abstractmethod
    def calculate(self, 
                 pdb_file: Union[str, Path],
                 donor_env: HemeEnvironment,
                 acceptor_env: HemeEnvironment) -> SASAResult:
        """Calculate SASA for donor and acceptor environments.
        
        Args:
            pdb_file: Path to PDB file
            donor_env: Donor heme environment
            acceptor_env: Acceptor heme environment
            
        Returns:
            SASAResult containing SASA values and Fe-Fe distance
            
        Raises:
            RuntimeError: If SASA calculation fails
            FileNotFoundError: If PDB file not found
        """
        pass

class VMDCalculator(SASACalculator):
    """VMD-based SASA calculator."""

    def __init__(self):
        """Initialize VMD calculator."""
        try:
            subprocess.run(['vmd', '-h'],
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
            self.vmd_available = True
        except FileNotFoundError:
            self.vmd_available = False
            logger.warning("VMD not found in system path")

    def generate_tcl_script(self,
                          pdb_file: str,
                          donor_env: HemeEnvironment,
                          acceptor_env: HemeEnvironment) -> str:
        """Generate TCL script for VMD SASA calculation."""

        donor_residues = donor_env.get_all_residues()
        acceptor_residues = acceptor_env.get_all_residues()

        script = f"""
mol new {pdb_file}
set allsel [atomselect top "all and not water and not ions"]

set donor [atomselect top "not water and not ions and resid {' '.join(map(str, donor_residues))} and not name N H CA HA C O"]
set acceptor [atomselect top "not water and not ions and resid {' '.join(map(str, acceptor_residues))} and not name N H CA HA C O"]

set dsasa [measure sasa 1.4 $allsel -restrict $donor]
set asasa [measure sasa 1.4 $allsel -restrict $acceptor]

set d_fe [atomselect top "resid {donor_env.heme_id} and name FE"]
set a_fe [atomselect top "resid {acceptor_env.heme_id} and name FE"]
set d_coords [lindex [$d_fe get {{x y z}}] 0]
set a_coords [lindex [$a_fe get {{x y z}}] 0]
set distance [veclength [vecsub $d_coords $a_coords]]

puts "Donor_SASA (Å^2) Acceptor_SASA (Å^2) Fe-Fe Distance (Å)"
puts " $dsasa $asasa $distance"

quit
"""
        return script

    def calculate(self,
                pdb_file: str,
                donor_env: HemeEnvironment,
                acceptor_env: HemeEnvironment) -> SASAResult:
        """Calculate SASA using VMD."""
        if not self.vmd_available:
            raise RuntimeError("VMD is not available")

        # Generate and save TCL script
        script = self.generate_tcl_script(pdb_file, donor_env, acceptor_env)
        script_file = "sasa_calc.tcl"
        with open(script_file, 'w') as f:
            f.write(script)

        # Run VMD
        result = subprocess.run(['vmd', '-dispdev', 'text', '-e', script_file],
                            capture_output=True,
                            text=True)

        # Parse output
        try:
            # Look for the line after "Donor_SASA"
            lines = result.stdout.split('\n')
            for i, line in enumerate(lines):
                if "Donor_SASA" in line and i + 1 < len(lines):
                    # Get the next line which contains our values
                    values = lines[i + 1].strip().split()
                    if len(values) == 3:
                        dsasa, asasa, distance = map(float, values)
                        # Only delete script after successful parsing
                        os.remove(script_file)
                        return SASAResult(dsasa, asasa, distance)
        
            raise ValueError("Could not find SASA values in VMD output")
            
        except (ValueError, IndexError) as e:
            # Keep script file for debugging if parsing fails
            raise RuntimeError(f"Failed to parse VMD output: {e}\nTCL script preserved in {script_file}")

#       return SASAResult(dsasa, asasa, distance)

class FreeSASACalculator(SASACalculator):
    """FreeSASA-based calculator."""

    def __init__(self):
        """Initialize FreeSASA calculator."""
        try:
            import freesasa
            self.freesasa = freesasa
        except ImportError:
            raise ImportError("FreeSASA not installed. Install with: pip install freesasa")

        # Set up parameters similar to VMD's defaults
        self.parameters = self.freesasa.Parameters({
            'algorithm': self.freesasa.ShrakeRupley,
            'probe-radius': 1.4,  # Å (same as VMD default)
            'n-points': 100,      # Resolution of calculation
            'n-threads': 4        # Use multiple threads for speed
        })

    def calculate(self,
                 pdb_file: str,
                 donor_env: HemeEnvironment,
                 acceptor_env: HemeEnvironment) -> SASAResult:
        """Calculate SASA using FreeSASA."""
        # Read structure and calculate total SASA
        structure = self.freesasa.Structure(pdb_file)
        result = self.freesasa.calc(structure, self.parameters)

        donor_residues = donor_env.get_all_residues()
        acceptor_residues = acceptor_env.get_all_residues()

        # Calculate SASA for selections
        donor_sasa = 0
        acceptor_sasa = 0
        fe_coords = {'donor': None, 'acceptor': None}

        for i in range(structure.nAtoms()):
            resid = structure.residueNumber(i)
            atom_name = structure.atomName(i)

            # Store Fe coordinates for distance calculation
            if atom_name == "FE":
                if resid == donor_env.heme_id:
                    fe_coords['donor'] = structure.coord(i)
                elif resid == acceptor_env.heme_id:
                    fe_coords['acceptor'] = structure.coord(i)

            # Sum up SASA for each atom in selections
            if (resid in donor_residues and
                atom_name not in ["N", "H", "CA", "HA", "C", "O"]):
                donor_sasa += result.atomArea(i)

            elif (resid in acceptor_residues and
                  atom_name not in ["N", "H", "CA", "HA", "C", "O"]):
                acceptor_sasa += result.atomArea(i)

        # Calculate Fe-Fe distance
        if None in fe_coords.values():
            raise ValueError("Could not find Fe atoms in both donor and acceptor")

        distance = self._calculate_distance(fe_coords['donor'], fe_coords['acceptor'])

        return SASAResult(donor_sasa, acceptor_sasa, distance)

    def _calculate_distance(self, coord1, coord2) -> float:
        """Calculate distance between two 3D coordinates."""
        import numpy as np
        return np.sqrt(sum((a - b) ** 2 for a, b in zip(coord1, coord2)))

@dataclass
class ReorganizationParameters:
    """Physical parameters for reorganization energy calculation."""
    alpha: float = 5.18    # Dielectric parameter
    beta: float = 0.016    # SASA scaling factor
    rd: float = 4.6       # Donor radius (Å)
    ra: float = 4.6       # Acceptor radius (Å)
    eopt: float = 1.84    # Optical dielectric constant

class LambdaCalculator:
    """Calculator for reorganization energy using SASA method."""
    
    def __init__(self, 
                 interaction_manager: InteractionManager,
                 launch_dir: Path,
                 parameters: Optional[ReorganizationParameters] = None,
                 sasa_backend: str = 'vmd'):
        """Initialize calculator.
        
        Args:
            interaction_manager: Manager for user interactions
            launch_dir: Project launch directory
            parameters: Optional custom parameters for calculation
            sasa_backend: SASA calculation method ('vmd' or 'freesasa')
        """
        self.interaction_manager = interaction_manager
        self.launch_dir = Path(launch_dir)
        self.ee_dir = self.launch_dir / "EE"
        self.params = parameters or ReorganizationParameters()
        
        # Initialize components
#       self.index_parser = ResidueIndexParser(launch_dir)
#       self.heme_environments = self.index_parser.parse_file()
        self._index_parser = None 
        self._heme_environments = None
        
        # Initialize SASA calculator
        if sasa_backend == 'vmd':
            self.sasa_calculator = VMDCalculator()
        elif sasa_backend == 'freesasa':
            self.sasa_calculator = FreeSASACalculator()
        else:
            raise ValueError(f"Unknown SASA backend: {sasa_backend}")

    def _load_heme_environments(self):
        """Lazily load heme environments only when needed."""
        if self._heme_environments is None:
            self._index_parser = ResidueIndexParser(self.launch_dir)
            self._heme_environments = self._index_parser.parse_file()
        return self._heme_environments

    def compute_reorganization_energy(self,
                                  sequence: List[int],
                                  pdb_file: str) -> Tuple[List[float], List[float]]:
        """
        Compute reorganization energy for a sequence of hemes.

        Args:
            sequence: List of heme residue IDs
            pdb_file: Path to input PDB file

        Returns:
            Tuple of (lambda_values, dielectric_constants) where:
            - lambda_values is a list of reorganization energies for each transfer step
            - dielectric_constants is a list of dielectric constants for each step
        """

        # Lazy load heme environments
        heme_environments = self._load_heme_environments()

        # Initialize for dielectric constant tracking
        dielectric_constants = []

        # Verify all hemes in sequence have environment info
        for heme_id in sequence:
            if heme_id not in heme_environments:
                raise ValueError(f"No environment information found for heme {heme_id}")

        # Check for existing calculations
        lambda_file = self.ee_dir / "Lambda.txt"
        previous_results = self._read_existing_results(lambda_file)

        if previous_results:
            es_values = previous_results.get('dielectric_constants', [])
            lambda_values = previous_results.get('lambda_values', [])

            if lambda_values and len(lambda_values) == len(sequence) - 1 and es_values:
                use_existing = self.interaction_manager.yes_no_prompt(
                    "use_existing_lambda",
                    "\nFound existing reorganization energy calculations. Use these values?"
                    )
                if use_existing:
                    print("\nUsing existing reorganization energies:")

                    for idx, value in enumerate(lambda_values):
                        print(
                            f"lambda_value_{idx}",
                            f"Step {idx+1}: λ = {value:.1f} meV",
                        )
                    return lambda_values, es_values

        # Ask calculation method
        method = self.interaction_manager.prompt(
            "lambda_method",
            "\nHow would you like to determine reorganization energies?\n"
            "1) Compute from structure using SASA analysis\n"
            "2) Enter values manually for each step\n"
            "Choice: ",
            choices=['1', '2']
        )

        if method == '2':
            lambda_values = self._get_manual_values(len(sequence) - 1)
            # For manual entry, ask for dielectric constant
            dielectric = self.interaction_manager.prompt(
                "manual_dielectric",
                "\nEnter the internal dielectric constant to use: ",
                input_type=float
            )
            dielectric_constants = [dielectric]
            self._write_to_file(sequence, lambda_values, None, dielectric_constants)
            return lambda_values, dielectric_constants

        print(f"\nComputing reorganization energies from structure {pdb_file}...")

        # Process each consecutive pair
        lambda_values = []
        sasa_results = []

        input_pdb = Path(pdb_file)
        if not input_pdb.exists():
            raise FileNotFoundError(f"Input PDB not found: {input_pdb}")

        for i in range(len(sequence) - 1):
            donor_id = sequence[i]
            acceptor_id = sequence[i + 1]

            print(f"\nAnalyzing transfer step {i+1}: Heme {donor_id} → Heme {acceptor_id}")

            # Get environments
            donor_env = heme_environments[donor_id]
            acceptor_env = heme_environments[acceptor_id]

            # Calculate SASA
            sasa_file = self.ee_dir / f"sasa_{donor_id}_{acceptor_id}.dat"
            result = None

            if sasa_file.exists():
                print(f"\nFound existing SASA results for step {i+1}")
                result = self._read_sasa_result(sasa_file)

            if result is None:
                try:
                    result = self.sasa_calculator.calculate(
                        str(input_pdb),
                        donor_env,
                        acceptor_env
                    )
                    self._write_sasa_result(sasa_file, result)
                except Exception as e:
                    print(f"\nError calculating SASA: {str(e)}",)

                    lambda_value = self.interaction_manager.prompt(
                        f"lambda_manual_{i}",
                        f"Enter reorganization energy for step {i+1} (meV): ",
                        input_type=float
                    )
                    lambda_values.append(lambda_value)
                    dielectric_constants.append(None)  # Placeholder for manual entry
                    continue

            # Calculate dielectric constant
            total_sasa = result.donor_sasa + result.acceptor_sasa
            es = self.params.alpha + (self.params.beta * total_sasa)
            dielectric_constants.append(es)

            sasa_results.append(result)
            lambda_value = self._compute_lambda(result)
            lambda_values.append(lambda_value)

            print(
                f"\nStep {i+1} Results:\n"
                f"  Donor SASA:     {result.donor_sasa:.1f} Å²\n"
                f"  Acceptor SASA:  {result.acceptor_sasa:.1f} Å²\n"
                f"  Distance:       {result.distance:.1f} Å\n"
                f"  Lambda:         {lambda_value:.3f} meV\n"
                f"  Dielectric:     {es:.2f}",
            )

        print(f"\nDielectric constants: {', '.join(f'{d:.3f}' if d is not None else 'N/A' for d in dielectric_constants)}")

        # Write results to file
        self._write_to_file(sequence, lambda_values, sasa_results, dielectric_constants)

        return lambda_values, dielectric_constants

    def _read_existing_results(self, filename: Path) -> Optional[Dict]:
        """Read existing results from Lambda.txt if it exists."""
        try:
            results = {'lambda_values': [], 'dielectric_constants': []}
            with open(filename) as f:
                for line in f:
                    if "Reorg. Eng. =" in line:
                        value = float(line.split("=")[-1].strip())
                        results['lambda_values'].append(value)
                    if "Es        =" in line:
                        value = float(line.split("=")[-1].strip())
                        results['dielectric_constants'].append(value)
            return results if results['lambda_values'] else None
        except (FileNotFoundError, ValueError):
            return None

    def _get_manual_values(self, num_steps: int) -> List[float]:
        """Get manually entered lambda values."""
        values = []
        for i in range(num_steps):
            value = self.interaction_manager.prompt(
                f"lambda_manual_{i}",
                f"Enter reorganization energy for step {i+1} (meV): ",
                input_type=float
            )
            values.append(value)
        return values

    def _compute_lambda(self, result: SASAResult) -> float:
        """Compute reorganization energy from SASA results using correct formula."""
        total_sasa = result.donor_sasa + result.acceptor_sasa
        es = self.params.alpha + (self.params.beta * total_sasa)

        # Dielectric term (M)
        M = (1 / self.params.eopt) - (1 / es)

        # Geometric term (R) with Bohr radius conversion
        bohr = 0.53  # Bohr radius in Angstroms
        R = (1 / ((2 * self.params.rd) / bohr)) + \
            (1 / ((2 * self.params.ra) / bohr)) - \
            (1 / (result.distance / bohr))

        # Final calculation with conversion to meV
        lambda_out = M * R * 27.2114 
        return lambda_out

    def _write_sasa_result(self, file_path: Path, result: SASAResult):
        """Write SASA results to file."""
        with open(file_path, 'w') as f:
            f.write(f"{result.donor_sasa} {result.acceptor_sasa} {result.distance}\n")

    def _read_sasa_result(self, file_path: Path) -> Optional[SASAResult]:
        """Read SASA results from file."""
        try:
            with open(file_path) as f:
                donor_sasa, acceptor_sasa, distance = map(float, f.readline().split())
                return SASAResult(donor_sasa, acceptor_sasa, distance)
        except (FileNotFoundError, ValueError):
            return None

    def _write_to_file(self,
                    sequence: List[int],
                    lambda_values: List[float],
                    sasa_results: Optional[List[SASAResult]] = None,
                    dielectric_constants: Optional[List[float]] = None) -> None:
        """Write reorganization energies and calculation details to file."""
        with open(self.ee_dir / "Lambda.txt", 'w') as f:
            for idx, (value, donor, acceptor) in enumerate(zip(lambda_values,
                                                        sequence[:-1],
                                                        sequence[1:])):
                f.write(f"\nHEM-{donor} -> HEM-{acceptor} ---------\n")
                if sasa_results and dielectric_constants and idx < len(sasa_results):
                    result = sasa_results[idx]
                    f.write(f"Dsasa     = {result.donor_sasa:.3f}\n")
                    f.write(f"Asasa     = {result.acceptor_sasa:.3f}\n")
                    f.write(f"Rda       = {result.distance:.3f}\n")
                    f.write(f"TotalSASA = {result.donor_sasa + result.acceptor_sasa:.3f}\n")
                    f.write(f"Es        = {dielectric_constants[idx]:.3f}\n")
                f.write(f"----------------------------\n")
                f.write(f"Reorg. Eng. = {value:.3f}\n")
