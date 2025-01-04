"""
Diffusion Calculator module for BioDC.
Handles computation of analytical charge diffusion constants.
"""
import math
from pathlib import Path
from typing import List, Tuple, Dict, Optional, Union
from dataclasses import dataclass
import re
import numpy as np
from rich.console import Console
from rich.table import Table

from biodc.utils.interaction import InteractionManager
from biodc.utils.structure_analyzer import PDBProcessor
from biodc.core.kineval_modules.derrida import VD

@dataclass
class DiffusionResult:
    """Container for diffusion calculation results."""
    velocity: float  # sites/s
    diffusion_coeff: float  # sites²/s
    diffusion_coeff_phys: float  # cm²/s
    avg_spacing: float  # Å
    chain_length: float  # Å
    sequence: List[int]
    geometry_types: List[str]
    conductivity: float  # S/cm
    current_subunit: float  # pA
    current_fixed: float  # pA
    charge_density: float  # charges/cm³

class DiffusionCalculator:
    """Calculator for analytical charge diffusion constants."""
    
    def __init__(self,
                interaction_manager: InteractionManager,
                launch_dir: Path,
                ee_dir: Path):
        """
        Initialize diffusion calculator.
        
        Args:
            interaction_manager: Manages user interactions
            launch_dir: Project launch directory
            ee_dir: Directory containing energetics results
        """
        self.interaction_manager = interaction_manager
        self.launch_dir = launch_dir
        self.ee_dir = ee_dir
        self.pdb_processor = PDBProcessor()
        
        # Set standard distance cutoff
        self.pdb_processor.distance_cutoff = 15.0
    
    def _parse_rates_file(self, filepath: Path) -> List[Tuple[str, str, float, float, str]]:
        """Parse rates file for donor-acceptor pairs and rates."""
        rates_data = []
        with open(filepath, 'r') as f:
            for line in f:
                pattern = r'(HEM-\d+) -> (HEM-\d+); kf = ([\d.E+-]+) s\^-1; kb = ([\d.E+-]+) s\^-1; geometry = ([STU])'
                match = re.match(pattern, line)
                if match:
                    donor, acceptor, kf, kb, geometry = match.groups()
                    rates_data.append((donor, acceptor, float(kf), float(kb), geometry))
        return rates_data

    def _parse_dg_file(self, filepath: Path) -> List[Tuple[str, str, float]]:
        """Parse DG.txt file to extract donor-acceptor pairs and DG values."""
        dg_data = []
        with open(filepath, 'r') as f:
            for line in f:
                pattern = r'\((HEM-\d+).+?\) -> \((HEM-\d+).+?\); DG =\s*([-\d.]+) eV'
                match = re.match(pattern, line)
                if match:
                    donor, acceptor, dg = match.groups()
                    dg_data.append((donor, acceptor, float(dg)))
        return dg_data

    def _parse_lambda_file(self, filepath: Path) -> List[Tuple[str, str, float]]:
        """Parse Lambda.txt file to extract donor-acceptor pairs and reorganization energies."""
        lambda_data = []
        with open(filepath, 'r') as f:
            lines = f.readlines()
            i = 0
            while i < len(lines):
                line = lines[i].strip()
                if '->' in line:
                    parts = line.split('->')
                    donor = 'HEM' + parts[0].strip()[3:]
                    acceptor = parts[1].strip().split()[0]
                    while i < len(lines):
                        if 'Reorg. Eng.' in lines[i]:
                            reorg = float(lines[i].split('=')[1].strip())
                            lambda_data.append((donor, acceptor, reorg))
                            break
                        i += 1
                i += 1
        return lambda_data

    def _parse_hda_file(self, filepath: Path) -> List[Tuple[str, str, float]]:
        """Parse Hda.txt file to extract donor-acceptor pairs and Hda values in eV."""
        hda_data = []
        with open(filepath, 'r') as f:
            for line in f:
                if 'Hda' in line:
                    pattern = r'Hda\((HEM-\d+) <-> (HEM-\d+)\).+?Hda = +?([\d.]+) meV'
                    match = re.match(pattern, line)
                    if match:
                        donor, acceptor, hda = match.groups()
                        hda_data.append((donor, acceptor, float(hda)/1000.0))  # Convert meV to eV
        return hda_data

    def _get_geometry(self, hda: float) -> str:
        """Determine geometry type based on Hda value."""
        if hda > 0.006:  # Close to 0.008 (8 meV)
            return 'S'
        elif hda < 0.004:  # Close to 0.002 (2 meV)
            return 'T'
        else:  # Around 0.005 (5 meV)
            return 'U'

    def _compute_marcus_rates(self, hda: float, lambda_reorg: float, 
                           dg: float, T: float = 300.0) -> Tuple[float, float]:
        """Calculate Marcus electron transfer rates."""
        PI = 3.141592654
        KB = 8.6173304E-5
        HBAR = 6.582119514E-16
        
        prefactor = (2 * PI * hda**2) / (HBAR * math.sqrt(4 * PI * lambda_reorg * KB * T))
        
        e_act_forward = ((dg + lambda_reorg)**2) / (4 * lambda_reorg)
        e_act_backward = ((-1 * dg + lambda_reorg)**2) / (4 * lambda_reorg)
        
        k_forward = prefactor * math.exp(-1 * e_act_forward / (KB * T))
        k_backward = prefactor * math.exp(-1 * e_act_backward / (KB * T))
        
        return k_forward, k_backward

    def _find_matching_value(self, pairs_data: List[Tuple[str, str, float]], 
                          donor: str, acceptor: str) -> Optional[float]:
        """Find matching value for a donor-acceptor pair in parsed data."""
        for d, a, val in pairs_data:
            if d == donor and a == acceptor:
                return val
        return None

    def _compute_conductivity(self, D_cm2_s: float, rho: float, T: float = 300.0) -> float:
        """
        Compute conductivity from diffusion coefficient using Einstein relation.

        Args:
            D_cm2_s: Diffusion coefficient in cm²/s
            rho: Charge density in charges/cm³
            T: Temperature in Kelvin

        Returns:
            Conductivity in S/cm
        """
        KB = 1.38E-23  # J/K
        E_CHARGE = 1.602E-19  # C
        return (E_CHARGE**2 * rho * D_cm2_s) / (KB * T)

    def _compute_current(self, D_cm2_s: float, rho: float, length_cm: float,
                       T: float = 300.0, V: float = 0.1) -> float:
        """
        Compute diffusive current for a given length.

        Args:
            D_cm2_s: Diffusion coefficient in cm²/s
            rho: Charge density in charges/cm³
            length_cm: Wire length in cm
            T: Temperature in Kelvin
            V: Applied voltage in V

        Returns:
            Current in Amperes
        """
        KB = 1.38E-23
        E_CHARGE = 1.602E-19
        R = 7.5E-8  # Wire radius in cm
        
        A = np.pi * R * R  # Cross-sectional area
        
        return ((A * E_CHARGE * E_CHARGE * rho * D_cm2_s) / (KB * T * length_cm)) * V

    def compute_diffusion_constant(self) -> DiffusionResult:
        """
        Compute analytical diffusion constant using rates from EE directory.
        
        Returns:
            DiffusionResult containing computed parameters
        """
        # Check for existence of min.pdb
        pdb_path = self.ee_dir / "min.pdb"
        if not pdb_path.exists():
            raise FileNotFoundError(f"Required PDB file not found: {pdb_path}")
            
        # First try to find rates file
        rates_path = self.ee_dir / "rates.txt"
        forward_rates = []
        backward_rates = []
        geometries = []
        
        if rates_path.exists():
            print("Using pre-computed rates from rates.txt")
            rates_data = self._parse_rates_file(rates_path)
            
            for _, _, kf, kb, geom in rates_data:
                forward_rates.append(kf)
                backward_rates.append(kb)
                geometries.append(geom)
                
        else:
            print("Computing rates from energetic parameters...")
            # Load required files
            dg_path = self.ee_dir / "DG.txt"
            lambda_path = self.ee_dir / "Lambda.txt"
            hda_path = self.ee_dir / "Hda.txt"
            
            if not all(p.exists() for p in [dg_path, lambda_path, hda_path]):
                raise FileNotFoundError(
                    "Missing required files. Need either rates.txt or "
                    "all of: DG.txt, Lambda.txt, Hda.txt"
                )
            
            # Parse all files first
            dg_data = self._parse_dg_file(dg_path)
            lambda_data = self._parse_lambda_file(lambda_path)
            hda_data = self._parse_hda_file(hda_path)
            
            # Compute rates for each step
            for donor, acceptor, dg in dg_data:
                # Find matching lambda and hda values
                lambda_val = self._find_matching_value(lambda_data, donor, acceptor)
                hda_val = self._find_matching_value(hda_data, donor, acceptor)
                
                if lambda_val is not None and hda_val is not None:
                    kf, kb = self._compute_marcus_rates(hda_val, lambda_val, dg)
                    forward_rates.append(kf)
                    backward_rates.append(kb)
                    geometries.append(self._get_geometry(hda_val))
        
        # Get chain information from PDB
        atoms_dict = self.pdb_processor.read_pdb_atoms(str(pdb_path))
        sequence = self.pdb_processor.detect_sequence(atoms_dict, topology='linear')
        
        # Calculate average Fe-Fe spacing
        distances = self.pdb_processor.measure_avg_dist(str(pdb_path), topology='linear')
        avg_spacing = np.mean([d['distance'] for d in distances])
        
        # Calculate full chain length
        chain_length = self.pdb_processor.measure_fe_fe_distance(
            str(pdb_path), sequence[0], sequence[-1]
        )
        
        # Compute Derrida parameters
        V, D = VD(forward_rates, backward_rates)
        if V is None or D is None:
            raise ValueError("Failed to compute Derrida parameters")
            
        # Convert D to physical units
        D_cm2_s = D * (avg_spacing * 1E-8)**2

        # Calculate wire geometry and charge density
        R = 7.5E-8  # Wire radius in cm
        A = np.pi * R * R
        chain_length_cm = chain_length * 1E-8
        V_element = A * chain_length_cm
        rho = (0.5 * len(sequence)) / V_element

        # Calculate conductivity
        conductivity = self._compute_conductivity(D_cm2_s, rho)

        # Calculate currents
        FIXED_L = 8.27E-5  # Fixed reference length in cm
        subunit_length_cm = chain_length * 1E-8  # Convert Å to cm
        
        current_subunit = self._compute_current(D_cm2_s, rho, subunit_length_cm)
        current_fixed = self._compute_current(D_cm2_s, rho, FIXED_L)
        
        return DiffusionResult(
            velocity=V,
            diffusion_coeff=D,
            diffusion_coeff_phys=D_cm2_s,
            avg_spacing=avg_spacing,
            chain_length=chain_length,
            sequence=sequence,
            geometry_types=geometries,
            conductivity=conductivity,
            current_subunit=current_subunit * 1E12,  # Convert to pA
            current_fixed=current_fixed * 1E12,  # Convert to pA
            charge_density=rho
        )

    def display_results(self, results: DiffusionResult) -> None:
        """Display calculation results in a formatted table."""
        console = Console()
        
        # Create main results table
        table = Table(title="Diffusion Constant Analysis Results")
        
        # Add main parameters
        table.add_column("Parameter", style="cyan")
        table.add_column("Value", justify="right")
        table.add_column("Units", style="green")
        
        table.add_row("Velocity", f"{results.velocity:.2E}", "sites/s")
        table.add_row("Diffusion Coefficient (lattice)", f"{results.diffusion_coeff:.2E}", "sites²/s")
        table.add_row("Diffusion Coefficient (physical)", f"{results.diffusion_coeff_phys:.2E}", "cm²/s")
        table.add_row("Average Fe-Fe Spacing", f"{results.avg_spacing:.2f}", "Å")
        table.add_row("Chain Length", f"{results.chain_length:.2f}", "Å")
        table.add_row("Conductivity", f"{results.conductivity:.2E}", "S/cm")
        table.add_row("Diffusive Current (Chain Length)", f"{results.current_subunit:.2E}", "pA")
        table.add_row("Diffusive Current (Fixed Length)", f"{results.current_fixed:.2E}", "pA")
        table.add_row("Charge Density", f"{results.charge_density:.2E}", "charges/cm³")
        
        console.print(table)
        
        # Print sequence information
        print("\nHeme Sequence:")
        print(" → ".join(map(str, results.sequence)))
        print("\nGeometry Types:")
        print(" → ".join(results.geometry_types))
        
    def compute_current_vs_length(self, D_cm2_s: float, n_hemes: int, 
                                   wire_length_ang: float,
                                   start_nm: float = 1.0,
                                   end_nm: float = 1000.0,
                                   steps: int = 1000) -> List[Tuple[float, float]]:
        """
        Compute diffusive current as a function of length.

        Args:
            D_cm2_s: Diffusion coefficient in cm²/s
            n_hemes: Number of hemes in the longest chain
            wire_length_ang: Length of the longest chain in Angstroms
            start_nm: Starting length in nanometers
            end_nm: Ending length in nanometers
            steps: Number of length points to calculate

        Returns:
            List of tuples (length_nm, current_pA)
        """
        if wire_length_ang <= 0:
            raise ValueError(f"Wire length must be positive, got {wire_length_ang} Å")

        # Constants
        KB = 1.38E-23  # J/K
        T = 300  # K
        V = 0.1  # V
        R = 7.5E-8  # cm (wire radius)
        E_CHARGE = 1.602E-19  # C

        # Calculate geometry
        A = np.pi * R * R  # cm²

        # Calculate charge density using actual wire length and chain heme count
        wire_length_cm = wire_length_ang * 1E-8
        V_element = A * wire_length_cm
        rho = (0.5 * n_hemes) / V_element  # charges/cm³

        # Generate length points in nm
        lengths_nm = np.linspace(start_nm, end_nm, steps)

        results = []
        for L_nm in lengths_nm:
            # Convert length to cm for calculations
            L_cm = L_nm * 1E-7

            # Calculate current using fixed charge density
            I_diff = ((A * E_CHARGE * E_CHARGE * rho * D_cm2_s) / (KB * T * L_cm)) * V

            # Convert to picoamps
            I_diff_pA = I_diff * 1E12

            results.append((L_nm, I_diff_pA))

        return results
        if output_dir is None:
            output_dir = self.launch_dir / "KE" / "diffusion"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Save main results
        results_file = output_dir / "diffusion_results.txt"
        with open(results_file, 'w') as f:
            f.write("Diffusion Constant Analysis Results\n")
            f.write("=" * 40 + "\n\n")
            
            # Write parameters
            f.write(f"Velocity: {results.velocity:.2E} sites/s\n")
            f.write(f"Diffusion Coefficient (lattice): {results.diffusion_coeff:.2E} sites²/s\n")
            f.write(f"Diffusion Coefficient (physical): {results.diffusion_coeff_phys:.2E} cm²/s\n")
            f.write(f"Average Fe-Fe Spacing: {results.avg_spacing:.2f} Å\n")
            f.write(f"Chain Length: {results.chain_length:.2f} Å\n")
            f.write(f"Conductivity: {results.conductivity:.2E} S/cm\n")
            f.write(f"Diffusive Current (Chain Length): {results.current_subunit:.2E} pA\n")
            f.write(f"Diffusive Current (Fixed Length): {results.current_fixed:.2E} pA\n")
            f.write(f"Charge Density: {results.charge_density:.2E} charges/cm³\n\n")
            
            # Write sequence information
            f.write("Heme Sequence:\n")
            f.write(" → ".join(map(str, results.sequence)) + "\n\n")
            f.write("Geometry Types:\n")
            f.write(" → ".join(results.geometry_types) + "\n")
            
        print(f"\nResults saved to: {results_file}")
