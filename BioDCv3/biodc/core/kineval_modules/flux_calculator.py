"""
Flux Calculator module for BioDC.
Handles computation of steady-state electron flux through multi-heme chains.
"""
from pathlib import Path
from typing import List, Tuple, Dict, Optional, Any
from dataclasses import dataclass
import re
import numpy as np

from biodc.utils.interaction import InteractionManager
from biodc.core.kineval_modules.hopping import solve_flux

@dataclass
class ElectronTransferStep:
    """Container for electron transfer step information."""
    donor: str
    acceptor: str
    dg: float           # Driving force (eV)
    lambda_reorg: float # Reorganization energy (eV)
    hda: float         # Electronic coupling (eV)
    geometry: str      # 'S' for slip-stacked, 'T' for T-shaped, 'U' for unknown

@dataclass
class RedoxState:
    """Container for redox state calculation results."""
    potentials: np.ndarray      # Redox potentials (eV)
    populations: np.ndarray     # Site populations
    forward_rates: np.ndarray   # Forward rates (s⁻¹)
    backward_rates: np.ndarray  # Backward rates (s⁻¹)
    forward_flux: float        # Net forward flux (s⁻¹)
    backward_flux: float       # Net backward flux (s⁻¹)
    convergence_iterations: int # Number of iterations to converge
    dgs: np.ndarray            # ΔG values for each step
    e_act_forward: np.ndarray  # Forward activation energies
    e_act_backward: np.ndarray # Backward activation energies

class FluxCalculator:
    """Calculator for steady-state electron flux."""
    
    def __init__(self,
                interaction_manager: InteractionManager,
                launch_dir: Path,
                ee_dir: Path):
        """
        Initialize flux calculator.
        
        Args:
            interaction_manager: Manages user interactions
            launch_dir: Project launch directory
            ee_dir: Directory containing energetics results
        """
        self.interaction_manager = interaction_manager
        self.launch_dir = launch_dir
        self.ee_dir = ee_dir
        
    def _parse_rates_file(self, filepath: Path) -> List[Tuple[str, str, float, float, str]]:
        """
        Parse rates file for donor-acceptor pairs and rates.
        Expected format: HEM-X -> HEM-Y; kf = X.XXE+XX s^-1; kb = X.XXE+XX s^-1; geometry = S/T/U
        """
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
        """Parse DG.txt file to extract donor-acceptor pairs and ΔG values."""
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
        """
        Parse Hda.txt file to extract donor-acceptor pairs and Hda values.
        Converts meV to eV.
        """
        hda_data = []
        with open(filepath, 'r') as f:
            for line in f:
                if 'Hda' in line:
                    pattern = r'Hda\((HEM-\d+) <-> (HEM-\d+)\).+?Hda = +?([\d.]+) meV'
                    match = re.match(pattern, line)
                    if match:
                        donor, acceptor, hda = match.groups()
                        hda_data.append((donor, acceptor, float(hda)/1000.0))  # meV to eV
        return hda_data
        
    def _get_geometry(self, hda: float) -> str:
        """Determine geometry type based on Hda value."""
        if hda > 0.006:    # Close to 0.008 (8 meV)
            return 'S'
        elif hda < 0.004:  # Close to 0.002 (2 meV)
            return 'T'
        else:  # Around 0.005 (5 meV)
            return 'U'

    def _read_E_matrix(self, filepath: Path) -> np.ndarray:
        """
        Read E matrix from StateEnergies file and convert from meV to eV.
        """
        matrix = []
        matrix_started = False

        with open(filepath, 'r') as f:
            for line in f:
                if not line.strip():
                    continue

                if line.strip().startswith('['):
                    matrix_started = True

                if matrix_started:
                    clean_line = line.strip('[] \n')
                    if clean_line:
                        # Convert meV to eV during reading
                        row = [float(x.strip())/1000.0 for x in clean_line.split(',')]
                        matrix.append(row)

        if not matrix:
            raise ValueError("No matrix data found in file")

        n = len(matrix)
        if not all(len(row) == n for row in matrix):
            raise ValueError("Input matrix must be square")

        return np.array(matrix)

    def _read_redox_potentials(self, filepath: Path) -> np.ndarray:
        """Extract unique redox potentials following chain structure from DG.txt"""
        heme_potentials = {}  # Dictionary to store unique heme:potential pairs

        with open(filepath, 'r') as f:
            for line in f:
                matches = re.findall(r'(HEM-\d+)\s*=\s*([-\d.]+)\s*eV', line)
                if matches:
                    for heme, potential in matches:
                        heme_potentials[heme] = float(potential)

        # Convert to ordered list following chain structure
        unique_potentials = []
        ordered_hemes = []

        with open(filepath, 'r') as f:
            # Get first donor from first line
            first_line = f.readline()
            first_match = re.search(r'\((HEM-\d+)', first_line)
            if first_match:
                first_heme = first_match.group(1)
                ordered_hemes.append(first_heme)
                unique_potentials.append(heme_potentials[first_heme])

        # Follow chain to maintain order
        for heme in heme_potentials:
            if heme not in ordered_hemes:
                ordered_hemes.append(heme)
                unique_potentials.append(heme_potentials[heme])

        return np.array(unique_potentials)

    def _compute_marcus_rates(self, hda: float, lambda_reorg: float, 
                                dg: float, T: float = 300.0) -> Tuple[float, float, float, float]:
        """
        Calculate Marcus electron transfer rates.
        
        Args:
            hda: Electronic coupling (eV)
            lambda_reorg: Reorganization energy (eV)
            dg: Driving force (eV)
            T: Temperature (K)
            
        Returns:
            Tuple of (k_forward, k_backward, E_act_forward, E_act_backward)
        """
        PI = 3.141592654
        KB = 8.6173304E-5  # eV/K
        HBAR = 6.582119514E-16  # eV·s
        
        prefactor = (2 * PI * hda**2) / (HBAR * math.sqrt(4 * PI * lambda_reorg * KB * T))
        
        # Calculate activation energies
        e_act_forward = ((dg + lambda_reorg)**2) / (4 * lambda_reorg)
        e_act_backward = ((-1 * dg + lambda_reorg)**2) / (4 * lambda_reorg)
        
        # Calculate rates
        k_forward = prefactor * math.exp(-1 * e_act_forward / (KB * T))
        k_backward = prefactor * math.exp(-1 * e_act_backward / (KB * T))
        
        return k_forward, k_backward, e_act_forward, e_act_backward

    def _modify_E_matrix(self, E: np.ndarray, 
                       diagonal_shift: float = 0.0,
                       offdiagonal_scale: float = 1.0) -> np.ndarray:
        """
        Modify energy matrix by shifting diagonal and scaling off-diagonal elements.
        
        Args:
            E: Original energy matrix (eV)
            diagonal_shift: Shift to apply to diagonal elements (eV)
            offdiagonal_scale: Scale factor for off-diagonal elements
            
        Returns:
            Modified energy matrix
        """
        E_modified = E.copy()
        n = len(E)

        # Shift diagonal elements
        np.fill_diagonal(E_modified, E_modified.diagonal() + diagonal_shift)
        
        # Scale off-diagonal elements
        off_diag_mask = ~np.eye(n, dtype=bool)
        E_modified[off_diag_mask] *= offdiagonal_scale
        
        return E_modified

    def _compute_redox_potentials(self, E: np.ndarray, case: str,
                               populations: Optional[np.ndarray] = None) -> np.ndarray:
        """
        Compute redox potentials based on case and populations.
        
        Args:
            E: Energy matrix (eV)
            case: One of ['reduced', 'oxidized', 'mixed']
            populations: Population array for mixed case
            
        Returns:
            Array of redox potentials
        """
        n_hemes = len(E)
        R = np.zeros(n_hemes)

        if case == 'reduced':
            R = E.diagonal().copy()
        elif case == 'oxidized':
            R = E.diagonal() + E.sum(axis=1) - E.diagonal()
        elif case == 'mixed':
            if populations is None:
                raise ValueError("Populations required for mixed case")
            R = E.diagonal() + np.sum(E * (1 - populations), axis=1) - E.diagonal()
        else:
            raise ValueError(f"Unknown case: {case}")
            
        return R

    def _adaptive_mixing_parameter(self, iteration: int,
                                max_diff: float,
                                prev_max_diff: float,
                                current_param: float,
                                max_diff_history: Optional[List[float]] = None,
                                window_size: int = 3) -> Tuple[float, List[float]]:
        """
        Adaptively adjust mixing parameter based on convergence behavior.
        
        Args:
            iteration: Current iteration number
            max_diff: Current maximum difference
            prev_max_diff: Previous maximum difference
            current_param: Current mixing parameter
            max_diff_history: List storing recent max_diff values
            window_size: Number of previous iterations to consider
            
        Returns:
            Tuple of (new_mixing_parameter, updated_history)
        """
        if max_diff_history is None:
            max_diff_history = []

        # Update history
        max_diff_history.append(max_diff)
        if len(max_diff_history) > window_size:
            max_diff_history.pop(0)

        # Don't adjust parameter until we have enough history
        if iteration < 2:
            return current_param, max_diff_history

        # Detect oscillations
        is_oscillating = False
        if len(max_diff_history) >= 3:
            diffs = [max_diff_history[i] - max_diff_history[i-1]
                    for i in range(1, len(max_diff_history))]
            sign_changes = sum(1 for i in range(1, len(diffs))
                            if diffs[i] * diffs[i-1] < 0)
            is_oscillating = sign_changes >= (len(diffs) - 1) / 2

        # Determine new mixing parameter
        if is_oscillating:
            new_param = max(0.1, current_param * 0.7)  # More aggressive reduction
        elif max_diff > prev_max_diff:
            new_param = max(0.1, current_param * 0.9)  # Normal reduction
        else:
            new_param = min(0.8, current_param * 1.05)  # Gentle increase

        return new_param, max_diff_history

    def _calculate_all_redox_states(self, 
                                  E: np.ndarray,
                                  H: np.ndarray,
                                  lambda_reorg: np.ndarray,
                                  adaptive_mixing: bool = True) -> Dict[str, RedoxState]:
        """
        Calculate all redox states (reduced, oxidized, mixed).

        Args:
            E: Energy matrix (eV)
            H: Electronic coupling values (eV)
            lambda_reorg: Reorganization energies (eV)
            adaptive_mixing: Whether to use adaptive mixing parameter

        Returns:
            Dictionary containing RedoxState objects for each case
        """
        states = {}
        for case in ['reduced', 'oxidized', 'mixed']:
            print(f"\n  Calculating {case} state...")
            states[case] = self._calculate_redox_state(
                E=E,
                H=H,
                lambda_reorg=lambda_reorg,
                case=case,
                adaptive_mixing=adaptive_mixing
            )
            print(f"    Converged in {states[case].convergence_iterations} iterations")
            print(f"    Forward flux: {states[case].forward_flux:.2E} s⁻¹")
            print(f"    Backward flux: {states[case].backward_flux:.2E} s⁻¹")

        return states

    def _compute_flux_from_rates(self, 
                              forward_rates: List[float], 
                              backward_rates: List[float]) -> Tuple[float, float]:
        """
        Compute forward and backward flux using hopping model.
        
        Args:
            forward_rates: List of forward rate constants
            backward_rates: List of backward rate constants
            
        Returns:
            Tuple of (forward_flux, backward_flux)
        """
        # Forward flux calculation
        forward_sol = solve_flux(forward_rates, backward_rates)
        forward_flux = forward_sol[-1]

        # Backward flux calculation (reverse the rates)
        backward_rates_rev = backward_rates[::-1]
        forward_rates_rev = forward_rates[::-1]
        backward_sol = solve_flux(backward_rates_rev, forward_rates_rev)
        backward_flux = backward_sol[-1]

        return forward_flux, backward_flux

    def compute_steady_state_flux(self) -> Dict[str, Any]:
        """
        Compute steady-state electron flux through the system.
        
        Returns:
            Dictionary containing calculation results:
            - forward_flux: Forward electron transfer flux (s⁻¹)
            - backward_flux: Backward electron transfer flux (s⁻¹)
            - net_flux: Net electron transfer flux (s⁻¹)
            - currents: Dictionary of currents in different units
            - steps: List of ElectronTransferStep objects
            - rates: List of (kf, kb, Ea_f, Ea_b) tuples
            - redox_states: Optional dictionary of RedoxState objects
        """
        # Check for required files
        self._check_required_files()
        print("\nStarting steady-state flux calculation...")
        
        # First check for rates file
        rates_file = self.ee_dir / "rates.txt"
        redox_states = None
        steps = []
        results = []
        
        if rates_file.exists():
            print("Using pre-computed rates from rates.txt")
            rates_data = self._parse_rates_file(rates_file)
            
            # Create electron transfer steps with placeholder energetics
            for donor, acceptor, kf, kb, geometry in rates_data:
                step = ElectronTransferStep(
                    donor=donor,
                    acceptor=acceptor,
                    dg=0.0,           # placeholder
                    lambda_reorg=0.0,  # placeholder
                    hda=0.0,          # placeholder
                    geometry=geometry
                )
                steps.append(step)
                results.append((kf, kb, 0.0, 0.0))  # forward rate, backward rate, placeholder activation energies
                
        else:
            print("Computing rates from energetic parameters...")
            state_energies_file = self.ee_dir / "StateEnergies.txt"
            dg_file = self.ee_dir / "DG.txt"
            lambda_file = self.ee_dir / "Lambda.txt"
            hda_file = self.ee_dir / "Hda.txt"
            
            if state_energies_file.exists():
                # PATH 2: Using energy matrix
                print("Using StateEnergies.txt pathway")
                E = self._read_E_matrix(state_energies_file)
                lambda_data = self._parse_lambda_file(lambda_file)
                hda_data = self._parse_hda_file(hda_file)
                
                # Convert to arrays
                lambda_values = np.array([l for _, _, l in lambda_data])
                hda_values = np.array([h for _, _, h in hda_data])
                
                # Calculate redox states
                redox_states = self._calculate_all_redox_states(E, hda_values, lambda_values)
                mixed_state = redox_states['mixed']
                
                # Generate steps and results
                n_sites = len(E)
                for i in range(n_sites - 1):
                    donor = f"HEM-{i+1}"
                    acceptor = f"HEM-{i+2}"
                    geometry = self._get_geometry(hda_values[i])
                    dg = mixed_state.dgs[i]
                    
                    step = ElectronTransferStep(
                        donor=donor,
                        acceptor=acceptor,
                        dg=dg,
                        lambda_reorg=lambda_values[i],
                        hda=hda_values[i],
                        geometry=geometry
                    )
                    
                    steps.append(step)
                    results.append((
                        mixed_state.forward_rates[i],
                        mixed_state.backward_rates[i],
                        mixed_state.e_act_forward[i],
                        mixed_state.e_act_backward[i]
                    ))
                
            elif dg_file.exists():
                # PATH 1: Using DG.txt
                print("Using DG.txt pathway")
                dg_data = self._parse_dg_file(dg_file)
                lambda_data = self._parse_lambda_file(lambda_file)
                hda_data = self._parse_hda_file(hda_file)
                
                # Create electron transfer steps
                for donor, acceptor, dg in dg_data:
                    # Find matching lambda and hda values
                    matching_lambda = next(
                        (l for d, a, l in lambda_data if d == donor and a == acceptor),
                        None
                    )
                    matching_hda = next(
                        (h for d, a, h in hda_data if d == donor and a == acceptor),
                        None
                    )
                    
                    if matching_lambda is not None and matching_hda is not None:
                        geometry = self._get_geometry(matching_hda)
                        step = ElectronTransferStep(
                            donor=donor,
                            acceptor=acceptor,
                            dg=dg,
                            lambda_reorg=matching_lambda,
                            hda=matching_hda,
                            geometry=geometry
                        )
                        steps.append(step)
                        results.append(self._compute_marcus_rates(
                            step.hda, step.lambda_reorg, step.dg))
            
            else:
                raise FileNotFoundError(
                    "Neither StateEnergies.txt nor DG.txt found. "
                    "Please provide rate parameters."
                )
        
        # Calculate fluxes
        forward_rates = [kf for kf, _, _, _ in results]
        backward_rates = [kb for _, kb, _, _ in results]
        
        forward_flux, backward_flux = self._compute_flux_from_rates(
            forward_rates, backward_rates)
        
        # Calculate currents
        E_CHARGE = 1.602E-19  # Elementary charge in Coulombs
        net_flux = forward_flux - backward_flux
        
        currents = {
            'forward_pA': forward_flux * E_CHARGE * 1E12,
            'backward_pA': backward_flux * E_CHARGE * 1E12,
            'net_pA': net_flux * E_CHARGE * 1E12
        }
        
        # Print summary
        print("\nFlux Analysis Results:")
        print("-" * 40)
        print(f"Forward flux: {forward_flux:.2E} s⁻¹")
        print(f"Backward flux: {backward_flux:.2E} s⁻¹")
        print(f"Net flux: {net_flux:.2E} s⁻¹")
        print(f"Net current: {currents['net_pA']:.2E} pA")
        
        return {
            'forward_flux': forward_flux,
            'backward_flux': backward_flux,
            'net_flux': net_flux,
            'currents': currents,
            'steps': steps,
            'rates': results,
            'redox_states': redox_states
        }
    
    def _save_rate_analysis(self, output_dir: Path,
                              steps: List[ElectronTransferStep],
                              results: List[tuple],
                              redox_states: Optional[Dict[str, RedoxState]] = None) -> Path:
        """
        Save rate calculations for all cases.
        
        Args:
            output_dir: Directory to save results
            steps: List of electron transfer steps
            results: List of (kf, kb, Ea_f, Ea_b) tuples
            redox_states: Optional dictionary of RedoxState objects
            
        Returns:
            Path to saved rate analysis file
        """
        rate_file = output_dir / "rate_analysis.txt"
        
        with open(rate_file, "w") as f:
            f.write("Rate Calculations from Different Energy Schemes\n")
            f.write("=" * 80 + "\n\n")

            if redox_states is not None:
                # Write results for each redox state
                for state in ['reduced', 'oxidized', 'mixed']:
                    redox_state = redox_states[state]
                    f.write(f"{state.capitalize()} State:\n")
                    f.write("-" * (len(state) + 7) + "\n")
                    f.write("Step,Geometry,DG(eV),Lambda(eV),Hda(meV),E_act_forward(eV),"
                           "E_act_backward(eV),k_forward(s⁻¹),k_backward(s⁻¹)\n")
                    
                    n_steps = len(redox_state.dgs)
                    for i in range(n_steps):
                        step = steps[i]
                        f.write(f"HEM-{i+1}->HEM-{i+2},{step.geometry},"
                               f"{redox_state.dgs[i]:.3f},{step.lambda_reorg:.3f},"
                               f"{step.hda*1000:.3f},{redox_state.e_act_forward[i]:.3f},"
                               f"{redox_state.e_act_backward[i]:.3f},"
                               f"{redox_state.forward_rates[i]:.2E},"
                               f"{redox_state.backward_rates[i]:.2E}\n")
                    f.write("\n")
            else:
                # Write direct rate results
                f.write("Direct Rate Results:\n")
                f.write("-" * 40 + "\n")
                if any(step.dg != 0.0 for step in steps):  # DG.txt pathway
                    f.write("Step,Geometry,DG(eV),Lambda(eV),Hda(meV),E_act_forward(eV),"
                            "E_act_backward(eV),k_forward(s⁻¹),k_backward(s⁻¹)\n")

                    for step, (kf, kb, eaf, eab) in zip(steps, results):
                        f.write(f"{step.donor}->{step.acceptor},{step.geometry},"
                               f"{step.dg:.3f},{step.lambda_reorg:.3f},"
                               f"{step.hda*1000:.3f},{eaf:.3f},{eab:.3f},"
                               f"{kf:.2E},{kb:.2E}\n")
                else:  # rates.txt pathway
                    f.write("Step,Geometry,k_forward(s⁻¹),k_backward(s⁻¹)\n")
                    for step, (kf, kb, _, _) in zip(steps, results):
                        f.write(f"{step.donor}->{step.acceptor},{step.geometry},"
                               f"{kf:.2E},{kb:.2E}\n")
                
        return rate_file

    def _save_flux_analysis(self, output_dir: Path,
                          forward_flux: float,
                          backward_flux: float,
                          steps: List[ElectronTransferStep],
                          results: List[tuple],
                          redox_states: Optional[Dict[str, RedoxState]] = None) -> Path:
        """
        Save flux analysis results.
        
        Args:
            output_dir: Directory to save results
            forward_flux: Forward electron transfer flux
            backward_flux: Backward electron transfer flux
            steps: List of electron transfer steps
            results: List of (kf, kb, Ea_f, Ea_b) tuples
            redox_states: Optional dictionary of RedoxState objects
            
        Returns:
            Path to saved flux analysis file
        """
        flux_file = output_dir / "flux_analysis.txt"
        E_CHARGE = 1.602E-19  # Elementary charge in Coulombs
        
        with open(flux_file, "w") as f:
            f.write("Electron Transfer Flux Analysis\n")
            f.write("=" * 80 + "\n\n")
            
            # Write flux results
            net_flux = forward_flux - backward_flux
            forward_current = forward_flux * E_CHARGE * 1E12
            backward_current = backward_flux * E_CHARGE * 1E12
            net_current = net_flux * E_CHARGE * 1E12
            
            f.write("Net Flux Results:\n")
            f.write("-" * 20 + "\n")
            f.write(f"Forward Flux: {forward_flux:.2E} s⁻¹\n")
            f.write(f"Forward Current: {forward_current:.2E} pA\n")
            f.write(f"Backward Flux: {backward_flux:.2E} s⁻¹\n")
            f.write(f"Backward Current: {backward_current:.2E} pA\n")
            f.write(f"Net Flux: {net_flux:.2E} s⁻¹\n")
            f.write(f"Net Current: {net_current:.2E} pA\n\n")
            
            # Write chain sequence
            f.write("Electron Transfer Chain:\n")
            f.write("-" * 20 + "\n")
            for i, step in enumerate(steps, 1):
                f.write(f"Step {i}: {step.donor} → {step.acceptor} ({step.geometry}-stack)\n")
            f.write("\n")
            
            # Write redox state information if available
            if redox_states is not None:
                f.write("Redox State Analysis:\n")
                f.write("-" * 20 + "\n")
                for state, data in redox_states.items():
                    f.write(f"\n{state.capitalize()} State:\n")
                    net = data.forward_flux - data.backward_flux
                    forward_current = data.forward_flux * E_CHARGE * 1E12
                    backward_current = data.backward_flux * E_CHARGE * 1E12
                    net_current = net * E_CHARGE * 1E12
                    
                    f.write(f"  Forward Flux: {data.forward_flux:.2E} s⁻¹\n")
                    f.write(f"  Forward Current: {forward_current:.2E} pA\n")
                    f.write(f"  Backward Flux: {data.backward_flux:.2E} s⁻¹\n")
                    f.write(f"  Backward Current: {backward_current:.2E} pA\n")
                    f.write(f"  Net Flux: {net:.2E} s⁻¹\n")
                    f.write(f"  Net Current: {net_current:.2E} pA\n")
            
        return flux_file

    def _save_dg_flux(self, output_dir: Path,
                    redox_states: Dict[str, RedoxState]) -> Path:
        """
        Save DG and flux information from mixed state.
        
        Args:
            output_dir: Directory to save results
            redox_states: Dictionary of RedoxState objects
            
        Returns:
            Path to saved DG flux file
        """
        dg_file = output_dir / "DG_flux.txt"
        mixed_state = redox_states['mixed']
        potentials = mixed_state.potentials
        dgs = mixed_state.dgs
        
        with open(dg_file, "w") as f:
            for i in range(len(dgs)):
                donor_pot = potentials[i]
                acceptor_pot = potentials[i+1]
                dg = dgs[i]
                f.write(f"(HEM-{i+1} = {donor_pot:.3f} eV) -> "
                       f"(HEM-{i+2} = {acceptor_pot:.3f} eV); "
                       f"DG = {dg:.3f} eV\n")
                
        return dg_file

    def save_results(self, results: Dict[str, Any], output_dir: Optional[Path] = None) -> None:
        """
        Save all flux calculation results.
        
        Args:
            results: Dictionary containing calculation results
            output_dir: Optional directory to save results (defaults to KE/flux directory)
        """
        if output_dir is None:
            output_dir = self.launch_dir / "KE" / "flux"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        print("\nSaving analysis results...")
        
        # Save rate analysis
        rate_file = self._save_rate_analysis(
            output_dir=output_dir,
            steps=results['steps'],
            results=results['rates'],
            redox_states=results.get('redox_states')
        )
        print(f"Rate analysis saved to: {rate_file}")
        
        # Save flux analysis
        flux_file = self._save_flux_analysis(
            output_dir=output_dir,
            forward_flux=results['forward_flux'],
            backward_flux=results['backward_flux'],
            steps=results['steps'],
            results=results['rates'],
            redox_states=results.get('redox_states')
        )
        print(f"Flux analysis saved to: {flux_file}")
        
        # Save DG_flux.txt if redox states are available
        if results.get('redox_states'):
            dg_file = self._save_dg_flux(
                output_dir=output_dir,
                redox_states=results['redox_states']
            )
            print(f"DG flux data saved to: {dg_file}")

    def _create_gradient_span(self, ax, y_min: float, y_max: float, y_mean: float,
                             x_min: float, x_max: float, n_steps: int = 1000) -> None:
        """
        Create spans with appropriate gradients and borders for rate comparison plots.
        
        Args:
            ax: Matplotlib axis
            y_min, y_max: Y-axis limits for span
            y_mean: Mean value for gradient peak
            x_min, x_max: X-axis limits for span
            n_steps: Number of steps for gradient
        """
        if y_min == 1e2 and y_max == 1e4:  # Enzymatic turnover range
            # Create array of y positions
            y_positions = np.concatenate([
                np.logspace(np.log10(y_min), np.log10(y_mean), n_steps//2),
                np.logspace(np.log10(y_mean), np.log10(y_max), n_steps//2)
            ])
            
            # Calculate alphas - highest at mean, lowest at extremes
            alphas = np.concatenate([
                np.linspace(0.2, 0.5, n_steps//2),  # Lower half
                np.linspace(0.5, 0.2, n_steps//2)   # Upper half
            ])
            
            # Create gradient spans
            for i in range(len(y_positions)-1):
                ax.axhspan(y_positions[i], y_positions[i+1],
                          xmin=x_min, xmax=x_max,
                          color='#F5DEB3', alpha=alphas[i], zorder=0)
        else:  # Experimental ranges
            # Solid gray span
            ax.axhspan(y_min, y_max, xmin=x_min, xmax=x_max,
                      color='gray', alpha=0.2, zorder=0)
            # Add black edges at extremes
            ax.axhline(y=y_min, xmin=x_min, xmax=x_max, 
                      color='black', linewidth=1, alpha=0.5, zorder=1)
            ax.axhline(y=y_max, xmin=x_min, xmax=x_max, 
                      color='black', linewidth=1, alpha=0.5, zorder=1)

    def create_rate_distribution_plot(self, 
                                   steps: List[ElectronTransferStep],
                                   results: List[tuple],
                                   output_dir: Optional[Path] = None) -> Path:
        """
        Create box plot of rate distributions with experimental comparisons.
        
        Args:
            steps: List of electron transfer steps
            results: List of (kf, kb, Ea_f, Ea_b) tuples
            output_dir: Optional output directory
            
        Returns:
            Path to saved plot
        """
        import matplotlib.pyplot as plt
        import matplotlib.patheffects as path_effects
        from matplotlib.patches import Patch
        
        if output_dir is None:
            output_dir = self.launch_dir / "KE" / "flux"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Organize data by geometry
        s_data = []  # Slip-stacked rates
        t_data = []  # T-shaped rates
        
        for step, (kf, kb, _, _) in zip(steps, results):
            if step.geometry == 'S':
                s_data.extend([kf, kb])
            elif step.geometry == 'T':
                t_data.extend([kf, kb])
        
        # Define experimental rates
        exp_rates_S = np.array([219E6, 14E6, 105E6, 1560E6])
        exp_rates_T = np.array([125E6, 8.7E6, 114E6, 87E6])
        exp_mean_S = np.mean(exp_rates_S)
        exp_mean_T = np.mean(exp_rates_T)
        
        # Create plot
        plt.style.use('default')
        fig, ax = plt.subplots(figsize=(6.5, 6.5))
        ax.set_yscale('log')
        ax.set_ylim(1e1, 1.1e13)
        
        # Layout parameters
        total_width = 10
        middle = total_width / 2
        
        # Add experimental ranges
        self._create_gradient_span(ax, np.min(exp_rates_S), np.max(exp_rates_S), 
                                exp_mean_S, 0, middle/total_width)
        self._create_gradient_span(ax, np.min(exp_rates_T), np.max(exp_rates_T), 
                                exp_mean_T, middle/total_width, 1)
        
        # Add enzymatic turnover range
        self._create_gradient_span(ax, 1e2, 1e4, np.sqrt(1e2 * 1e4), 0, 1)
        
        # Create box plots
        positions = [total_width * 0.25, total_width * 0.75]
        data = [s_data, t_data]
        bp = ax.boxplot(data, patch_artist=True, positions=positions,
                       widths=total_width * 0.2)
        
        # Set colors
        colors = ['#4169E1', '#2E8B57']  # Royal blue for S, Sea green for T
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.6)
        
        # Add dividing line and labels
        ax.axvline(x=middle, color='black', linestyle='-', linewidth=1.5)
        
        plt.text(middle/2, ax.get_ylim()[1]*1.1, 'Slip-stacked',
                horizontalalignment='center', fontsize=14, weight='bold')
        plt.text(middle + middle/2, ax.get_ylim()[1]*1.1, 'T-shaped',
                horizontalalignment='center', fontsize=14, weight='bold')
        
        # Add experimental reference labels
        plt.text(total_width*0.25, 1e8, 'Exp. MtrC', 
                color='black', fontweight='bold', fontsize=12,
                horizontalalignment='center')
        
        plt.text(total_width*0.75, 1e7, 'Exp. STC', 
                color='black', fontweight='bold', fontsize=12,
                horizontalalignment='center')
        
        # Add enzymatic turnover text
        plt.text(middle, np.sqrt(1e2 * 1e4), 'Typical Enzymatic Turnover', 
                color='black', fontweight='bold', fontsize=12,
                horizontalalignment='center')
        
        # Configure axes
        ax.set_ylabel('Rate (s⁻¹)')
        ax.set_xticks([])  # Remove x-axis ticks
        ax.tick_params(axis='y', direction='in')
        ax.set_xlim(0, total_width)
        
        # Save plot
        plot_file = output_dir / "rate_distribution.png"
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        return plot_file

    def create_redox_plot(self,
                        redox_states: Dict[str, RedoxState],
                        output_dir: Optional[Path] = None,
                        label_step: int = 1) -> Path:
        """
        Create redox potential and population plot.
        
        Args:
            redox_states: Dictionary of RedoxState objects
            output_dir: Optional output directory
            label_step: Step size for x-axis labels
            
        Returns:
            Path to saved plot
        """
        import matplotlib.pyplot as plt
        
        if output_dir is None:
            output_dir = self.launch_dir / "KE" / "flux"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Create figure with two subplots
        plt.style.use('default')
        fig = plt.figure(figsize=(3.3, 3.3))
        gs = plt.GridSpec(2, 1, height_ratios=[1, 1], hspace=0)
        ax1 = plt.subplot(gs[0])  # Potentials
        ax2 = plt.subplot(gs[1])  # Populations
        
        # Get x values
        x = np.arange(1, len(redox_states['oxidized'].potentials) + 1)
        
        # Plot potentials
        ax1.plot(x, redox_states['oxidized'].potentials, ':', marker='s', mfc='none',
                color='blue', label='Oxidized')
        ax1.plot(x, redox_states['reduced'].potentials, ':', marker='s', mfc='none',
                color='red', label='Reduced')
        ax1.plot(x, redox_states['mixed'].potentials, ':', marker='s', mfc='none',
                color='green', label='Mixed')
        
        # Set potential plot limits
        pot_all = np.concatenate([state.potentials for state in redox_states.values()])
        pot_min, pot_max = np.min(pot_all), np.max(pot_all)
        pot_range = pot_max - pot_min
        ax1.set_ylim(pot_min - 0.1*pot_range, pot_max + 0.1*pot_range)
        
        # Plot populations
        ax2.plot(x, redox_states['mixed'].populations, ':k')
        for i, pop in enumerate(redox_states['mixed'].populations):
            gray_value = 1.0 - pop  # Darker for higher population
            ax2.plot(x[i], pop, 'o',
                    markerfacecolor=f'{gray_value:.3f}',
                    markeredgecolor='black',
                    markersize=6)
        
        # Configure axis limits and labels
        x_min = x[0] - 0.2
        x_max = x[-1] + 0.2
        for ax in [ax1, ax2]:
            ax.set_xlim(x_min, x_max)
            ax.grid(True, alpha=0.2)
            ax.tick_params(direction='in', which='both', top=True)
        
        # Configure specific subplot settings
        ax1.set_ylabel('Potential (eV)')
        plt.setp(ax1.get_xticklabels(), visible=False)
        
        ax2.set_xlabel('Heme Index')
        ax2.set_ylabel('Population')
        ax2.set_ylim(-0.05, 1.05)
        ax2.set_xticks(x[::label_step])
        
        if label_step > 1:
            current_ticks = list(ax2.get_xticks())
            if 1 not in current_ticks:
                current_ticks = [1] + current_ticks
                ax2.set_xticks(current_ticks)
        
        # Save plot
        plot_file = output_dir / "redox_potentials.png"
        plt.savefig(plot_file, dpi=150, bbox_inches='tight', pad_inches=0.1)
        plt.close()
        
        return plot_file

    def create_visualizations(self, 
                           results: Dict[str, Any],
                           output_dir: Optional[Path] = None) -> Dict[str, Path]:
        """
        Create all available visualizations for the results.
        
        Args:
            results: Dictionary of calculation results
            output_dir: Optional output directory
            
        Returns:
            Dictionary mapping plot types to file paths
        """
        if output_dir is None:
            output_dir = self.launch_dir / "KE" / "flux" / "plots"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        plot_files = {}
        
        # Create rate distribution plot
        plot_files['rates'] = self.create_rate_distribution_plot(
            steps=results['steps'],
            results=results['rates'],
            output_dir=output_dir
        )
        print(f"Rate distribution plot saved to: {plot_files['rates']}")
        
        # Create redox plot if data available
        if results.get('redox_states'):
            plot_files['redox'] = self.create_redox_plot(
                redox_states=results['redox_states'],
                output_dir=output_dir
            )
            print(f"Redox potentials plot saved to: {plot_files['redox']}")
        
        return plot_files

    def _check_required_files(self) -> None:
        """
        Verify that all required input files exist.
        Raises FileNotFoundError if any required file is missing.
        """
        required_files = {
            'min.pdb': 'minimized structure',
            'DG.txt': 'driving forces',
            'Lambda.txt': 'reorganization energies',
            'Hda.txt': 'electronic couplings'
        }
        
        # First check for rates file as alternative
        rates_file = self.ee_dir / "rates.txt"
        if rates_file.exists():
            return
            
        # If no rates file, check for individual energy files
        for filename, description in required_files.items():
            file_path = self.ee_dir / filename
            if not file_path.exists():
                raise FileNotFoundError(
                    f"Required file '{filename}' containing {description} "
                    f"not found in {self.ee_dir}. Either provide individual "
                    f"energy files or rates.txt."
                )

    def _replicate_array(self, arr: np.ndarray, n_replications: int) -> np.ndarray:
        """
        Replicate a 1D array n times by repeating all values.
        
        Args:
            arr: Original array
            n_replications: Number of times to replicate
            
        Returns:
            Replicated array
        """
        result = arr.copy()
        for _ in range(n_replications):
            result = np.concatenate([result, arr])
        return result

    def _replicate_matrix(self, matrix: np.ndarray, n_replications: int) -> np.ndarray:
        """
        Replicate a square matrix n times by sliding the original matrix.
        This maintains the chain-like connectivity in the replicated system.
        
        Args:
            matrix: Original square matrix
            n_replications: Number of times to replicate
            
        Returns:
            Replicated matrix
        """
        orig_size = matrix.shape[0]

        for n in range(n_replications):
            current_size = matrix.shape[0]
            new_size = current_size + (orig_size - 1)
            new_matrix = np.zeros((new_size, new_size))

            # Copy existing matrix to top-left corner
            new_matrix[:current_size, :current_size] = matrix

            # Fill in the slid matrix maintaining connectivity
            for i in range(orig_size):
                for j in range(orig_size):
                    if i == 0 and j == 0:
                        continue

                    new_i = i + (current_size - 1)
                    new_j = j + (current_size - 1)

                    if new_i < new_size and new_j < new_size:
                        new_matrix[new_i, new_j] = matrix[i, j]

            matrix = new_matrix

        return matrix

    def _modify_system(self, 
                    E: np.ndarray,
                    lambda_values: np.ndarray,
                    hda_values: np.ndarray,
                    diagonal_shift: float = 0.0,
                    offdiagonal_scale: float = 1.0,
                    n_replications: int = 0) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Modify and optionally replicate the system matrices and arrays.
        
        Args:
            E: Energy matrix (eV)
            lambda_values: Reorganization energies (eV)
            hda_values: Electronic coupling values (eV)
            diagonal_shift: Shift to apply to diagonal elements (eV)
            offdiagonal_scale: Scale factor for off-diagonal elements
            n_replications: Number of times to replicate the system
            
        Returns:
            Tuple of (modified_E, modified_lambda, modified_hda)
        """
        # First modify the energy matrix
        E_modified = self._modify_E_matrix(E, diagonal_shift, offdiagonal_scale)
        
        # Handle replication if requested
        if n_replications > 0:
            E_modified = self._replicate_matrix(E_modified, n_replications)
            lambda_values = self._replicate_array(lambda_values, n_replications)
            hda_values = self._replicate_array(hda_values, n_replications)
            
        return E_modified, lambda_values, hda_values

    def _calculate_redox_state(self, 
                            E: np.ndarray,
                            H: np.ndarray,
                            lambda_reorg: np.ndarray,
                            case: str,
                            initial_mixing_param: float = 0.3,
                            max_iterations: int = 1000,
                            convergence_threshold: float = 0.001,
                            adaptive_mixing: bool = True) -> RedoxState:
        """
        Calculate redox state properties using self-consistent iteration.
        
        Args:
            E: Energy matrix (eV)
            H: Electronic coupling values (eV)
            lambda_reorg: Reorganization energies (eV)
            case: One of ['reduced', 'oxidized', 'mixed']
            initial_mixing_param: Initial mixing parameter for convergence
            max_iterations: Maximum number of iterations
            convergence_threshold: Convergence criterion (eV)
            adaptive_mixing: Whether to use adaptive mixing parameter
            
        Returns:
            RedoxState object containing results
        """
        n_hemes = len(E)
        populations = np.zeros(n_hemes)
        mixing_param = initial_mixing_param
        prev_max_diff = float('inf')
        max_diff_history = []

        # Initialize redox potentials
        R = self._compute_redox_potentials(E, case, populations)

        print(f"  Starting potential updates for {case} state...")

        # Iteration loop for convergence
        for iteration in range(max_iterations):
            prev_R = R.copy()

            # Calculate intermediate rates
            kfor = np.zeros(n_hemes-1)
            kback = np.zeros(n_hemes-1)

            for i in range(n_hemes-1):
                deltaA = (R[i] - R[i+1])
                lambda_i = lambda_reorg[i]
                h_i = H[i]

                kfor[i], kback[i], _, _ = self._compute_marcus_rates(h_i, lambda_i, deltaA)

            # Calculate flux and populations
            sol = solve_flux(kfor.tolist(), kback.tolist(), verbose=False)
            populations = np.array(sol[:-1])

            # Calculate new potentials
            R_new = self._compute_redox_potentials(E, case, populations)
            max_diff = np.max(np.abs(R_new - prev_R))

            if iteration > 0:
                if adaptive_mixing:
                    mixing_param, max_diff_history = self._adaptive_mixing_parameter(
                        iteration, max_diff, prev_max_diff, mixing_param, max_diff_history)
                R = mixing_param * R_new + (1 - mixing_param) * prev_R
            else:
                R = R_new

            print(f"  Update {iteration + 1}: maximum change in potentials = {max_diff:.6f} eV")
            
            if max_diff < convergence_threshold:
                print(f"  Converged after {iteration + 1} potential updates")
                break

            prev_max_diff = max_diff

        # Calculate final values for output
        n_steps = n_hemes - 1
        final_kfor = np.zeros(n_steps)
        final_kback = np.zeros(n_steps)
        dgs = np.zeros(n_steps)
        e_act_forward = np.zeros(n_steps)
        e_act_backward = np.zeros(n_steps)

        for i in range(n_steps):
            dgs[i] = R[i] - R[i+1]
            lambda_i = lambda_reorg[i]
            h_i = H[i]

            kf, kb, eaf, eab = self._compute_marcus_rates(h_i, lambda_i, dgs[i])
            final_kfor[i] = kf
            final_kback[i] = kb
            e_act_forward[i] = eaf
            e_act_backward[i] = eab

        # Calculate final fluxes
        forward_sol = solve_flux(final_kfor.tolist(), final_kback.tolist(), verbose=False)
        backward_sol = solve_flux(final_kback[::-1].tolist(), final_kfor[::-1].tolist(), verbose=False)

        return RedoxState(
            potentials=R,
            populations=populations,
            forward_rates=final_kfor,
            backward_rates=final_kback,
            forward_flux=forward_sol[-1],
            backward_flux=backward_sol[-1],
            convergence_iterations=iteration + 1,
            dgs=dgs,
            e_act_forward=e_act_forward,
            e_act_backward=e_act_backward
        )












