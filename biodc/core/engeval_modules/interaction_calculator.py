"""
Module for calculating heme-heme interaction energies using PBSA method.
Supports all heme types and computes pairwise interaction energies.
"""

import os
import re
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Set
from dataclasses import dataclass
import logging
import itertools
from concurrent.futures import ThreadPoolExecutor
import numpy as np

try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None

from biodc.utils.interaction import InteractionManager
from biodc.utils.state_selector import RedoxStateManager, RedoxState

logger = logging.getLogger(__name__)

# ANSI color codes for colorful console printing
COLORS = {
    'HEADER': '\033[95m',
    'BLUE': '\033[94m',
    'GREEN': '\033[92m',
    'YELLOW': '\033[93m',
    'RED': '\033[91m',
    'ENDC': '\033[0m',
    'BOLD': '\033[1m',
}

@dataclass
class PBSAParameters:
    """Parameters for PBSA calculations."""
    epsin: float          # Internal dielectric constant
    epsout: float        # External dielectric constant
    istrng: float        # Ionic strength (mM)
    membraneopt: int     # Membrane option (0/1)
    epsmem: float        # Membrane dielectric constant
    mthick: float        # Membrane thickness
    poretype: int        # Pore detection (0/1)
    
    # Calculation method specific parameters
    ipb: int = 2         # PB method option
    inp: int = 2         # Non-polar solvation option
    ivalence: int = 0    # Ion valence option
    bcopt: int = 5       # Boundary condition option
    eneopt: int = 2      # Energy calculation option
    maxitn: int = 100    # Maximum iterations
    smoothopt: int = 1   # Smoothing option
    nfocus: int = 2      # Number of focusing steps

class HemeInteractionCalculator:
    """Calculator for heme-heme interaction energies using PBSA method."""

    def __init__(self, 
                 interaction_manager: InteractionManager,
                 launch_dir: Path,
                 forcefield_dir: Path,
                 pdb_file: str):
        """
        Initialize calculator.
        
        Args:
            interaction_manager: Manager for user interactions
            launch_dir: Project launch directory
            forcefield_dir: Directory containing forcefield files
            pdb_file: Path to input PDB file
        """
        self.interaction_manager = interaction_manager
        self.launch_dir = Path(launch_dir)
        self.forcefield_dir = Path(forcefield_dir)
        self.pdb_file = str(Path(pdb_file))
        self.ee_dir = self.launch_dir / "EE"
        self.ee_dir.mkdir(exist_ok=True)

    def read_and_display_energy_matrix(self) -> Optional[Dict[Tuple[int, int], float]]:
        """
        Read the energy matrix from EnergyMatrix.txt and display it.

        Returns:
            Dictionary of energies if file exists and is valid, None otherwise
        """
        matrix_file = self.ee_dir / "EnergyMatrix.txt"
        if not matrix_file.exists():
            return None

        try:
            with open(matrix_file, 'r') as f:
                # Read first line containing heme IDs
                heme_ids = [int(x) for x in f.readline().strip().split()]

                # Read matrix values
                matrix = []
                for line in f:
                    row = [float(x) for x in line.strip().split()]
                    matrix.append(row)

                # Convert to dictionary format
                energies = {}
                for i, heme1_id in enumerate(heme_ids):
                    for j, heme2_id in enumerate(heme_ids):
                        if heme1_id <= heme2_id:
                            energies[(heme1_id, heme2_id)] = matrix[i][j]

                # Display the matrix
                print("\nLoaded existing energy matrix from file:")
                self.print_energy_matrix(heme_ids, energies)

                return energies

        except Exception as e:
            print(f"Error reading energy matrix file: {e}")
            return None

    def _is_nearest_neighbor(self, i: int, j: int, sequence: List[int]) -> bool:
        """Check if two hemes are nearest neighbors in the sequence."""
        idx_i = sequence.index(i)
        idx_j = sequence.index(j)
        return abs(idx_i - idx_j) == 1

    def _format_energy_value(self, value: float, heme1_id: int, heme2_id: int, sequence: List[int]) -> str:
        """
        Format energy value with color coding based on heme relationships:
        - RED: Diagonal elements (site energies)
        - YELLOW: First nearest neighbors
        - GREEN: More distant neighbors
        """
        if heme1_id == heme2_id:  # Diagonal elements (site energies)
            return f"{COLORS['RED']}{value:8.3f}{COLORS['ENDC']}"
        elif self._is_nearest_neighbor(heme1_id, heme2_id, sequence):  # Nearest neighbors
            return f"{COLORS['YELLOW']}{value:8.3f}{COLORS['ENDC']}"
        else:  # More distant neighbors
            return f"{COLORS['GREEN']}{value:8.3f}{COLORS['ENDC']}"

    def print_energy_matrix(self, sequence: List[int], energies: Dict[Tuple[int, int], float]):
        """
        Print a nicely formatted energy matrix to the console.

        Args:
            sequence: List of heme IDs in order
            energies: Dictionary mapping heme pairs to their energies
        """
        n = len(sequence)
        matrix = [[0.0 for _ in range(n)] for _ in range(n)]

        # Fill matrix with values
        for i, heme1_id in enumerate(sequence):
            for j, heme2_id in enumerate(sequence):
                if heme1_id <= heme2_id:
                    value = energies[(heme1_id, heme2_id)]
                else:
                    value = energies[(heme2_id, heme1_id)]
                matrix[i][j] = value

        # Convert to numpy array for easier processing
        matrix_np = np.array(matrix)
        max_abs_value = np.max(np.abs(matrix_np))

        # Print header
        print(f"\n{COLORS['BOLD']}Energy Matrix (eV){COLORS['ENDC']}")
        print("Color coding:")
        print(f"{COLORS['RED']}■{COLORS['ENDC']} Site energies (diagonal elements)")
        print(f"{COLORS['YELLOW']}■{COLORS['ENDC']} First nearest neighbor interactions")
        print(f"{COLORS['GREEN']}■{COLORS['ENDC']} More distant neighbor interactions")
        print()

        # Print column headers
        print("      " + "".join(f"{heme_id:8d}" for heme_id in sequence))
        print("      " + "--------" * n)

        # Print matrix with row labels and color coding
        for i, heme_id in enumerate(sequence):
            print(f"{heme_id:4d} |", end=" ")
            for j in range(n):
                formatted_value = self._format_energy_value(matrix[i][j], sequence[i], sequence[j], sequence)
                print(formatted_value, end=" ")
            print()  # New line after each row

        # Print summary statistics
        print(f"\n{COLORS['BOLD']}Summary Statistics:{COLORS['ENDC']}")
        print(f"Maximum interaction energy: {np.max(matrix_np[~np.eye(n, dtype=bool)]):.3f} eV")
        print(f"Minimum interaction energy: {np.min(matrix_np[~np.eye(n, dtype=bool)]):.3f} eV")
        print(f"Average interaction magnitude: {np.mean(np.abs(matrix_np[~np.eye(n, dtype=bool)])):.3f} eV")
        print(f"Maximum site energy: {np.max(np.diag(matrix_np)):.3f} eV")
        print(f"Minimum site energy: {np.min(np.diag(matrix_np)):.3f} eV")

    def plot_energy_matrix(self,
                        sequence: List[int],
                        energies: Dict[Tuple[int, int], float],
                        output_file: Optional[Path] = None) -> None:
        """
        Visualize the energy matrix with color-coded diagonal and off-diagonal elements.

        Args:
            sequence: List of heme IDs in order
            energies: Dictionary mapping heme pairs to their energies
            output_file: Optional path to save the plot
        """
        if plt is None:
            print("Matplotlib is not installed. Cannot create visualization.")
            return

        # Ask user about periodicity
        use_periodic = self.interaction_manager.yes_no_prompt(
            "use_periodic_numbering",
            "\nUse periodic numbering (last heme as '1'')? [If no, will use sequential numbering]"
        )

        # Ask whether to shift diagonal elements
        shift_diagonal = self.interaction_manager.yes_no_prompt(
            "make_diagonal_relative_to_min",
            "\nShould the minimum diagonal element be set as the zero reference point?"
        )

        # Convert dictionary to matrix form
        n = len(sequence)
        matrix = np.zeros((n, n))
        for i, heme1_id in enumerate(sequence):
            for j, heme2_id in enumerate(sequence):
                if heme1_id <= heme2_id:
                    matrix[i, j] = energies[(heme1_id, heme2_id)]
                else:
                    matrix[i, j] = energies[(heme2_id, heme1_id)]

        # Convert to meV and store original min_diagonal if needed
        matrix = matrix * 1000 #convert to meV 
        min_diagonal = np.min(np.diag(matrix)) if shift_diagonal else None

        # Shift diagonal if requested
        if shift_diagonal:
            matrix = matrix.copy()
            np.fill_diagonal(matrix, matrix.diagonal() - min_diagonal)

        # Create masks for diagonal and off-diagonal elements
        diagonal_mask = ~np.eye(n, dtype=bool)
        off_diagonal_mask = np.eye(n, dtype=bool)

        # Calculate value ranges for color scales
        diagonal_elements = np.diag(matrix)
        off_diagonal_elements = matrix[~np.eye(n, dtype=bool)]
        vmin_diag, vmax_diag = np.min(diagonal_elements), np.max(diagonal_elements)
        vmin_off, vmax_off = np.min(off_diagonal_elements), np.max(off_diagonal_elements)

        # Create figure and axis with appropriate size
        plt.figure(figsize=(3.3, 3.3), dpi=300)
        ax = plt.gca()

        # Plot off-diagonal elements
        plt.pcolormesh(np.ma.array(matrix, mask=off_diagonal_mask),
                    cmap='YlOrBr_r',
                    vmin=vmin_off,
                    vmax=vmax_off)

        # Plot diagonal elements
        plt.pcolormesh(np.ma.array(matrix, mask=diagonal_mask),
                    cmap='Blues_r',
                    vmin=vmin_diag,
                    vmax=vmax_diag)

        # Add text annotations with dynamic color
        for i in range(n):
            for j in range(n):
                val = matrix[i, j]
                is_diagonal = (i == j)

                # Determine text color based on background
                if is_diagonal:
                    if vmax_diag == vmin_diag:
                        text_color = 'white' if vmin_diag < 0 else 'black'
                    else:
                        normalized_val = (val - vmin_diag) / (vmax_diag - vmin_diag)
                        text_color = 'white' if normalized_val < 0.5 else 'black'
                else:
                    if vmax_off == vmin_off:
                        text_color = 'white' if vmax_off < 0 else 'black'
                    else:
                        normalized_val = (val - vmin_off) / (vmax_off - vmin_off)
                        text_color = 'white' if normalized_val < 0.6 else 'black'

                plt.text(j + 0.5, i + 0.5, f'{val:.0f}',
                        ha='center', va='center',
                        color=text_color, fontsize=8)

        # Create tick labels based on user choice
        if use_periodic:
            tick_labels = [str(i) for i in sequence[:-1]]
            tick_labels.append("1'")  # Replace last number with "1'" for periodic
        else:
            tick_labels = [str(i) for i in sequence]  # Use sequential numbering

        # Set ticks and labels
        plt.xticks(np.arange(n) + 0.5, tick_labels)
        plt.yticks(np.arange(n) + 0.5, tick_labels)

        # Add labels and title
        plt.xlabel('Heme Index')
        plt.ylabel('Heme Index')
        title = f'Energy Matrix (meV)'
        if shift_diagonal:
            title += f'\n(diagonal terms shifted by {min_diagonal:.2f})'
        plt.title(title)

        # Make plot square 
        ax.set_aspect('equal')

        # Adjust layout and save/show
        plt.tight_layout()
        
        if output_file:
            # Ensure the parent directory exists
            output_file.parent.mkdir(parents=True, exist_ok=True)
            # Save with high DPI and tight borders
            plt.savefig(str(output_file), dpi=300, bbox_inches='tight')
            print(f"\nPlot saved to: {output_file}")
        else:
            plt.show()
        
        plt.close()

    def compute_heme_interactions(self,
                                sequence: List[int],
                                dielectric_constants: List[float],
                                n_parallel: Optional[int] = None,
                                plot: bool = True,
                                shift_diagonal: bool = False
    ) -> Dict[Tuple[int, int], float]:

        # First check for existing energy matrix
        matrix_file = self.ee_dir / "EnergyMatrix.txt"

        if matrix_file.exists():
            use_existing = self.interaction_manager.yes_no_prompt(
                "use_existing_interactions",
                "\nFound existing energy matrix. Use these values?"
            )
            if use_existing:
                energies = self.read_and_display_energy_matrix()

                print("\n Plotting Energy Matrix...")
                output_path = self.ee_dir / "EnergyMatrix.png" 
                self.plot_energy_matrix(sequence, energies, output_path)

                if energies:
                    return energies
                print("\nWarning: Could not read existing energy matrix.")

        # Check existing PBSA calculations and prepare needed ones
        print(f"\nAnalyzing calculations needed for {len(sequence)} hemes...")
        pairs = [(i, j) for i in sequence for j in sequence if i <= j]

        energies = {}
        existing_calcs = 0
        needed_calcs = []

        for heme1_id, heme2_id in pairs:
            if heme1_id == heme2_id:
                # Diagonal elements - need ox and red states
                for state in ['o', 'r']:
                    output_file = self.ee_dir / f"pbsa_{state}_{heme1_id}_{heme1_id}.out"
                    if output_file.exists():
                        existing_calcs += 1
                    else:
                        needed_calcs.append((heme1_id, heme1_id, state, None))
            else:
                # Off-diagonal elements - need all four states
                for state in ['oo', 'or', 'ro', 'rr']:
                    output_file = self.ee_dir / f"pbsa_{state}_{heme1_id}_{heme2_id}.out"
                    if output_file.exists():
                        existing_calcs += 1
                    else:
                        needed_calcs.append((heme1_id, heme2_id, state, None))

        total_calcs = sum(4 if i != j else 2 for i, j in pairs)

        # Report status
        if not needed_calcs:
            print(f"\nAll {total_calcs} PBSA calculations already exist.")
        else:
            print(
                f"\nFound {existing_calcs} existing calculations. "
                f"Need to run {len(needed_calcs)} new calculations."
            )

        # Get reference state choice
        ref_state = RedoxState.OXIDIZED if self.interaction_manager.prompt(
            "ref_state",
            "\nShould the reference state be all-hemes oxidized (ox) "
            "or all-hemes reduced (red)? ",
            choices=['ox', 'red']
        ) == 'ox' else RedoxState.REDUCED

        # Initialize state manager and generate needed states
        state_manager = RedoxStateManager(
            input_pdb=self.pdb_file,
            heme_ids=sequence,
            launch_dir=self.launch_dir,
            forcefield_dir=self.forcefield_dir,
            reference_state=ref_state
        )

        # Generate all needed states
        state_manager.generate_all_states()

        # Compute local dielectric constants for all hemes
        local_dielectrics = self._compute_local_dielectrics(sequence, dielectric_constants)

        # Setup parallel processing if requested
        if n_parallel is None and needed_calcs:
            parallel = self.interaction_manager.yes_no_prompt(
                "run_parallel",
                "\nRun calculations in parallel?"
            )
            if parallel:
                n_parallel = self.interaction_manager.prompt(
                    "n_parallel",
                    "How many calculations to run in parallel? (Enter for maximum): ",
                    input_type=int,
                    allow_empty=True
                ) or len(pairs)

        # Add PBSA parameters to needed calculations
        for i, calc in enumerate(needed_calcs):
            heme1_id, heme2_id, state, _ = calc
            params = self._get_pbsa_parameters(heme1_id, heme2_id, local_dielectrics)
            needed_calcs[i] = (heme1_id, heme2_id, state, params)

        # Run needed calculations
        if needed_calcs:
            if n_parallel and n_parallel > 1:
                self._run_parallel_pbsa(needed_calcs, n_parallel)
            else:
                self._run_serial_pbsa(needed_calcs)

        # Process results and generate energies dictionary
        for heme1_id, heme2_id in pairs:
            if heme1_id == heme2_id:
                # Calculate site energy (E_ox - E_red)
                ox_energy = self._extract_pbsa_energy(
                    self.ee_dir / f"pbsa_o_{heme1_id}_{heme1_id}.out")
                red_energy = self._extract_pbsa_energy(
                    self.ee_dir / f"pbsa_r_{heme1_id}_{heme1_id}.out")
                energies[(heme1_id, heme2_id)] = ox_energy - red_energy
            else:
                # Calculate interaction energy (E_OO - E_OR - E_RO + E_RR)
                oo_energy = self._extract_pbsa_energy(
                    self.ee_dir / f"pbsa_oo_{heme1_id}_{heme2_id}.out")
                or_energy = self._extract_pbsa_energy(
                    self.ee_dir / f"pbsa_or_{heme1_id}_{heme2_id}.out")
                ro_energy = self._extract_pbsa_energy(
                    self.ee_dir / f"pbsa_ro_{heme1_id}_{heme2_id}.out")
                rr_energy = self._extract_pbsa_energy(
                    self.ee_dir / f"pbsa_rr_{heme1_id}_{heme2_id}.out")
                energies[(heme1_id, heme2_id)] = (
                    oo_energy - or_energy - ro_energy + rr_energy)

        # Write results to file
        self._write_results(sequence, energies)

        # Display the matrix
        self.print_energy_matrix(sequence, energies)

        # Plot the matrix if requested
        if plot:
            print("\n Plotting Energy Matrix...")
            # Use default or provided plot name
            output_path = self.ee_dir / "EnergyMatrix.png" 
            self.plot_energy_matrix(sequence, energies, output_path)

        return energies

    def _parse_lambda_file(self, sequence: List[int]) -> List[float]:
        """
        Parse Lambda.txt file to extract dielectric constants for hemes.

        Args:
            sequence: List of heme IDs to process

        Returns:
            List of dielectric constants matching the sequence order
        """
        lambda_file = self.ee_dir / "Lambda.txt"

        if not lambda_file.exists():
            print("No Lambda.txt file found in the EE directory.")
            return []

        dielectrics = {}  # Store all heme pairs and their Es values
        current_pair = None

        with open(lambda_file, 'r') as f:
            for line in f:
                line = line.strip()

                # Check for new heme pair
                pair_match = re.match(r'HEM-(\d+)\s*->\s*HEM-(\d+)', line)
                if pair_match:
                    current_pair = (int(pair_match.group(1)), int(pair_match.group(2)))
                    continue

                # Look for Es value if we have a current pair
                if current_pair and 'Es' in line:
                    es_match = re.search(r'Es\s*=\s*([\d.]+)', line)
                    if es_match:
                        # Store the Es value for both hemes in the pair
                        es_value = float(es_match.group(1))
                        dielectrics[current_pair[0]] = es_value
                        dielectrics[current_pair[1]] = es_value

        # Build the result list following the sequence order
        result = []
        for heme_id in sequence:
            if heme_id in dielectrics:
                result.append(dielectrics[heme_id])
            else:
                print(f"Warning: No Es value found for HEM-{heme_id}")
                result.append(0.0)  # or some default value

        return result

    def _compute_local_dielectrics(self,
                                    sequence: List[int],
                                    pair_dielectrics: Optional[List[float]] = None
    ) -> Dict[Tuple[int, int], float]:
        """
        Compute local dielectric constants for all heme pairs, with user confirmation.
        Handles cases where dielectric constants are not available.
        """
        # If no pair dielectrics are provided, try to parse from Lambda.txt first
        if not pair_dielectrics:
            lambda_file = self.ee_dir / "Lambda.txt"
            if lambda_file.exists():
                use_lambda = self.interaction_manager.yes_no_prompt(
                    "use_lambda",
                    "\nFound Lambda.txt file. Would you like to use dielectric constants from it?"
                )
                if use_lambda:
                    pair_dielectrics = self._parse_lambda_file(sequence)
        
            # If still no pair dielectrics (no file or user declined), prompt user
            if not pair_dielectrics:
                print("\nNo dielectric constants available from previous calculation.")

                # Ask user for dielectric constant for each heme
                heme_dielectrics = {}
                for heme_id in sequence:
                    manual_eps = self.interaction_manager.prompt(
                        f"manual_eps_{heme_id}",
                        f"Enter dielectric constant for heme {heme_id}: ",
                        input_type=float
                    )
                    heme_dielectrics[heme_id] = manual_eps
            
                # Compute pair dielectrics as average of adjacent hemes
                pair_eps = {}
                for i in sequence:
                    for j in sequence:
                        if i <= j:
                            if i == j:
                                pair_eps[(i, j)] = heme_dielectrics[i]
                            else:
                                pair_eps[(i, j)] = (heme_dielectrics[i] + heme_dielectrics[j]) / 2
            
                return pair_eps

        # Logic for when pair dielectrics are available
        # First compute local average for each heme
        heme_dielectrics = {}
        for idx, heme_id in enumerate(sequence):
            local_eps = []
            
            # Add dielectric from previous pair if it exists
            if idx > 0:
                local_eps.append(pair_dielectrics[idx-1])
                
            # Add dielectric from next pair if it exists
            if idx < len(sequence) - 1:
                local_eps.append(pair_dielectrics[idx])
                
            # Compute local average and ask user
            suggested_eps = sum(local_eps) / len(local_eps)
            
            use_suggested = self.interaction_manager.yes_no_prompt(
                f"use_eps_{heme_id}",
                f"\nUse suggested dielectric constant {suggested_eps:.3f} for heme {heme_id}?"
            )
            
            if use_suggested:
                heme_dielectrics[heme_id] = suggested_eps
            else:
                manual_eps = self.interaction_manager.prompt(
                    f"manual_eps_{heme_id}",
                    f"Enter dielectric constant for heme {heme_id}: ",
                    input_type=float
                )
                heme_dielectrics[heme_id] = manual_eps

        # Compute pair dielectrics as average of adjacent hemes
        pair_eps = {}
        for i in sequence:
            for j in sequence:
                if i <= j:
                    if i == j:
                        pair_eps[(i, j)] = heme_dielectrics[i]
                    else:
                        # Show and confirm pair dielectric
                        suggested_eps = (heme_dielectrics[i] + heme_dielectrics[j]) / 2
                        use_suggested = self.interaction_manager.yes_no_prompt(
                            f"use_eps_{i}_{j}",
                            f"\nUse suggested dielectric constant {suggested_eps:.3f} for pair {i}-{j}?"
                        )
                        
                        if use_suggested:
                            pair_eps[(i, j)] = suggested_eps
                        else:
                            manual_eps = self.interaction_manager.prompt(
                                f"manual_eps_{i}_{j}",
                                f"Enter dielectric constant for pair {i}-{j}: ",
                                input_type=float
                            )
                            pair_eps[(i, j)] = manual_eps
        
        return pair_eps

    def _get_pbsa_parameters(
        self, 
        heme1_id: int, 
        heme2_id: int, 
        local_dielectrics: Dict[Tuple[int, int], float]
    ) -> PBSAParameters:
        """
        Get PBSA parameters with focused user input.
        
        Args:
            heme1_id: First heme ID
            heme2_id: Second heme ID
            local_dielectrics: Dictionary mapping heme pairs to their dielectric constants
            
        Returns:
            PBSAParameters object
        """
        # Global parameter file for consistent settings
        global_params_file = self.ee_dir / "global_pbsa_params.txt"
        
        # Get pair-specific dielectric constant
        epsin = local_dielectrics[(min(heme1_id, heme2_id), max(heme1_id, heme2_id))]
        
        # Default parameters (these will be set once and reused)
        default_params = {
            'istrng': 150.0,
            'membraneopt': 0,
            'epsmem': 1.0,
            'mthick': 40.0,
            'poretype': 0,
            'ipb': 2,
            'inp': 2,
            'ivalence': 0,
            'bcopt': 5,
            'eneopt': 2,
            'maxitn': 100,
            'nfocus': 2,
            'smoothopt': 1
        }
        
        # Try to load global parameters if they exist
        if global_params_file.exists():
            try:
                with open(global_params_file, 'r') as f:
                    saved_global_params = eval(f.read())
                default_params.update(saved_global_params)
            except Exception as e:
                print(f"Warning: Could not read global parameters: {e}")
        else:
            # First-time setup: prompt for global parameters
            print("\n--- Initial PBSA Global Parameter Setup ---")
            for param, value in default_params.items():
                new_value = self.interaction_manager.prompt(
                    f"global_{param}",
                    f"\nEnter value for {param}",  # Simplified message
                    input_type=type(value),
                    default=value,  # Set the default value
                    allow_empty=True
                )
                
                # Update parameter if a new value was provided
                if new_value is not None:
                    default_params[param] = new_value

            # Save global parameters
            with open(global_params_file, 'w') as f:
                f.write(str(default_params))
       
        # Handle external dielectric constant
        epsout = self.interaction_manager.prompt(
            f"epsout_{heme1_id}_{heme2_id}",
            f"\nExternal dielectric constant for heme pair {heme1_id}-{heme2_id}",  # Remove the colon
            input_type=float,
            default=default_params.get('epsout', 80.0),  # Set explicit default
            allow_empty=True
        )

        # Save epsout to global params if not already present
        if 'epsout' not in default_params:
            default_params['epsout'] = epsout
            with open(global_params_file, 'w') as f:
                f.write(str(default_params))
        
        # Create PBSAParameters
        return PBSAParameters(
            epsin=epsin,
            epsout=epsout,
            istrng=default_params['istrng'],
            membraneopt=default_params['membraneopt'],
            epsmem=default_params['epsmem'],
            mthick=default_params['mthick'],
            poretype=default_params['poretype'],
            ipb=default_params['ipb'],
            inp=default_params['inp'],
            ivalence=default_params['ivalence'],
            bcopt=default_params['bcopt'],
            eneopt=default_params['eneopt'],
            maxitn=default_params['maxitn'],
            nfocus=default_params['nfocus'],
            smoothopt=default_params['smoothopt']
        )

    def _generate_pbsa_input(self, 
                            heme1_id: int,
                            heme2_id: int,
                            state: str,
                            params: PBSAParameters) -> Path:
        """Generate PBSA input file."""
        input_file = self.ee_dir / f"pbsa_{state}_{heme1_id}_{heme2_id}.in"
        
        with open(input_file, 'w') as f:
            f.write(f"""# Single point PB calculation for heme pair {heme1_id}-{heme2_id} ({state})
    &cntrl
    ipb={params.ipb},        ! PB method option
    inp={params.inp},        ! Non-polar solvation option
    ntx=1,            ! Read coordinates only
    imin=1,           ! Single-point energy
    /

    &pb
    pbtemp=300,       ! Temperature
    ivalence={params.ivalence},     ! Ion valence option
    istrng={params.istrng},      ! Ionic strength (mM)
    epsin={params.epsin},       ! Internal dielectric
    epsout={params.epsout},      ! External dielectric
    epsmem={params.epsmem},      ! Membrane dielectric
    membraneopt={params.membraneopt},    ! Membrane present
    mthick={params.mthick},   ! Membrane thickness
    mctrdz=0,         ! Membrane center
    poretype={params.poretype},       ! Pore detection
    radiopt=0,        ! Use topology radii
    dprob=1.4,        ! Solvent probe radius
    iprob=2.0,        ! Ion probe radius
    sasopt=0,         ! Surface calculation option
    bcopt={params.bcopt},          ! Boundary condition option
    eneopt={params.eneopt},         ! Energy calculation option
    maxitn={params.maxitn},       ! Maximum iterations
    nfocus={params.nfocus},         ! Focusing steps
    fscale=8,         ! Focus scaling
    smoothopt={params.smoothopt},      ! Smoothing option
    /
    """)
        return input_file

    def _run_parallel_pbsa(self,
                            calculations: List[Tuple[int, int, str, PBSAParameters]], 
                            n_parallel: int):
        """Run PBSA calculations in parallel."""
        from math import ceil
        
        total = len(calculations)
        batch_size = min(n_parallel, total)
        n_batches = ceil(total / batch_size)
        
        print(
            f"\nRunning {total} calculations in {n_batches} batches "
            f"of {batch_size}...",
        )
        
        for batch_num in range(n_batches):
            start = batch_num * batch_size
            end = min(start + batch_size, total)
            batch = calculations[start:end]
            
            # Start all processes in this batch
            processes = []
            for heme1_id, heme2_id, state, params in batch:
                print(
                    f"Starting PBSA for pair {heme1_id}-{heme2_id} {state} state...",
                )
                
                # Generate input file
                input_file = self._generate_pbsa_input(heme1_id, heme2_id, state, params)
                output_file = self.ee_dir / f"pbsa_{state}_{heme1_id}_{heme2_id}.out"
                
                cmd = (f"pbsa -O -i {input_file} -o {output_file} "
                    f"-p {state}_{heme1_id}_{heme2_id}.prmtop "
                    f"-c {state}_{heme1_id}_{heme2_id}.rst7")
                
                processes.append(subprocess.Popen(cmd, shell=True))
            
            # Wait for all processes in this batch
            for p in processes:
                p.wait()
                if p.returncode != 0:
                    raise RuntimeError(f"PBSA calculation failed with code {p.returncode}")

    def _run_serial_pbsa(self,
                        calculations: List[Tuple[int, int, str, PBSAParameters]]):
        """Run PBSA calculations serially."""
        for heme1_id, heme2_id, state, params in calculations:
            print(
                f"\nRunning PBSA for state {state} of pair {heme1_id}-{heme2_id} ...",
            )
            
            # Generate input file
            input_file = self._generate_pbsa_input(heme1_id, heme2_id, state, params)
            output_file = self.ee_dir / f"pbsa_{state}_{heme1_id}_{heme2_id}.out"
            
            cmd = (f"pbsa -O -i {input_file} -o {output_file} "
                f"-p {state}_{heme1_id}_{heme2_id}.prmtop "
                f"-c {state}_{heme1_id}_{heme2_id}.rst7")
            
            try:
                subprocess.run(cmd, shell=True, check=True)
            except subprocess.CalledProcessError as e:
                raise RuntimeError(f"PBSA calculation failed for pair "
                                f"{heme1_id}-{heme2_id} {state} state: {e}")
      
    def _check_files_exist(self, heme1_id: int, heme2_id: int, state: str) -> Dict[str, bool]:
        """
        Check if files exist for a given heme pair and state.
        
        Args:
            heme1_id: First heme ID
            heme2_id: Second heme ID
            state: Redox state (OO, OR, RO, or RR)
            
        Returns:
            Dictionary of file existence status
        """
        return {
            'prmtop': (self.ee_dir / f"{state}_{heme1_id}_{heme2_id}.prmtop").exists(),
            'rst7': (self.ee_dir / f"{state}_{heme1_id}_{heme2_id}.rst7").exists(),
            'pbsa_out': (self.ee_dir / f"pbsa_{state}_{heme1_id}_{heme2_id}.out").exists()
        }
    
    def _extract_pbsa_energy(self, output_file: Path) -> float:
        """Extract total energy from PBSA output file."""
        try:
            with open(output_file, 'r', encoding='utf-8') as f:
                for line in f:
                    if 'Etot' in line:
                        return float(line.split()[2]) * 0.043  # Convert to eV
        except (FileNotFoundError, ValueError, IndexError) as e:
            raise RuntimeError(f"Failed to extract energy from {output_file}: {e}")

    def _read_existing_results(self, file_path: Path) -> Dict[Tuple[int, int], float]:
        """
        Read energies from existing StateEnergies.txt.
        
        Args:
            file_path: Path to existing results file
            
        Returns:
            Dictionary mapping heme pairs to their energies
            (includes both site energies and interaction energies)
        """
        energies = {}
        try:
            with open(file_path) as f:
                # Find matrix section
                for line in f:
                    if line.strip() == "Complete Energy Matrix (meV):":
                        # Next line is header with heme IDs
                        header = next(f).strip().split()[1:]  # Skip first empty space
                        heme_ids = [int(x) for x in header]
                        
                        # Read matrix
                        for i, line in enumerate(f):
                            if not line.strip():  # Empty line
                                break
                                
                            # Parse row
                            values = [float(x)/1000.0 for x in line.strip().split()[1:]]  # Skip row label
                            heme1_id = heme_ids[i]
                            
                            for j, value in enumerate(values):
                                heme2_id = heme_ids[j]
                                if heme1_id <= heme2_id:  # Only store upper triangle and diagonal
                                    energies[(heme1_id, heme2_id)] = value
                        break
                        
        except (FileNotFoundError, ValueError, IndexError, StopIteration):
            return {}
            
        return energies

    def _write_results(self, 
                    sequence: List[int],
                    energies: Dict[Tuple[int, int], float]) -> None:
        """
        Write results to StateEnergies.txt.
        
        Args:
            sequence: List of heme IDs in order
            energies: Dictionary mapping heme pairs to their energies
                (includes both site energies on diagonal and interaction energies off-diagonal)
        """
        # Generate matrix first (we'll need it for both files)
        n = len(sequence)
        matrix = [[0.0 for _ in range(n)] for _ in range(n)]
        
        for i, heme1_id in enumerate(sequence):
            for j, heme2_id in enumerate(sequence):
                if heme1_id <= heme2_id:
                    # Direct from energies dict
                    value = energies[(heme1_id, heme2_id)]
                else:
                    # Symmetric for interaction energies
                    value = energies[(heme2_id, heme1_id)]
                matrix[i][j] = value  # Store in eV
            
        with open(self.ee_dir / "StateEnergies.txt", 'w') as f:
            f.write("Heme-Heme Energetics\n\n")
            
            # Write method description
            f.write("Energies calculated using PBSA method:\n")
            f.write("Diagonal elements (site energies):\n")
            f.write("  E_site = E_ox - E_red\n")
            f.write("Off-diagonal elements (interaction energies):\n")
            f.write("  E_int = E_OO - E_OR - E_RO + E_RR\n")
            f.write("where O = oxidized, R = reduced\n\n")
            
            # Write site energies
            f.write("Site Energies:\n")
            for heme_id in sequence:
                energy = energies[(heme_id, heme_id)]
                f.write(f"Heme {heme_id}: {energy:.3f} eV\n")
            
            # Write pair-wise interaction energies
            f.write("\nPair-wise Interaction Energies:\n")
            for (heme1_id, heme2_id), energy in energies.items():
                if heme1_id != heme2_id:
                    f.write(f"Heme {heme1_id} - Heme {heme2_id}: {energy:.3f} eV\n")
                        
            f.write("\nComplete Energy Matrix (meV):\n")
            # Write heme IDs as column headers
            f.write("      " + "".join(f"{heme_id:8d}" for heme_id in sequence) + "\n")
            
            # Write matrix with row labels
            for i, heme_id in enumerate(sequence):
                f.write(f"{heme_id:4d}  ")
                f.write(" ".join(f"{value:8.3f}" for value in matrix[i]))
                f.write("\n")

        # Write computer-parsable matrix
        with open(self.ee_dir / "EnergyMatrix.txt", 'w') as f:
            # First line: heme IDs
            f.write(" ".join(str(heme_id) for heme_id in sequence) + "\n")
            # Matrix values (in eV)
            for i in range(n):
                f.write(" ".join(f"{value:.6f}" for value in matrix[i]) + "\n")

