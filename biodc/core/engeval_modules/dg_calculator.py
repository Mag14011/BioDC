"""
Module for calculating reaction free energy using PBSA method.
Supports both automated PBSA calculations and manual entry of values.
"""

import os
import re
import logging
import subprocess
from enum import Enum
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
from dataclasses import dataclass

from biodc.utils.interaction import InteractionManager
from biodc.utils.state_selector import RedoxStateManager, RedoxState

logger = logging.getLogger(__name__)

@dataclass
class PBSAParameters:
    """Container for PBSA calculation parameters."""
    epsin: float         # Internal dielectric constant
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

class CalculationType(Enum):
    """Types of PBSA calculations."""
    STANDARD = "standard"
    DELPHI = "delphi"
    MEMBRANE = "membrane"

class DeltaGCalculator:
    """Calculator for reaction free energy using PBSA method."""
    
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
        self.pdb_file = str(Path(pdb_file))  # Convert to string after validating
        self.ee_dir = launch_dir / "EE"
        self.ee_dir.mkdir(exist_ok=True)

    def compute_reaction_free_energy(
        self,
        sequence: List[int],
        dielectric_constants: Optional[List[float]] = None,
        is_cyclic: bool = False,
        n_parallel: Optional[int] = None
    ) -> List[float]:
        """
        Compute reaction free energy for sequence of hemes.
        Handles both multi-heme and single-heme cases.
        
        Args:
            sequence: List of heme residue IDs
            dielectric_constants: Optional list of dielectric constants from lambda calculation
            is_cyclic: Whether to treat sequence as cyclic
            n_parallel: Number of parallel PBSA calculations (None for serial)
            
        Returns:
            List of reaction free energies for each transfer step, or a list with a single
            value representing the redox potential for a single heme
        """
        # Check if we're dealing with a single heme
        is_single_heme = len(sequence) == 1
        
        if is_single_heme:
            print(f"\nSingle heme detected (ID: {sequence[0]})")
            print("Computing redox potential for this heme...")
            
        # Method selection
        method = self.interaction_manager.prompt(
            "dg_method",
            "\nSelect method for DG calculation:\n"
            "1) Compute DG with PBSA approach\n"
            "2) Enter DG values manually\n"
            "Choice: ",
            choices=['1', '2']
        )
        
        if method == '2':
            if is_single_heme:
                return self._get_manual_redox_potential(sequence[0])
            else:
                return self._get_manual_dg_values(sequence, is_cyclic)
            
        # Get reference state
        ref_state = self._get_reference_state()
        
        # Initialize state manager and generate files
        state_manager = self._setup_state_manager(sequence, ref_state)
        
        # Process dielectric constants
        heme_dielectrics = self._process_dielectric_constants(
            sequence, dielectric_constants, is_cyclic)
            
        # Get calculation type and parameters
        calc_type = self._get_calculation_type()
        
        # Run calculations
        energies = self._run_pbsa_calculations(
            sequence, heme_dielectrics, calc_type, n_parallel)
            
        # Compute and return DG values
        if is_single_heme:
            results = self._compute_single_heme_redox_potential(energies, sequence[0])
        else:
            results = self._compute_dg_values(energies, sequence, is_cyclic)

        # Display results automatically
        if is_single_heme:
            self.display_single_heme_results(sequence[0], results[0])
        else:
            self.display_results(sequence, results, is_cyclic)

        return results

    def _get_manual_redox_potential(self, heme_id: int) -> List[float]:
        """
        Get manually entered redox potential for a single heme.
        
        Args:
            heme_id: ID of the single heme
            
        Returns:
            List containing the single redox potential value
        """
        value = self.interaction_manager.prompt(
            f"redox_potential_manual_{heme_id}",
            f"Enter redox potential (eV) for HEM-{heme_id}: ",
            input_type=float
        )
        
        # Write to DG.txt
        with open(self.ee_dir / 'DG.txt', 'w') as f:
            f.write(f"Manual entry: HEM-{heme_id} redox potential = {value:.3f} eV\n")
        
        return [value]

    def _compute_single_heme_redox_potential(self, 
                                            energies: Dict[int, Tuple[float, float]],
                                            heme_id: int) -> List[float]:
        """
        Compute redox potential for a single heme.
        
        Args:
            energies: Dict mapping heme_id to (E_ox, E_red) tuple
            heme_id: ID of the single heme
            
        Returns:
            List containing the single redox potential value
        """
        # Calculate energy difference between oxidized and reduced states
        # E_oxidized - E_reduced
        redox_potential = energies[heme_id][0] - energies[heme_id][1]
        
        # Write to DG.txt
        with open(self.ee_dir / 'DG.txt', 'w') as f:
            f.write(f"HEM-{heme_id} redox potential = {redox_potential:.3f} eV\n")
        
        return [redox_potential]

    def display_single_heme_results(self, heme_id: int, redox_potential: float):
        """
        Display the redox potential results for a single heme.
        
        Args:
            heme_id: ID of the single heme
            redox_potential: Calculated redox potential value
        """
        print("\n" + "=" * 50)
        print(f"Redox Potential for Heme {heme_id}")
        print("=" * 50)
        print(f"\nHEM-{heme_id} redox potential = {redox_potential:.3f} eV")
        print("\nThis value represents the energy difference between")
        print("the oxidized and reduced states (E_oxidized - E_reduced).")
        print("\nResult also saved in EE/DG.txt")

    def _get_reference_state(self) -> RedoxState:
        """Get user's choice of reference state."""
        choice = self.interaction_manager.prompt(
            "ref_state",
            "\nShould the reference state be all-hemes oxidized (ox) "
            "or all-hemes reduced (red)? ",
            choices=['ox', 'red']
        )
        return RedoxState.OXIDIZED if choice == 'ox' else RedoxState.REDUCED

    def _setup_state_manager(
        self, 
        sequence: List[int],
        ref_state: RedoxState
    ) -> RedoxStateManager:
        """Initialize state manager and generate needed files."""
        state_manager = RedoxStateManager(
            input_pdb=str(self.launch_dir / self.pdb_file,),
            heme_ids=sequence,
            launch_dir=self.launch_dir,
            forcefield_dir=self.forcefield_dir,
            reference_state=ref_state
        )
        
        # Generate reference state and single heme files
        state_manager.generate_all_states()
        return state_manager

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
                        dielectrics[current_pair[1]] = es_value  # Add this line to store for second heme
    
        # Build the result list following the sequence order
        result = []
        for heme_id in sequence:
            if heme_id in dielectrics:
                result.append(dielectrics[heme_id])
            else:
                print(f"Warning: No Es value found for HEM-{heme_id}")
                result.append(0.0)  # or some default value
            
        return result

    def _process_dielectric_constants(
        self,
        sequence: List[int],
        pair_dielectrics: Optional[List[float]],
        is_cyclic: bool
    ) -> Dict[int, float]:
        """
        Process dielectric constants for each heme.

        Args:
            sequence: List of heme IDs
            pair_dielectrics: Optional list of dielectric constants from lambda calculation
            is_cyclic: Whether sequence is cyclic

        Returns:
            Dict mapping heme ID to its dielectric constant
        """
        # First, check for dielectrics from Lambda.txt if no pair_dielectrics
        if not pair_dielectrics:
            pair_dielectrics = self._parse_lambda_file(sequence)

        heme_dielectrics = {}

        for idx, heme_id in enumerate(sequence):
            # Determine dielectric constant logic
            if idx == 0:
                # First heme - use first pair dielectric
                eps = pair_dielectrics[0] if pair_dielectrics else None

            elif idx == len(sequence) - 1:
                # Last heme - depends on cyclic/acyclic
                if is_cyclic and pair_dielectrics:
                    # For cyclic, average last and first pair dielectrics
                    eps = (pair_dielectrics[-1] + pair_dielectrics[0]) / 2
                else:
                    # For linear, use last available pair dielectric
                    eps = pair_dielectrics[-1] if pair_dielectrics else None

            else:
                # Middle hemes - average of adjacent pair dielectrics
                if pair_dielectrics and idx < len(pair_dielectrics):
                    eps = (pair_dielectrics[idx-1] + pair_dielectrics[idx]) / 2
                else:
                    eps = None

            # Allow user to accept/reject suggested value
            if eps is not None:
                use_suggested = self.interaction_manager.yes_no_prompt(
                    f"use_eps_{heme_id}",
                    f"\nUse suggested dielectric constant {eps:.3f} for heme-{heme_id}?"
                )
                if not use_suggested:
                    eps = None

            # If no value or rejected, prompt for manual entry
            if eps is None:
                eps = self.interaction_manager.prompt(
                    f"eps_{heme_id}",
                    f"Enter dielectric constant for heme-{heme_id}: ",
                    input_type=float
                )

            heme_dielectrics[heme_id] = eps

        return heme_dielectrics

    def _get_calculation_type(self) -> CalculationType:
        """Get user's choice of calculation type."""
        choice = self.interaction_manager.prompt(
            "pbsa_type",
            "\nSelect PBSA calculation type:\n"
            "1) Standard calculation\n"
            "2) Delphi-like calculation\n"
            "Choice: ",
            choices=['1', '2']
        )
        
        if choice == '1':
            # Check if membrane should be used
            use_membrane = self.interaction_manager.yes_no_prompt(
                "use_membrane",
                "\nShould an implicit membrane be included?"
            )
            return CalculationType.MEMBRANE if use_membrane else CalculationType.STANDARD
        else:
            return CalculationType.DELPHI

    def _get_pbsa_parameters(
        self,
        heme_id: int,
        dielectric: float,
        calc_type: CalculationType
    ) -> PBSAParameters:
        """Get PBSA parameters for a specific heme."""
        # Check if external dielectric is already in input_dict (from single heme workflow)
        if "epsout" in self.interaction_manager.input_dict:
            epsout = float(self.interaction_manager.input_dict["epsout"])
            print(f"Using previously entered external dielectric constant: {epsout}")
        else:
            # Get common parameters
            epsout = self.interaction_manager.prompt(
                f"epsout_{heme_id}",
                f"Enter external dielectric constant for heme-{heme_id}: ",
                input_type=float
            )
        
        istrng = self.interaction_manager.prompt(
            f"istrng_{heme_id}",
            "Enter ionic strength (mM): ",
            input_type=float
        )
        
        # Base parameters
        params = PBSAParameters(
            epsin=dielectric,
            epsout=epsout,
            istrng=istrng,
            membraneopt=0,
            epsmem=1.0,
            mthick=40.0,
            poretype=0
        )
        
        # Adjust based on calculation type
        if calc_type == CalculationType.DELPHI:
            params.ipb = 1
            params.inp = 0
            params.ivalence = 1
            params.bcopt = 5
            params.eneopt = 2
            params.smoothopt = 2
            params.maxitn = 100
            params.nfocus = 1
            
        elif calc_type == CalculationType.MEMBRANE:
            params.membraneopt = 1
            params.epsmem = self.interaction_manager.prompt(
                f"epsmem_{heme_id}",
                "Enter membrane dielectric constant: ",
                input_type=float
            )
            params.mthick = self.interaction_manager.prompt(
                f"mthick_{heme_id}",
                "Enter membrane thickness (Å): ",
                input_type=float
            )
            params.poretype = 1 if self.interaction_manager.yes_no_prompt(
                f"use_pore_{heme_id}",
                "Should solvent-filled channels be detected automatically?"
            ) else 0
            params.ipb = 1
            params.inp = 0
            params.bcopt = 10
            params.eneopt = 1
            params.smoothopt = 1
            params.maxitn = 200
            params.nfocus = 1
            
        return params

    def _generate_pbsa_input(
        self,
        heme_id: int,
        params: PBSAParameters
    ) -> Path:
        """Generate PBSA input file."""
        input_file = self.ee_dir / f"pbsa_{heme_id}.in"
        
        with open(input_file, 'w') as f:
            f.write(f"""# Single point PB calculation for heme {heme_id}
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

    def _run_pbsa_calculations(
        self,
        sequence: List[int],
        dielectrics: Dict[int, float],
        calc_type: CalculationType,
        n_parallel: Optional[int]
    ) -> Dict[int, Tuple[float, float]]:
        """Run PBSA calculations for all hemes."""
        # Ask about parallel execution if not specified
        if n_parallel is None:
            parallel = self.interaction_manager.yes_no_prompt(
                "run_parallel",
                "\nRun PBSA calculations in parallel?"
            )
            if parallel:
                n_parallel = self.interaction_manager.prompt(
                    "n_parallel",
                    "How many calculations to run in parallel? (Enter for maximum): ",
                    input_type=int,
                    allow_empty=True
                ) or len(sequence) * 2
        
        # Prepare all calculations
        calculations = []
        total_needed = len(sequence) * 2  # Each heme needs ox and red states
        existing_files = 0

        for heme_id in sequence:
            params = self._get_pbsa_parameters(heme_id, dielectrics[heme_id], calc_type)
            input_file = self._generate_pbsa_input(heme_id, params)
            
            for state in ['o', 'r']:
                output_file = self.ee_dir / f"pbsa_{state}_{heme_id}_{heme_id}.out"
                if output_file.exists():
                    existing_files += 1
                    if state == 'o':
                        print(
                        f"\nFound existing PBSA output for oxidized heme-{heme_id}.")
                    if state == 'r':
                        print(
                        f"\nFound existing PBSA output for reduced heme-{heme_id}.")
                        
                cmd = (f"pbsa -O -i {input_file} -o {output_file} "
                    f"-p {state}_{heme_id}_{heme_id}.prmtop -c {state}_{heme_id}_{heme_id}.rst7")
                calculations.append((heme_id, state, cmd))
        
        # Run calculations with more informative status
        if not calculations:
            print(f"\nAll {total_needed} PBSA calculations already exist. Using existing results.")
        elif existing_files > 0:
            print(
                f"\nFound {existing_files} existing calculations. "
                f"Need to run {len(calculations)} new calculations.")
            
            if n_parallel and n_parallel > 1:
                self._run_parallel(calculations, n_parallel, sequence)
            else:
                self._run_serial(calculations)
        else:
            if n_parallel and n_parallel > 1:
                self._run_parallel(calculations, n_parallel, sequence)
            else:
                self._run_serial(calculations)

        # Process all results
        energies = {}
        for heme_id in sequence:
            try:
                ox_energy = self._extract_pbsa_energy(f"pbsa_o_{heme_id}_{heme_id}.out")
                red_energy = self._extract_pbsa_energy(f"pbsa_r_{heme_id}_{heme_id}.out")
                energies[heme_id] = (ox_energy, red_energy)
            except Exception as e:
                raise RuntimeError(f"Failed to process PBSA results for heme {heme_id}: {e}")
                
        return energies

    def _run_serial(self, calculations: List[Tuple[int, str, str]]):
        """Run PBSA calculations serially."""
        for heme_id, state, cmd in calculations:
            if state == 'o':
                print(
                    f"\nRunning PBSA for oxidized heme-{heme_id} ...")
            if state == 'r':
                print(f"\nRunning PBSA for reduced heme-{heme_id} ...")
            try:
                subprocess.run(cmd, shell=True, check=True)
            except subprocess.CalledProcessError as e:
                raise RuntimeError(f"PBSA calculation failed for heme {heme_id} {state} state: {e}")

    def _run_parallel(self, 
                    calculations: List[Tuple[int, str, str]], 
                    n_parallel: int,
                    sequence: List[int]):  # Add sequence parameter
        """
        Run PBSA calculations in parallel batches.
        
        Args:
            calculations: List of (heme_id, state, command) tuples
            n_parallel: Number of parallel processes to run
            sequence: List of heme IDs to process
    """
        from math import ceil

        total = len(calculations)
        batch_size = min(n_parallel, total)
        n_batches = ceil(total / batch_size)
        
        print(f"\nRunning {total} calculations in {n_batches} batches of {batch_size}...")
        
        for batch_num in range(n_batches):
            start = batch_num * batch_size
            end = min(start + batch_size, total)
            batch = calculations[start:end]
            
            # Start all processes in this batch
            processes = []
            for heme_id, state, cmd in batch:
                if state == 'o':
                    print(f"Starting calculation for oxidized heme-{heme_id} ...")
                if state == 'r':
                    print(f"Starting calculation for reduced heme-{heme_id} ...")
                processes.append(subprocess.Popen(cmd, shell=True))
            
            # Wait for all processes in this batch
            for p in processes:
                p.wait()
                if p.returncode != 0:
                    raise RuntimeError(f"PBSA calculation failed with code {p.returncode}")
                    
        energies = {}
        
        # Process all results
        for heme_id in sequence:
            try:
                ox_energy = self._extract_pbsa_energy(f"pbsa_o_{heme_id}_{heme_id}.out")
                red_energy = self._extract_pbsa_energy(f"pbsa_r_{heme_id}_{heme_id}.out")
                energies[heme_id] = (ox_energy, red_energy)
            except Exception as e:
                raise RuntimeError(f"Failed to process PBSA results for heme {heme_id}: {e}")
                
        return energies

    def _extract_pbsa_energy(self, output_file: Path) -> float:
        """Extract total energy from PBSA output file."""
        output_file = self.ee_dir / output_file
        try:
            with open(output_file, 'r', encoding='utf-8') as f:
                for line in f:
                    if 'Etot' in line:
                        return float(line.split()[2]) * 0.043  # Convert to eV
        except (FileNotFoundError, ValueError, IndexError) as e:
            raise RuntimeError(f"Failed to extract energy from {output_file}: {e}")
            
        raise RuntimeError(f"Could not find total energy in {output_file}")

    def _compute_dg_values(self,
                          energies: Dict[int, Tuple[float, float]],
                          sequence: List[int],
                          is_cyclic: bool) -> List[float]:
        """
        Compute DG values from heme energies.
        
        Args:
            energies: Dict mapping heme_id to (E_ox, E_red) tuple
            sequence: Ordered list of heme IDs
            is_cyclic: Whether to compute DG between last and first heme
            
        Returns:
            List of DG values for each transfer step
        """
        dg_values = []
        
        # Calculate regular sequence DGs
        for i in range(len(sequence) - 1):
            donor = sequence[i]
            acceptor = sequence[i + 1]
            
            # Calculate DE = E_ox - E_red for each heme
            donor_de = energies[donor][0] - energies[donor][1]
            acceptor_de = energies[acceptor][0] - energies[acceptor][1]
            
            # DG = -(-DE_donor + DE_acceptor)
            dg = -1 * (-donor_de + acceptor_de)
            dg_values.append(dg)
            
            # Write to DG.txt
            mode = 'w' if i == 0 else 'a'
            with open(self.ee_dir / 'DG.txt', mode) as f:
                f.write(f"(HEM-{donor} = {donor_de:.3f} eV) -> "
                       f"(HEM-{acceptor} = {acceptor_de:.3f} eV); "
                       f"DG = {dg:.3f} eV\n")
        
        # Add cyclic DG if requested
        if is_cyclic:
            donor = sequence[-1]
            acceptor = sequence[0]
            donor_de = energies[donor][0] - energies[donor][1]
            acceptor_de = energies[acceptor][0] - energies[acceptor][1]
            dg = -1 * (-donor_de + acceptor_de)
            dg_values.append(dg)
            
            # Write cyclic DG to file
            with open(self.ee_dir / 'DG.txt', 'a') as f:
                f.write(f"(HEM-{donor} = {donor_de:.3f} eV) -> "
                       f"(HEM-{acceptor} = {acceptor_de:.3f} eV); "
                       f"DG = {dg:.3f} eV [Cyclic]\n")
        
        return dg_values

    def _get_manual_dg_values(self, sequence: List[int], is_cyclic: bool) -> List[float]:
        """Get manually entered DG values."""
        dg_values = []
        n_steps = len(sequence) if is_cyclic else len(sequence) - 1
        
        for i in range(n_steps):
            donor = sequence[i]
            acceptor = sequence[(i + 1) % len(sequence)]
            
            value = self.interaction_manager.prompt(
                f"dg_manual_{donor}_{acceptor}",
                f"Enter DG value (eV) for transfer {donor}->{acceptor}: ",
                input_type=float
            )
            dg_values.append(value)
            
            # Write to DG.txt
            mode = 'w' if i == 0 else 'a'
            with open(self.ee_dir / 'DG.txt', mode) as f:
                f.write(f"Manual entry: HEM-{donor} -> HEM-{acceptor}; "
                       f"DG = {value:.3f} eV{' [Cyclic]' if is_cyclic and i == n_steps-1 else ''}\n")
                
        return dg_values

    def display_results(self, sequence: List[int], dg_values: List[float], is_cyclic: bool):
        """
        Nicely format and print reaction free energy results to console.
    
        Args:
            sequence: List of heme residue IDs
            dg_values: Calculated free energy values
            is_cyclic: Whether the sequence is cyclic
        """
        # Print header
        print("\n" + "=" * 50)
        print(f"{'Cyclic' if is_cyclic else 'Linear'} Electron Transfer Reaction Free Energies")
        print("=" * 50)
    
        # Prepare transfer direction arrows
        transfer_arrows = []
        for i in range(len(sequence) - 1):
            transfer_arrows.append(f"HEM-{sequence[i]} → HEM-{sequence[i+1]}")
    
        if is_cyclic:
            transfer_arrows.append(f"HEM-{sequence[-1]} → HEM-{sequence[0]} (Cyclic Closure)")
    
        # Print transfer steps and corresponding DG values
        print("\nTransfer Steps:")
        print("-" * 50)
        max_transfer_length = max(len(arrow) for arrow in transfer_arrows)
    
        for i, (arrow, dg) in enumerate(zip(transfer_arrows, dg_values), 1):
            print(f"{i:>2}. {arrow:<{max_transfer_length+5}} ΔG = {dg:>7.3f} eV")
    
        # Compute and print summary statistics
        if len(dg_values) > 0:
            print("\nSummary Statistics:")
            print("-" * 50)
            print(f"Mean ΔG:      {sum(dg_values)/len(dg_values):7.3f} eV")
            print(f"Min ΔG:       {min(dg_values):7.3f} eV")
            print(f"Max ΔG:       {max(dg_values):7.3f} eV")
    
        print("\nResults also saved in EE/DG.txt")
