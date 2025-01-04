"""
Monte Carlo Calculator module for BioDC.
Handles parameter exploration for achieving target diffusion constants.
"""
from pathlib import Path
from typing import List, Tuple, Dict, Optional, Any, Union
from dataclasses import dataclass
import numpy as np
import time
from datetime import datetime
import multiprocessing
from rich.console import Console
from rich.table import Table

from biodc.utils.interaction import InteractionManager
from biodc.utils.structure_analyzer import PDBProcessor
from biodc.core.kineval_modules.derrida import VD

@dataclass
class ParameterSet:
    """Container for a set of electron transfer parameters."""
    couplings: List[float]      # meV
    lambdas: List[float]        # eV
    deltaG: Optional[List[float]]  # eV, None if free energy optimized
    geometry: List[str]         # 'S' or 'T' for each pair
    diffusion_coeff: float      # Physical D value (cm²/s)

@dataclass
class MonteCarloResult:
    """Container for Monte Carlo exploration results."""
    parameter_sets: List[ParameterSet]
    acceptance_rate: float
    total_attempts: int
    total_time: float
    spacing_factor: float  # For converting lattice D to physical D
    sequence: List[int]    # Heme residue IDs
    geometry: str         # String of S and T characters


class MonteCarloOptimizer:
    """Calculator for Monte Carlo parameter exploration."""
    
    # Physical constants and parameter ranges as before...
    T = 300.0
    PI = 3.141592654
    KB = 8.6173304E-5  # eV/K
    HBAR = 6.582119514E-16  # eV·s
    
    COUPLING_RANGES = {
        'S': (3, 13),    # meV
        'T': (1, 4)      # meV
    }
    LAMBDA_RANGE = (0.24, 1.1)  # eV
    DG_RANGE = (-0.3, 0.3)      # eV
    
    def __init__(self,
                interaction_manager: InteractionManager,
                launch_dir: Path,
                ee_dir: Path):
        """
        Initialize Monte Carlo optimizer.
        """
        self.interaction_manager = interaction_manager
        self.launch_dir = launch_dir
        self.ee_dir = ee_dir
        self.pdb_processor = PDBProcessor()
        
    def select_heme_sequence(self, pdb_file: str) -> List[int]:
        """
        Interactively select heme sequence from PDB file.
        
        Args:
            pdb_file: Path to PDB structure file
        
        Returns:
            Selected sequence of heme residue IDs
        """
        atoms_dict = self.pdb_processor.read_pdb_atoms(pdb_file)
        
        # Initial distance cutoff
        current_cutoff = 20.0
        self.pdb_processor.distance_cutoff = current_cutoff
        
        print("\nAutomatically detecting possible heme sequences. Please wait! ...")

        while True:
            try:
                # Try linear topology first
                linear_sequence = self.pdb_processor.detect_linear_sequence(atoms_dict)
                
                # Try branched topology
                branched_sequences = []
                try:
                    branched_sequences = self.pdb_processor.detect_branched_sequence(atoms_dict)
                except ValueError:
                    pass
                
                # Build numbered menu options
                menu_options = []
                option_num = 1
                
                # Add detected sequences with numbers
                if linear_sequence:
                    menu_options.append(f"{option_num}) Use linear sequence: {' → '.join(map(str, linear_sequence))}")
                    option_num += 1
                
                for i, branch in enumerate(branched_sequences, 1):
                    menu_options.append(f"{option_num}) Use branch {i}: {' → '.join(map(str, branch))}")
                    option_num += 1
                
                # Add manual entry and threshold options
                menu_options.extend([
                    f"{option_num}) Enter sequence manually",
                    f"{option_num + 1}) Adjust detection threshold (Current: {current_cutoff}Å)"
                ])
                
                # Create menu prompt with clear layout
                menu_prompt = "\nSelect heme sequence:\n" + "\n".join(menu_options) + "\n\nEnter choice number: "
                
                # Get user selection
                choice = self.interaction_manager.prompt(
                    "sequence_selection",
                    menu_prompt,
                    choices=[str(i) for i in range(1, len(menu_options) + 1)]
                )
                
                # Process selection based on number chosen
                choice_num = int(choice)
                selected_option = menu_options[choice_num - 1]
                
                if "linear sequence" in selected_option.lower():
                    sequence = linear_sequence
                
                elif "branch" in selected_option.lower():
                    branch_index = int(selected_option.split("branch")[1].split(":")[0].strip()) - 1
                    sequence = branched_sequences[branch_index]
                
                elif "manually" in selected_option.lower():
                    while True:
                        manual_input = self.interaction_manager.prompt(
                            "manual_sequence",
                            "Enter heme residue IDs (space-separated):",
                            input_type=str
                        )
                        
                        try:
                            sequence = [int(x) for x in manual_input.split()]
                        except ValueError:
                            print("Invalid input. Please enter integer residue IDs.")
                            continue
                        
                        # Validate all hemes exist in the structure
                        if not all(heme_id in atoms_dict for heme_id in sequence):
                            print("Error: Some heme IDs not found in the PDB structure.")
                            continue
                        
                        # Confirm and display sequence
                        confirm = self.interaction_manager.yes_no_prompt(
                            "confirm_manual_sequence",
                            f"Confirm sequence: {' → '.join(map(str, sequence))}?"
                        )
                        
                        if confirm:
                            break
                
                elif "threshold" in selected_option.lower():
                    # Get new distance threshold from user
                    current_cutoff = self.interaction_manager.prompt(
                        "distance_threshold",
                        f"Enter new distance threshold (current: {current_cutoff}Å):",
                        input_type=float
                    )
                    
                    # Update processor's distance cutoff and continue loop
                    self.pdb_processor.distance_cutoff = current_cutoff
                    continue
                
                # Structural analysis of the selected sequence
                if len(sequence) > 1:
                    console = Console()
                    table = Table(title="Heme Sequence Structural Analysis")
                    table.add_column("Heme Pair", style="cyan")
                    table.add_column("Edge-to-Edge Distance (Å)", justify="right")
                    table.add_column("Plane Angle (°)", justify="right")
                    table.add_column("Stacking Type", justify="center")
                    
                    # Analyze consecutive heme pairs
                    for i in range(len(sequence) - 1):
                        heme1 = sequence[i]
                        heme2 = sequence[i+1]
                        
                        # Calculate minimum distance
                        min_distance = self.pdb_processor.calculate_min_distance(
                            atoms_dict[heme1], 
                            atoms_dict[heme2]
                        )
                        
                        # Calculate plane angle
                        plane_angle = self.pdb_processor.calculate_plane_angle(
                            atoms_dict[heme1], 
                            atoms_dict[heme2]
                        )
                        
                        # Classify stacking
                        stacking_type = self.pdb_processor.classify_stacking(plane_angle)
                        
                        # Add row to table
                        table.add_row(
                            f"{heme1} → {heme2}", 
                            f"{min_distance:.2f}", 
                            f"{plane_angle:.2f}", 
                            stacking_type
                        )
                    
                    # Print the analysis table
                    console.print(table)
                
                return sequence
                
            except ValueError as e:
                # If both linear and branched detection fail
                print(f"Sequence detection error: {e}")
                
                # Fallback to manual entry
                fallback = self.interaction_manager.yes_no_prompt(
                    "manual_entry_fallback",
                    "Automatic sequence detection failed. Enter sequence manually?"
                )
                
                if fallback:
                    while True:
                        manual_input = self.interaction_manager.prompt(
                            "manual_sequence_fallback",
                            "Enter heme residue IDs (space-separated):",
                            input_type=str
                        )
                        
                        try:
                            sequence = [int(x) for x in manual_input.split()]
                        except ValueError:
                            print("Invalid input. Please enter integer residue IDs.")
                            continue
                        
                        # Validate all hemes exist in the structure
                        if not all(heme_id in atoms_dict for heme_id in sequence):
                            print("Error: Some heme IDs not found in the PDB structure.")
                            continue
                        
                        return sequence
                else:
                    # If user doesn't want manual entry, exit
                    raise ValueError("No valid heme sequence selected.")
                    
    def _get_geometry_sequence(self, pdb_file: str, sequence: List[int]) -> str:
        """
        Determine sequence of geometry types (S/T) from structural analysis.
        
        Args:
            pdb_file: Path to PDB file
            sequence: List of heme residue IDs
            
        Returns:
            String of S and T characters representing geometry types
        """
        atoms_dict = self.pdb_processor.read_pdb_atoms(pdb_file)
        geometry_sequence = []
        
        for i in range(len(sequence) - 1):
            heme1, heme2 = sequence[i], sequence[i+1]
            
            # Calculate plane angle
            angle = self.pdb_processor.calculate_plane_angle(
                atoms_dict[heme1], 
                atoms_dict[heme2]
            )
            
            # Map stacking type to S/T
            stacking_type = self.pdb_processor.classify_stacking(angle)
            geometry = 'S' if stacking_type == "slip-stacked" else 'T'
            geometry_sequence.append(geometry)
            
        return ''.join(geometry_sequence)
        
    def _calculate_average_spacing(self, pdb_file: str, sequence: List[int]) -> float:
        """
        Calculate average spacing between consecutive hemes in sequence.
        
        Args:
            pdb_file: Path to PDB file
            sequence: List of heme residue IDs
            
        Returns:
            Average spacing in Angstroms
        """
        distances = self.pdb_processor.measure_avg_dist(
            pdb_file=pdb_file,
            topology='linear'  # Using linear since we have a specific sequence
        )
        
        # Extract distances from results
        spacings = [d['distance'] for d in distances]
        
        if not spacings:
            raise ValueError("Could not calculate heme spacings from structure")
            
        return np.mean(spacings)

    def _get_parameters(self) -> Dict[str, Any]:
        """Get Monte Carlo parameters through interaction manager."""
        # Number of parameter sets
        num_sets = self.interaction_manager.prompt(
            "mc_num_sets",
            "Number of parameter sets to generate (1-1000):",
            input_type=int,
            choices=[str(x) for x in range(1, 1001)]
        )
        
        # Optimization approach
        free_energy_opt = self.interaction_manager.yes_no_prompt(
            "mc_free_energy_opt",
            "Use free energy optimization?",
            default=False
        )
        
        # Target diffusion coefficient if not using free energy optimization
        if not free_energy_opt:
            d_target = float(self.interaction_manager.prompt(
                "mc_d_target",
                "Target diffusion coefficient (cm²/s):",
                input_type=str
            ))
            
            d_tolerance = float(self.interaction_manager.prompt(
                "mc_d_tolerance", 
                "Relative tolerance for target matching (0-1):",
                input_type=str,
                default="0.1"
            ))
        else:
            d_target = 0.0
            d_tolerance = 0.1
            
        # Maximum attempts
        max_attempts = self.interaction_manager.prompt(
            "mc_max_attempts",
            "Maximum number of attempts:",
            input_type=int,
            default=1000000
        )
        
        # Number of processes
        max_cores = multiprocessing.cpu_count()
        num_processes = self.interaction_manager.prompt(
            "mc_num_processes",
            f"Number of processes (1-{max_cores}):",
            input_type=int,
            choices=[str(x) for x in range(1, max_cores + 1)],
            default=str(max_cores)
        )
        
        return {
            'num_sets': num_sets,
            'free_energy_opt': free_energy_opt,
            'd_target': d_target,
            'd_tolerance': d_tolerance,
            'max_attempts': max_attempts,
            'num_processes': num_processes
        }

    def _generate_parameters(self, sequence: str) -> ParameterSet:
        """Generate random parameters for each pair in the sequence."""
        num_pairs = len(sequence)
        
        couplings = [np.random.uniform(*self.COUPLING_RANGES[geom]) 
                    for geom in sequence]
        
        lambdas = [np.random.uniform(*self.LAMBDA_RANGE) 
                  for _ in range(num_pairs)]
        
        deltaG = [np.random.uniform(*self.DG_RANGE) 
                 for _ in range(num_pairs)]
        
        return ParameterSet(
            couplings=couplings,
            lambdas=lambdas,
            deltaG=deltaG,
            geometry=list(sequence),
            diffusion_coeff=0.0  # Will be set later
        )
    
    def _compute_marcus_rates(self, coupling: float, lambda_val: float,
                           deltaG: Optional[float] = None) -> Union[float, Tuple[float, float]]:
        """
        Compute Marcus electron transfer rates.
        
        Args:
            coupling: Electronic coupling in meV
            lambda_val: Reorganization energy in eV
            deltaG: Optional driving force in eV
            
        Returns:
            Either single optimized rate or (forward_rate, backward_rate) tuple
        """
        # Convert coupling from meV to eV
        coupling_ev = coupling / 1000.0
        
        prefactor = (2 * self.PI * coupling_ev**2) / (self.HBAR * 
                   np.sqrt(4 * self.PI * lambda_val * self.KB * self.T))
        
        if deltaG is None:  # free energy optimized case
            return prefactor
        
        # Forward rate
        e_act_forward = ((deltaG + lambda_val)**2) / (4 * lambda_val)
        k_forward = prefactor * np.exp(-e_act_forward / (self.KB * self.T))
        
        # Backward rate
        e_act_backward = ((-deltaG + lambda_val)**2) / (4 * lambda_val)
        k_backward = prefactor * np.exp(-e_act_backward / (self.KB * self.T))
        
        return k_forward, k_backward
    
    def _compute_diffusion_coeff(self, params: ParameterSet,
                              free_energy_opt: bool) -> float:
        """Compute diffusion coefficient for a parameter set."""
        if free_energy_opt:
            rates = [self._compute_marcus_rates(c, l)
                    for c, l in zip(params.couplings, params.lambdas)]
            V, D = VD(rates, rates)
        else:
            forward_rates = []
            backward_rates = []
            for c, l, dg in zip(params.couplings, params.lambdas, params.deltaG):
                kf, kb = self._compute_marcus_rates(c, l, dg)
                forward_rates.append(kf)
                backward_rates.append(kb)
            V, D = VD(forward_rates, backward_rates)
            
        return D if D is not None else 0.0

    def _worker_process(self, sequence: str, free_energy_opt: bool,
                      spacing_factor: float, d_target: float,
                      d_tolerance: float, worker_id: int,
                      target_sets: int, max_attempts: int) -> Tuple[List[ParameterSet], int, int]:
        """Worker function for parallel parameter exploration."""
        # Set unique random seed for this worker
        max_seed = 2**32 - 1
        base_seed = int(time.time()) % 1000000
        worker_seed = (base_seed + worker_id) % max_seed
        np.random.seed(worker_seed)
        
        local_results = []
        local_attempts = 0
        local_accepted = 0
        
        while local_attempts < max_attempts and local_accepted < target_sets:
            local_attempts += 1
            
            # Generate parameters and compute diffusion coefficient
            params = self._generate_parameters(sequence)
            D = self._compute_diffusion_coeff(params, free_energy_opt)
            
            # Scale D to physical units
            D_scaled = D * spacing_factor
            params.diffusion_coeff = D_scaled
            
            # Check if result should be accepted
            accept = False
            if free_energy_opt:
                accept = D > 0
            elif d_target > 0:
                relative_error = abs(D_scaled - d_target) / d_target
                accept = relative_error <= d_tolerance
            else:
                accept = D > 0
                
            if accept:
                local_results.append(params)
                local_accepted += 1
                
            if local_attempts % 10000 == 0:
                print(f"Worker {worker_id}: Attempts={local_attempts}, "
                      f"Accepted={local_accepted}")
                
        return local_results, local_attempts, local_accepted
    
    def _are_parameters_similar(self, params1: ParameterSet,
                               params2: ParameterSet) -> bool:
        """
        Check if two parameter sets are similar within tolerances.
        
        Uses parameter-specific tolerances:
        - Couplings: 1.0 meV
        - DeltaG: 0.025 eV
        - Lambda: 0.025 eV
        """
        coupling_tol = 1.0    # meV
        deltaG_tol = 0.025   # eV
        lambda_tol = 0.025   # eV
        
        # Compare couplings
        for v1, v2 in zip(params1.couplings, params2.couplings):
            if abs(v1 - v2) > coupling_tol:
                return False
                
        # Compare lambdas
        for l1, l2 in zip(params1.lambdas, params2.lambdas):
            if abs(l1 - l2) > lambda_tol:
                return False
                
        # Compare deltaGs if present
        if params1.deltaG is not None and params2.deltaG is not None:
            for dg1, dg2 in zip(params1.deltaG, params2.deltaG):
                if abs(dg1 - dg2) > deltaG_tol:
                    return False
                    
        return True
    
    def _remove_redundant_sets(self, parameter_sets: List[ParameterSet]) -> List[ParameterSet]:
        """Remove redundant parameter sets using KD-tree for efficiency."""
        if not parameter_sets:
            return []
            
        # Convert parameters to numpy array for KD-tree
        data = []
        for params in parameter_sets:
            values = params.couplings + params.lambdas
            if params.deltaG is not None:
                values += params.deltaG
            data.append(values)
            
        data = np.array(data)
        
        # Normalize data for better comparison
        data_normalized = np.zeros_like(data)
        for i in range(data.shape[1]):
            column = data[:, i]
            if np.any(column < 0):  # For values that can be negative
                max_abs = np.max(np.abs(column))
                if max_abs > 0:
                    data_normalized[:, i] = column / max_abs
            else:  # For strictly positive values
                if np.any(column > 0):
                    min_positive = np.min(column[column > 0])
                    column_log = np.log10(column + min_positive)
                    data_normalized[:, i] = (
                        (column_log - np.min(column_log)) / 
                        (np.max(column_log) - np.min(column_log))
                    )
        
        # Build KD-tree and find similar pairs
        tree = cKDTree(data_normalized)
        pairs = tree.query_pairs(0.001, output_type='ndarray')
        
        # Keep sets with higher diffusion coefficients
        indices_to_remove = set()
        for i, j in pairs:
            if parameter_sets[i].diffusion_coeff < parameter_sets[j].diffusion_coeff:
                indices_to_remove.add(i)
            else:
                indices_to_remove.add(j)
                
        indices_to_keep = sorted(set(range(len(parameter_sets))) - indices_to_remove)
        unique_sets = [parameter_sets[i] for i in indices_to_keep]
        
        return unique_sets
    
    def _save_parameters(self, result: MonteCarloResult, 
                      free_energy_opt: bool,
                      d_target: float,
                      output_dir: Optional[Path] = None) -> Path:
        """Save parameter sets to file."""
        if output_dir is None:
            output_dir = self.launch_dir / "KE" / "monte_carlo"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Create filename
        timestamp = datetime.now().strftime("%d-%m-%Y_%H-%M-%S")
        output_file = output_dir / f"parameters_{timestamp}.txt"
        
        with open(output_file, 'w') as f:
            # Write header information
            f.write("# Monte Carlo Parameter Exploration Results\n")
            f.write("=" * 80 + "\n\n")
            
            # Write calculation parameters
            f.write(f"# Calculation Statistics:\n")
            f.write(f"# Total attempts: {result.total_attempts:,d}\n")
            f.write(f"# Parameter sets found: {len(result.parameter_sets):,d}\n")
            f.write(f"# Acceptance rate: {result.acceptance_rate*100:.2f}%\n")
            f.write(f"# Total time: {result.total_time:.1f}s\n")
            f.write(f"# Average spacing factor: {result.spacing_factor:.2e} cm²\n")
            if d_target > 0:
                f.write(f"# Target diffusion coefficient: {d_target:.2e} cm²/s\n")
            f.write("\n")
            
            # Write headers
            headers = ["Index"]
            for i, params in enumerate(result.parameter_sets[0].geometry):
                headers.append(f"V{i+1}-{i+2}({params})")
            for i in range(len(result.parameter_sets[0].lambdas)):
                headers.append(f"L{i+1}-{i+2}")
            if not free_energy_opt:
                for i in range(len(result.parameter_sets[0].lambdas)):
                    headers.append(f"dG{i+1}-{i+2}")
            headers.append("D")
            
            # Write header line
            f.write("  ".join(f"{h:>8s}" for h in headers) + "\n")
            
            # Write parameter sets
            for i, params in enumerate(result.parameter_sets, 1):
                row = [f"{i:8d}"]
                row.extend(f"{v:8.3f}" for v in params.couplings)
                row.extend(f"{l:8.3f}" for l in params.lambdas)
                if not free_energy_opt and params.deltaG is not None:
                    row.extend(f"{dg:8.3f}" for dg in params.deltaG)
                row.append(f"{params.diffusion_coeff:8.2e}")
                f.write("  ".join(row) + "\n")
                
        return output_file
    
    def explore_parameters(self) -> MonteCarloResult:
        """
        Run Monte Carlo exploration of parameter space.
        
        Returns:
            MonteCarloResult containing exploration results
        """
        # First, get heme sequence from structure
        pdb_path = self.ee_dir / "min.pdb"
        if not pdb_path.exists():
            raise FileNotFoundError(f"Required PDB file not found: {pdb_path}")
        
        # Get sequence and determine geometry
        sequence = self.select_heme_sequence(str(pdb_path))
        geometry_sequence = self._get_geometry_sequence(str(pdb_path), sequence)
        
        print(f"\nSelected sequence: {' → '.join(map(str, sequence))}")
        print(f"Geometry sequence: {geometry_sequence}")
        
        # Calculate average spacing
        avg_spacing = self._calculate_average_spacing(str(pdb_path), sequence)
        spacing_factor = (avg_spacing * 1E-8)**2  # Convert to cm²
        
        print(f"\nStructural parameters:")
        print(f"Average heme spacing: {avg_spacing:.2f} Å")
        print(f"Spacing factor: {spacing_factor:.2e} cm²")
        
        # Get Monte Carlo parameters
        params = self._get_parameters()
        
        print(f"\nStarting Monte Carlo parameter exploration...")
        print(f"Using {params['num_processes']} processes")
        print(f"Target parameter sets: {params['num_sets']}")
        
        # Calculate parameters using multiprocessing
        with multiprocessing.Pool(params['num_processes']) as pool:
            tasks = []
            base_sets_per_worker = params['num_sets'] // params['num_processes']
            extra_sets = params['num_sets'] % params['num_processes']
            attempts_per_worker = params['max_attempts'] // params['num_processes']
            
            # Create tasks for each worker
            for i in range(params['num_processes']):
                worker_target = base_sets_per_worker + (1 if i < extra_sets else 0)
                tasks.append((
                    geometry_sequence,  # Now using determined geometry sequence
                    params['free_energy_opt'],
                    spacing_factor,
                    params['d_target'],
                    params['d_tolerance'],
                    i,
                    worker_target,
                    attempts_per_worker
                ))
            
            # Start workers and monitor progress
            async_results = [pool.apply_async(self._worker_process, t) for t in tasks]
            
            # Monitor progress
            while any(not r.ready() for r in async_results):
                time.sleep(1)
                finished = sum(1 for r in async_results if r.ready())
                current_results = sum(len(r.get()[0]) for r in async_results if r.ready())
                print(f"Progress: {finished}/{params['num_processes']} workers finished, "
                      f"{current_results}/{params['num_sets']} results found", end='\r')
            
            # Collect all results
            print("\nCollecting results from all workers...")
            all_results = []
            total_attempts = 0
            total_accepted = 0
            
            for r in async_results:
                worker_results, worker_attempts, worker_accepted = r.get()
                all_results.extend(worker_results)
                total_attempts += worker_attempts
                total_accepted += worker_accepted
        
        # Sort results by diffusion coefficient
        all_results.sort(key=lambda x: x.diffusion_coeff, reverse=True)
        
        # Remove redundant sets
        print("\nRemoving redundant parameter sets...")
        initial_sets = len(all_results)
        unique_results = self._remove_redundant_sets(all_results)
        final_sets = len(unique_results)
        
        # Calculate statistics
        elapsed_time = time.time() - start_time
        acceptance_rate = total_accepted / total_attempts if total_attempts > 0 else 0.0
        
        print(f"\nExploration complete:")
        print(f"Total attempts: {total_attempts:,d}")
        print(f"Initial accepted sets: {initial_sets:,d}")
        print(f"Unique sets after deduplication: {final_sets:,d}")
        print(f"Acceptance rate: {acceptance_rate*100:.2f}%")
        print(f"Total time: {elapsed_time:.1f}s")
        
        return MonteCarloResult(
            parameter_sets=unique_results,
            acceptance_rate=acceptance_rate,
            total_attempts=total_attempts,
            total_time=elapsed_time,
            spacing_factor=spacing_factor,
            sequence=sequence,          # Now including selected sequence
            geometry=geometry_sequence  # and geometry sequence
        ) 
