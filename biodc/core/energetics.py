"""
Energetic Evaluation module for BioDC.
Handles computation of various energetic parameters for electron transfer.
"""
import sys
from pathlib import Path
from typing import List, Optional, Tuple, Dict
from dataclasses import dataclass, field
from rich.console import Console
from rich.table import Table

from biodc.utils.interaction import InteractionManager
from biodc.utils.structure_analyzer import PDBProcessor
from biodc.core.engeval_modules.lambda_calculator import LambdaCalculator
from biodc.core.engeval_modules.dg_calculator import DeltaGCalculator
from biodc.core.engeval_modules.interaction_calculator import HemeInteractionCalculator
from biodc.core.engeval_modules.cooperativity_analyzer import HemeCooperativityAnalyzer
from biodc.core.engeval_modules.hda_calculator import CouplingCalculator
from biodc.core.engeval_modules.rate_calculator import RateCalculator

@dataclass
class EnergeticParameters:
    """Container for energetic calculation parameters."""
    lambda_values: List[float] = field(default_factory=list)
    dielectric_constants: List[float] = field(default_factory=list)
    delta_g_values: List[float] = field(default_factory=list)  # From standard DG calculation
    dg_values_dict: Dict[str, List[float]] = field(default_factory=dict)  # Named DG sets from various sources
    coupling_values: List[float] = field(default_factory=list)
    rates_forward: List[float] = field(default_factory=list)
    rates_backward: List[float] = field(default_factory=list)
    interaction_energies: Optional[List[float]] = None
    cooperativity_results: Optional[Dict] = None  # Full cooperativity results if needed

class EnergeticEvaluation:
    """
    Orchestrates the computation of energetic parameters for electron transfer.
    
    Coordinates calculation of:
    - Reorganization energy (λ)
    - Reaction free energy (ΔG)
    - Heme-heme interaction energies
    - Electronic coupling
    - Electron transfer rates
    """

    def __init__(self, 
                interaction_manager: InteractionManager,
                launch_dir: Path,
                forcefield_dir: Path,
                pdb_file: str):
        """
        Initialize the energetic evaluation workflow.
        
        Args:
            interaction_manager: Manages user interactions and input tracking
            launch_dir: Project launch directory
            forcefield_dir: Directory containing forcefield files
            pdb_file: Path to input PDB file
        """
        # Validate PDB file exists
        pdb_path = Path(pdb_file)
        if not pdb_path.exists():
            raise FileNotFoundError(f"PDB file not found: {pdb_file}")
        
        self.interaction_manager = interaction_manager
        self.launch_dir = launch_dir
        self.ee_dir = launch_dir / "EE"
        self.forcefield_dir = forcefield_dir
        self.pdb_file = str(pdb_path)
        
        # Initialize computed parameters container
        self.computed_params = EnergeticParameters()
        
        # Initialize calculators
        self.lambda_calculator = LambdaCalculator(
            interaction_manager, launch_dir)

        self.delta_g_calculator = DeltaGCalculator(
            interaction_manager, launch_dir, forcefield_dir, self.pdb_file)

        self.heme_interaction_calculator = HemeInteractionCalculator(
            interaction_manager, launch_dir, forcefield_dir, self.pdb_file)

        self.cooperativity_analyzer = HemeCooperativityAnalyzer(
            interaction_manager=interaction_manager,
            ee_dir=self.ee_dir
        )

        self.coupling_calculator = CouplingCalculator(
            interaction_manager, launch_dir)

        self.rate_calculator = RateCalculator(
            interaction_manager=interaction_manager,
            launch_dir=launch_dir
        )

    def _analyze_sequence_structure(self, 
                                sequence: List[int], 
                                pdb_file: str) -> None:
        """
        Present detailed structural analysis of the selected sequence.
        
        Args:
            sequence: Selected heme sequence
            pdb_file: Path to PDB file
        """
        processor = PDBProcessor()
        atoms_dict = processor.read_pdb_atoms(pdb_file)
        
        # Prepare console for rich formatting
        console = Console()
        
        # Create a table for sequence analysis
        table = Table(title="Heme Sequence Structural Analysis")
        table.add_column("Heme Pair", style="cyan")
        table.add_column("Edge-to-Edge Distance (Å)", justify="right")
        table.add_column("Plane Angle (°)", justify="right")
        table.add_column("Vertical Separation (Å)", justify="right")
        table.add_column("Stacking Type", justify="center")
        
        # Analyze consecutive heme pairs
        for i in range(len(sequence) - 1):
            heme1 = sequence[i]
            heme2 = sequence[i+1]
            
            # Calculate minimum distance
            min_distance = processor.calculate_min_distance(
                atoms_dict[heme1], 
                atoms_dict[heme2]
            )
            
            # Calculate plane angle and separation
            plane_angle = processor.calculate_plane_angle(
                atoms_dict[heme1], 
                atoms_dict[heme2]
            )
            
            vertical_separation = processor.calculate_plane_separation(
                atoms_dict[heme1], 
                atoms_dict[heme2]
            )
            
            # Classify stacking
            stacking_type = processor.classify_stacking(plane_angle, vertical_separation)
            
            # Add row to table
            table.add_row(
                f"{heme1} → {heme2}", 
                f"{min_distance:.2f}", 
                f"{plane_angle:.2f}",
                f"{vertical_separation:.2f}", 
                stacking_type
            )
        
        # Print the analysis table
        console.print(table)

    def select_heme_sequence(self, pdb_file: str) -> List[int]:
        """
        Interactively select heme sequence from PDB file.
        
        Args:
            pdb_file: Path to PDB structure file
        
        Returns:
            Selected sequence of heme residue IDs
        """
        import time

        start_time = time.time()
        print("\nStarting heme sequence detection process...")
        sys.stdout.flush()

        processor = PDBProcessor()

        print("About to read PDB atoms...")
        sys.stdout.flush()
        start_read = time.time()
        atoms_dict = processor.read_pdb_atoms(pdb_file)
        print(f"PDB atoms read. Time taken: {time.time() - start_read:.2f} seconds")
        sys.stdout.flush()

##########
        # Check for single heme case
        if len(atoms_dict) == 1:
            single_heme_id = next(iter(atoms_dict.keys()))
            print(f"\nDetected a single heme with ID: {single_heme_id}")
        
            # Ask user to confirm
            confirm = self.interaction_manager.yes_no_prompt(
                "confirm_single_heme",
                f"Found only one heme (ID: {single_heme_id}). Proceed with this single heme?"
            )
        
            if confirm:
                return [single_heme_id]
            else:
                # If user doesn't confirm, let them enter manually
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
                        return sequence

        # Initial distance cutoff and retry tracking
        current_cutoff = 13.0
        processor.distance_cutoff = current_cutoff
        max_retries = 5
        retry_count = 0
        tried_thresholds = set()
        
        while retry_count < max_retries:
            print(f"\nAttempting automatic detection with {current_cutoff}Å threshold...")
            sys.stdout.flush()

            start_detect = time.time()
            linear_sequence = None
            branched_sequences = []
            detection_failed = False

            try:
                # Try linear topology
                print("Detecting linear sequence...")
                sys.stdout.flush()
                linear_sequence = processor.detect_linear_sequence(atoms_dict)
                
                # Try branched topology only if linear fails
                if not linear_sequence:
                    print("Detecting branched sequences...")
                    sys.stdout.flush()
                    try:
                        branched_sequences = processor.detect_branched_sequence(atoms_dict)
                    except ValueError:
                        pass

                print(f"Sequence detection complete. Total time: {time.time() - start_time:.2f} seconds")
                sys.stdout.flush()

            except ValueError as e:
                print(f"\nAutomatic detection failed: {e}")
                detection_failed = True

            # If we found any valid sequence, show the full menu
            if linear_sequence or branched_sequences:
                menu_options = []
                option_num = 1

                if linear_sequence:
                    menu_options.append(f"{option_num}) Use linear sequence: {' → '.join(map(str, linear_sequence))}")
                    option_num += 1

                for i, branch in enumerate(branched_sequences, 1):
                    menu_options.append(f"{option_num}) Use branch {i}: {' → '.join(map(str, branch))}")
                    option_num += 1

                menu_options.extend([
                    f"{option_num}) Adjust detection threshold (Current: {current_cutoff}Å)",
                    f"{option_num + 1}) Enter sequence manually"
                ])

                # Clear cached input before showing menu
                self.interaction_manager.input_dict.pop('sequence_selection', None)
                
                # Show menu and get choice
                menu_prompt = "\nSelect heme sequence:\n" + "\n".join(menu_options) + "\n\nEnter choice number: "
                choice = self.interaction_manager.prompt(
                    "sequence_selection",
                    menu_prompt,
                    choices=[str(i) for i in range(1, len(menu_options) + 1)]
                )
                choice_num = int(choice)
                
                # Process choice
                if "linear sequence" in menu_options[choice_num - 1].lower():
                    sequence = linear_sequence
                    break
                elif "branch" in menu_options[choice_num - 1].lower():
                    branch_index = int(menu_options[choice_num - 1].split("branch")[1].split(":")[0].strip()) - 1
                    sequence = branched_sequences[branch_index]
                    break
                elif "threshold" in menu_options[choice_num - 1].lower():
                    retry_count += 1
                    if retry_count >= max_retries:
                        print("\nWarning: Multiple threshold adjustments have not helped.")
                        print("The structure may not be suitable for automatic detection.")
                        
                        # Clear cached input before prompt
                        self.interaction_manager.input_dict.pop('manual_entry_fallback', None)
                        
                        if self.interaction_manager.yes_no_prompt(
                            "manual_entry_fallback",
                            "\nWould you like to enter the sequence manually?"
                        ):
                            sequence = self._handle_manual_entry(atoms_dict)
                            if sequence:
                                break
                        else:
                            print("\nResetting retry count for one final attempt...")
                            retry_count = 0
                            tried_thresholds.clear()
                            continue
                    
                    # Get new threshold with guidance
                    current_cutoff = self._get_new_threshold(current_cutoff, tried_thresholds)
                    tried_thresholds.add(current_cutoff)
                    processor.distance_cutoff = current_cutoff
                    continue
                    
                else:  # Manual entry
                    sequence = self._handle_manual_entry(atoms_dict)
                    if sequence:
                        break

            else:  # No sequences found
                print("\nNo valid sequences found. Options:")
                print(f"1) Adjust detection threshold (Current: {current_cutoff}Å)")
                print("2) Enter sequence manually")
                
                # Clear cached input before showing options
                self.interaction_manager.input_dict.pop('failed_detection_choice', None)
                
                choice = self.interaction_manager.prompt(
                    "failed_detection_choice",
                    "\nEnter choice number: ",
                    choices=['1', '2']
                )
                
                if choice == '1':
                    retry_count += 1
                    if retry_count >= max_retries:
                        print("\nWarning: Multiple threshold adjustments have not helped.")
                        print("The structure may not be suitable for automatic detection.")
                        
                        # Clear cached input before prompt
                        self.interaction_manager.input_dict.pop('manual_entry_fallback', None)
                        
                        if self.interaction_manager.yes_no_prompt(
                            "manual_entry_fallback",
                            "\nWould you like to enter the sequence manually?"
                        ):
                            sequence = self._handle_manual_entry(atoms_dict)
                            if sequence:
                                break
                        else:
                            print("\nResetting retry count for one final attempt...")
                            retry_count = 0
                            tried_thresholds.clear()
                            continue
                    
                    # Get new threshold with guidance
                    current_cutoff = self._get_new_threshold(current_cutoff, tried_thresholds)
                    tried_thresholds.add(current_cutoff)
                    processor.distance_cutoff = current_cutoff
                    continue
                    
                else:  # Manual entry
                    sequence = self._handle_manual_entry(atoms_dict)
                    if sequence:
                        break
        
        # If we've exceeded max retries without finding a sequence, fall back to manual entry
        if retry_count >= max_retries and not sequence:
            print("\nAutomatic detection unsuccessful after maximum attempts.")
            sequence = self._handle_manual_entry(atoms_dict)
        
        # Analyze final sequence if it exists
        if sequence and len(sequence) > 1:
            self._analyze_sequence_structure(sequence, pdb_file)
        
        return sequence

    def _get_new_threshold(self, current_cutoff: float, tried_thresholds: set) -> float:
        """Helper method to get a new threshold value with guidance."""
        while True:
            try:
                suggestion = ""
                if tried_thresholds:
                    min_tried = min(tried_thresholds)
                    max_tried = max(tried_thresholds)
                    if not any(t for t in tried_thresholds if min_tried < t < max_tried):
                        suggestion = f"\nSuggestion: Try a value between {min_tried}Å and {max_tried}Å"
                
                print(f"\nPreviously tried thresholds: {', '.join(f'{t}Å' for t in sorted(tried_thresholds))}")
                if suggestion:
                    print(suggestion)
                
                # Clear cached input before threshold prompt
                self.interaction_manager.input_dict.pop('new_threshold', None)
                
                new_cutoff = self.interaction_manager.prompt(
                    "new_threshold",
                    f"Enter new distance threshold (current: {current_cutoff}Å):",
                    input_type=float
                )
                
                if new_cutoff <= 0:
                    print("Threshold must be positive")
                    continue
                    
                if new_cutoff in tried_thresholds:
                    print(f"\nWarning: {new_cutoff}Å has already been tried. Consider a different value.")
                    continue
                
                return new_cutoff
                
            except ValueError:
                print("Please enter a valid number")

    def _handle_manual_entry(self, atoms_dict: Dict) -> Optional[List[int]]:
        """Helper method to handle manual sequence entry."""
        while True:
            # Clear cached inputs before manual entry prompts
            self.interaction_manager.input_dict.pop('manual_sequence', None)
            self.interaction_manager.input_dict.pop('confirm_manual_sequence', None)
            
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
            
            # Validate all hemes exist in structure
            if not all(heme_id in atoms_dict for heme_id in sequence):
                print("Error: Some heme IDs not found in the PDB structure.")
                continue
            
            # Confirm sequence
            if self.interaction_manager.yes_no_prompt(
                "confirm_manual_sequence",
                f"Confirm sequence: {' → '.join(map(str, sequence))}?"
            ):
                return sequence
                
            return None

    def run(self, pdb_file: str) -> EnergeticParameters:
        """
        Execute the energetic evaluation workflow with flexible calculator selection.
        Handles both multi-heme and single-heme cases.
        """
        # Select heme sequence
        sequence = self.select_heme_sequence(pdb_file)
        
        # Initialize parameters container
        computed_params = EnergeticParameters()
        dielectric_constants = None
        
        # Check if we have a single heme or multiple hemes
        is_single_heme = len(sequence) == 1
        
        if is_single_heme:
            print(f"\nWorking with a single heme (ID: {sequence[0]})")
            print("Only reaction free energy calculation is applicable for a single heme.")

        while True:
            # Only clear cached selections if they're not from the input file
            # This modification preserves values from input.txt
            input_file_keys = set()
            if hasattr(self.interaction_manager, 'input_file') and self.interaction_manager.input_file.exists():
                try:
                    with open(self.interaction_manager.input_file, 'r') as f:
                        for line in f:
                            if '=' in line:
                                key = line.split('=')[0].strip()
                                input_file_keys.add(key)
                except Exception:
                    pass
                    
            # Only remove if not from input file
            if 'calculator_selection' not in input_file_keys:
                self.interaction_manager.input_dict.pop('calculator_selection', None)
            if 'use_existing_interactions' not in input_file_keys:
                self.interaction_manager.input_dict.pop('use_existing_interactions', None)

            # Define calculation options based on whether we have single or multiple hemes
            if is_single_heme:
                calc_options = [
                    "Reaction Free Energy",
                    "Exit Program"
                ]
            else:
                calc_options = [
                    "Reorganization Energy",
                    "Reaction Free Energy",
                    "Interaction Energies",
                    "Cooperativity Analysis",
                    "Electronic Couplings",
                    "Marcus Theory Rates",
                    "Exit Program"
                ]
            
            # Present numbered list of options
            print("\nSelect calculators to run (enter numbers separated by spaces):")
            for i, option in enumerate(calc_options, 1):
                print(f"{i}) {option}")
                
            if not is_single_heme:
                print("0) Compute All Quantities")
            
            # Get user selection one at a time
            valid_choices = [str(i) for i in range(1, len(calc_options) + 1)]
            if not is_single_heme:
                valid_choices.append("0")
                
            selection = self.interaction_manager.prompt(
                "calculator_selection",
                "Enter calculator number:",
                choices=valid_choices
            )

            # Convert to int
            selected_num = int(selection)

            # For multi-heme: Check for exit selection (last option)
            if selected_num == len(calc_options):  # If last option (Exit) is selected
                print("\nExiting program...")
                return computed_params

            # Handle "Compute All" selection (only for multi-heme case)
            if not is_single_heme and selected_num == 0:
                selected_nums = list(range(1, len(calc_options)))
            else:
                selected_nums = [selected_num]

            try:
                # Single heme case: Map the selected number to the appropriate calculator
                if is_single_heme:
                    # Map option 1 to reaction free energy for single heme
                    if selected_num == 1:  # Reaction Free Energy
                        print("\nCalculating reaction free energy...")
                        # For single heme, we need to handle dielectric constants differently
                        # since we can't compute them from reorganization energy
                        if not computed_params.dielectric_constants:
                            computed_params.dielectric_constants = self._get_manual_dielectric_constants()
                        
                        computed_params.delta_g_values = self.compute_delta_g(
                            sequence=sequence,
                            dielectric_constants=computed_params.dielectric_constants
                        )
                        print("Reaction free energy calculation completed.")
                
                # Multi-heme case: Original logic
                else:
                    # Compute reorganization energy
                    if 1 in selected_nums:
                        print("\nCalculating reorganization energy...")
                        lambda_vals, diel_consts = self.compute_lambda(sequence, pdb_file)
                        computed_params.lambda_values = lambda_vals
                        computed_params.dielectric_constants = diel_consts
                        print("Reorganization energy calculation completed.")

                    # Compute reaction free energy
                    if 2 in selected_nums:
                        if not computed_params.dielectric_constants:
                            print(
                                "Warning: Dielectric constants are needed."
                                "You can enter them manually, or estimate them by computing the reorganization energy.")
                            computed_params.dielectric_constants = self._get_manual_dielectric_constants()
                    
                        print("\nCalculating reaction free energy...")
                        computed_params.delta_g_values = self.compute_delta_g(
                            sequence=sequence,
                            dielectric_constants=computed_params.dielectric_constants
                        )
                        print("Reaction free energy calculation completed.")

                    # Compute interaction energy matrix
                    if 3 in selected_nums:
                        if not computed_params.dielectric_constants:
                            print(
                                "Warning: Dielectric constants are needed. "
                                "You can enter them manually, or estimate them by computing the reorganization energy."
                            )
                            computed_params.dielectric_constants = self._get_manual_dielectric_constants()
                    
                        print("\nCalculating interaction energies...")
                        computed_params.interaction_energies = self.compute_interactions(
                            sequence=sequence, 
                            dielectric_constants=computed_params.dielectric_constants
                        )
                        print("Interaction energies calculation completed.")

                    # Analyze Cooperativities
                    if 4 in selected_nums:
                        print("\nLaunching cooperativity analysis...")
                        cooperativity_results = self.analyze_cooperativity(
                            sequence=sequence,
                            existing_interactions=computed_params.interaction_energies
                        )
                        
                        if cooperativity_results:
                            # Store full results
                            computed_params.cooperativity_results = cooperativity_results
                            
                            # Extract and store DG values for potential rate calculations
                            if 'sequential' in cooperativity_results:
                                computed_params.sequential_dg_values = [
                                    level_results['delta_G']
                                    for level_results in cooperativity_results['sequential']
                                ]
                            
                            if 'geometric' in cooperativity_results:
                                computed_params.geometric_dg_values = cooperativity_results['geometric']['delta_G']
                            
                            print("\nCooperativity analysis complete!")
                            print(f"Results saved in: {self.ee_dir}/cooperativity_analysis/")

                    # Compute electronic couplings
                    if 5 in selected_nums:
                        print("\nCalculating electronic couplings...")
                        computed_params.coupling_values = self.compute_coupling(
                            sequence, pdb_file)
                        print("Electronic coupling calculation completed.")

                    # Compute Marcus Theory Rates
                    if 6 in selected_nums:
                        if (not computed_params.lambda_values or 
                            not computed_params.coupling_values):
                            print("\nWarning: Marcus rates calculation requires:")
                            print("- Reorganization energy (Option 1)")
                            print("- Electronic coupling (Option 5)")
                            print("Please calculate these quantities first.")
                            continue
                        
                        print("\nCalculating Marcus theory rates...")
                        try:
                            computed_params.rates_forward, computed_params.rates_backward = self.compute_rates(
                                sequence=sequence,
                                lambda_values=computed_params.lambda_values,
                                delta_g_values=[],  # Placeholder - will be selected during computation
                                coupling_values=computed_params.coupling_values
                            )
                            print("Marcus rates calculation completed.")
                        except Exception as e:
                            print(f"\nError during rate calculation: {str(e)}")
                            continue

            except Exception as e:
                print(f"\nError during calculation: {str(e)}")
                print("Returning to calculator selection menu...")
                continue

            # After successful calculations, ask if user wants to perform more
            more_calcs = self.interaction_manager.yes_no_prompt(
                "continue_calculations",
                "\nWould you like to perform additional calculations?"
            )
            
            if not more_calcs:
                return computed_params

    def _get_manual_dielectric_constants(self) -> List[float]:
        """
        Get dielectric constants manually from user input.
        This is used when we can't calculate them from reorganization energy (single heme case).
        
        Returns:
            List of dielectric constants [protein_dielectric, solvent_dielectric]
        """
        print("\nEntering dielectric constants manually:")
        
        # Get protein dielectric constant
        protein_diel = self.interaction_manager.prompt(
            "protein_dielectric",
            "Enter protein dielectric constant:",
            input_type=float
        )
        
        # Get solvent dielectric constant
        solvent_diel = self.interaction_manager.prompt(
            "solvent_dielectric",
            "Enter solvent dielectric constant:",
            input_type=float
        )
        
        # Store in interaction manager for reuse in PBSA calculations
        self.interaction_manager.input_dict["epsout"] = str(solvent_diel)
        
        return [protein_diel, solvent_diel]

    def compute_lambda(self, sequence: List[int], pdb_file: str) -> Tuple[List[float], List[float]]:
        """
        Compute reorganization energies.
        
        Returns:
            Tuple of (lambda_values, dielectric_constants)
        """
        lambda_values, dielectric_constants = self.lambda_calculator.compute_reorganization_energy(
            sequence=sequence,
            pdb_file=pdb_file
        )
        return lambda_values, dielectric_constants

    def compute_delta_g(
        self,
        sequence: List[int],
        dielectric_constants: Optional[List[float]] = None
    ) -> List[float]:
        """
        Compute reaction free energies.

        Args:
            sequence: List of heme IDs
            dielectric_constants: Optional list of dielectric constants from lambda calculation

        Returns:
            List of delta G values
        """
        return self.delta_g_calculator.compute_reaction_free_energy(
            sequence=sequence,
            dielectric_constants=dielectric_constants
        )

    def compute_interactions(
        self, 
        sequence: List[int], 
        dielectric_constants: List[float]
    ) -> Dict[Tuple[int, int], float]:
        """
        Compute heme-heme interaction energies.
        
        Returns:
            Energy matrix as a dictionary mapping (heme_i, heme_j) to interaction energy
        """
        return self.heme_interaction_calculator.compute_heme_interactions(
            sequence=sequence,
            dielectric_constants=dielectric_constants
        )

    def analyze_cooperativity(self,
                            sequence: List[int]
        ) -> Optional[Dict]:
        """
        Launch cooperativity analysis for heme interactions.
        
        Args:
            sequence: List of heme IDs
            
        Returns:
            Dictionary containing analysis results or None if analysis fails
        """
        try:
            print("\nLaunching cooperativity analysis...")
            results = self.cooperativity_analyzer.analyze_cooperativity(sequence=sequence)
            
            if results:
                print("\nCooperativity analysis complete!")
                print(f"Results saved in: {self.ee_dir}/cooperativity_analysis/")
                return results
            
            return None
            
        except Exception as e:
            print(f"\nError during analysis: {str(e)}")
            return None

    def _show_analysis_summary(self, results: Dict) -> None:
        """Show brief summary of analysis results."""
        from rich.console import Console
        from rich.table import Table
        
        console = Console()
        
        # Create summary table
        table = Table(title="Cooperativity Analysis Summary")
        
        if 'sequential' in results:
            table.add_row("Sequential Model", "")
            table.add_row("Oxidation Order", 
                        " → ".join(map(str, results['sequential']['oxidation_order'])))
            table.add_row("Total ΔG", 
                        f"{sum(dg[2] for dg in results['sequential']['delta_G']):.3f} eV")
                        
        if 'geometric' in results:
            table.add_row("Geometric Model", "")
            table.add_row("Oxidation Order",
                        " → ".join(map(str, results['geometric']['oxidation_order'])))
            table.add_row("Total ΔG",
                        f"{sum(dg[2] for dg in results['geometric']['delta_G']):.3f} eV")
        
        console.print(table)

    def compute_coupling(self, sequence: List[int], pdb_file: str) -> List[float]:
        """
        Compute electronic coupling values.
        
        Returns:
            List of coupling values in meV
        """
        return self.coupling_calculator.compute_electronic_coupling(
            sequence=sequence,
            pdb_file=pdb_file
        )

    def _preview_rates(self,
                    available_dgs: List[List[float]],
                    descriptions: List[str],
                    sequence: List[int],
                    lambda_values: List[float],
                    coupling_values: List[float]) -> None:
        """
        Calculate and display preview of rates for all available DG sets.
        
        Args:
            available_dgs: List of available DG value sets
            descriptions: Descriptions for each DG set
            sequence: List of heme IDs
            lambda_values: List of reorganization energies
            coupling_values: List of coupling values
        """
        print("\nPreviewing rates for each DG set...")
        print("=" * 50)

        # Store results for each set
        preview_results = []
        
        for dg_set, desc in zip(available_dgs, descriptions):
            try:
                # Calculate rates using this DG set without writing report
                forward_rates, reverse_rates = self.rate_calculator.compute_marcus_rates(
                    lambda_values=lambda_values,
                    delta_g_values=dg_set,
                    coupling_values=coupling_values,
                    write_report=False  # Don't write report for previews
                )
                
                preview_results.append({
                    'description': desc,
                    'dg_values': dg_set,
                    'forward_rates': forward_rates,
                    'reverse_rates': reverse_rates
                })
                
            except Exception as e:
                print(f"\nError calculating rates for {desc}: {str(e)}")
                continue
        
        # Display comparison table
        if preview_results:
            from rich.console import Console
            from rich.table import Table
            
            console = Console()
            table = Table(title="Rate Preview Comparison")
            
            # Add columns
            table.add_column("Step", style="cyan")
            table.add_column("DG Source", style="green")
            table.add_column("ΔG (eV)", justify="right")
            table.add_column("Forward Rate (s⁻¹)", justify="right")
            table.add_column("Reverse Rate (s⁻¹)", justify="right")
            
            # Add data rows
            for result in preview_results:
                desc = result['description']
                for step, (dg, kf, kr) in enumerate(zip(
                    result['dg_values'],
                    result['forward_rates'],
                    result['reverse_rates']
                )):
                    heme1, heme2 = sequence[step], sequence[step + 1]
                    table.add_row(
                        f"{heme1}→{heme2}",
                        desc,
                        f"{dg:.3f}",
                        f"{kf:.2E}",
                        f"{kr:.2E}"
                    )
                # Add blank row between DG sets
                table.add_row("", "", "", "", "")
            
            # Print table
            console.print(table)
            
            # Add summary analysis
            print("\nSummary Analysis:")
            print("-" * 20)
            for result in preview_results:
                desc = result['description']
                kf_avg = sum(result['forward_rates']) / len(result['forward_rates'])
                kr_avg = sum(result['reverse_rates']) / len(result['reverse_rates'])
                print(f"\n{desc}:")
                print(f"  Average forward rate: {kf_avg:.2E}")
                print(f"  Average reverse rate: {kr_avg:.2E}")
                print(f"  Forward/Reverse ratio: {kf_avg/kr_avg:.2f}")
        else:
            print("\nNo valid rate calculations to preview.")

    def _select_dg_values(self,
                            sequence: List[int],
                            lambda_values: List[float],
                            coupling_values: List[float]) -> Optional[List[float]]:
        """Present interface for selecting DG values to use in rate calculations."""
        available_dgs = []
        descriptions = []

        # Collect all available DG sources
        if self.computed_params.delta_g_values:
            available_dgs.append(self.computed_params.delta_g_values)
            descriptions.append("Standard DG calculator")

        # Add cooperativity DGs if available
        if hasattr(self.computed_params, 'dg_values_dict'):
            for name, dg_values in self.computed_params.dg_values_dict.items():
                if dg_values:  # Only add if we have values
                    available_dgs.append(dg_values)
                    descriptions.append(f"Cooperativity analysis: {name}")

        def handle_manual_entry():
            print("\nEnter ΔG values to be used for rate calculations:")
            num_steps = len(sequence) - 1
            manual_dgs = []

            for i in range(num_steps):
                while True:
                    dg = self.interaction_manager.prompt(
                        f"manual_dg_{i}",
                        f"Step {i+1} ({sequence[i]} → {sequence[i+1]}) ΔG (eV): ",
                        input_type=float
                    )

                    # Confirm value
                    confirm = self.interaction_manager.yes_no_prompt(
                        f"confirm_dg_{i}",
                        f"Confirm ΔG = {dg:.3f} eV?"
                    )
                    if confirm:
                        manual_dgs.append(dg)
                        break

            return manual_dgs

        while True:
            if not available_dgs:
                print("\nNo DG values available from previous calculations.")
                return handle_manual_entry()
            
            # Show all available options
            print("\nAvailable sources for ΔG values:")
            valid_choices = []

            # Show numbered options for DG sets
            for i, desc in enumerate(descriptions, 1):
                dg_set = available_dgs[i-1]
                print(f"\n{i}) {desc}")
                print(f"   Number of ΔGs: {len(dg_set)}")
                print(f"   Values (eV): {', '.join(f'{dg:.3f}' for dg in dg_set)}")
                valid_choices.append(str(i))

            # Add preview and manual options
            print("\nP) Preview rates for all DG sets")
            print("M) Enter DG values manually (for computing rates)")
            valid_choices.extend(['P', 'M'])

            # Get user selection
            choice = self.interaction_manager.prompt(
                "dg_source_selection",
                "\nSelect option:",
                choices=valid_choices,
                input_type=str
            ).upper()

            if choice == 'P':
                # Show rates for all sets
                self._preview_rates(
                    available_dgs,
                    descriptions,
                    sequence,
                    lambda_values,
                    coupling_values
                )
                
                # After preview, prompt for selection
                print("\nAfter preview, please select a DG set to use:")
                select_choice = self.interaction_manager.prompt(
                    "post_preview_selection",
                    "Enter number (1-{}) or 'M' for manual entry: ".format(len(descriptions)),
                    choices=[str(i) for i in range(1, len(descriptions) + 1)] + ['M'],
                    input_type=str
                ).upper()
                
                if select_choice == 'M':
                    return handle_manual_entry()
                else:
                    return available_dgs[int(select_choice) - 1]

            elif choice == 'M':
                return handle_manual_entry()

            else:
                # Return selected DG set directly
                return available_dgs[int(choice) - 1]

    def _display_final_rates(self,
                            sequence: List[int],
                            dg_values: List[float],
                            lambda_values: List[float],
                            coupling_values: List[float],
                            forward_rates: List[float],
                            reverse_rates: List[float],
                            source_description: str = "Selected DG set") -> None:
        """
        Display a formatted table of the final calculated rates and parameters 
        and generate plots.

        Args:
            sequence: List of heme IDs
            dg_values: List of delta G values used
            lambda_values: List of reorganization energies
            coupling_values: List of electronic coupling values
            forward_rates: List of calculated forward rates
            reverse_rates: List of calculated reverse rates
            source_description: Description of the DG source used
        """
        from rich.console import Console
        from rich.table import Table

        console = Console()
        table = Table(title="Final Calculated Rates and Parameters")

        # Add columns
        table.add_column("Step", style="cyan")
        table.add_column("Hda (meV)", justify="right")
        table.add_column("λ (eV)", justify="right")
        table.add_column("ΔG (eV)", justify="right")
        table.add_column("Forward Rate (s⁻¹)", justify="right")
        table.add_column("Reverse Rate (s⁻¹)", justify="right")

        # Add data rows
        for i, (dg, lam, hda, kf, kr) in enumerate(zip(
            dg_values, lambda_values, coupling_values, forward_rates, reverse_rates)):
            heme1, heme2 = sequence[i], sequence[i + 1]

            table.add_row(
                f"{heme1}→{heme2}",
                f"{hda:.3f}",
                f"{lam:.3f}",
                f"{dg:.3f}",
                f"{kf:.2E}",
                f"{kr:.2E}"
            )

        # Print table with source info
        print(f"\nSource of ΔG values: {source_description}")
        console.print(table)

        # Generate plots
        print("\nGenerating analysis plots...")
        try:
            # Store the values in the rate calculator for plotting
            self.rate_calculator.lambda_values = lambda_values
            self.rate_calculator.delta_g_values = dg_values
            self.rate_calculator.coupling_values = coupling_values
            self.rate_calculator.forward_rates = forward_rates
            self.rate_calculator.reverse_rates = reverse_rates

            # Generate plots
            plot_path = self.ee_dir / "rate_analysis.png"
            self.rate_calculator.plot_analysis(
                pdb_file=self.pdb_file,  # We have this from initialization
                save_path=str(plot_path),
                dpi=300
            )
            print(f"Analysis plots saved to: {plot_path}")

        except Exception as e:
            print(f"Warning: Could not generate analysis plots: {str(e)}")


    def compute_rates(self,
                sequence: List[int],
                lambda_values: List[float],
                delta_g_values: List[float],
                coupling_values: List[float]) -> Tuple[List[float], List[float]]:
        """
        Compute electron transfer rates.

        Now includes generation of comprehensive analysis plots.
        """
        # Get user-selected DG values and source description
        selected_dgs = self._select_dg_values(
            sequence=sequence,
            lambda_values=lambda_values,
            coupling_values=coupling_values
        )

        if selected_dgs is None:
            raise ValueError("No DG values selected for rate calculation.")

        # Validate lengths match
        if not (len(selected_dgs) == len(lambda_values) == len(coupling_values)):
            raise ValueError(
                f"Mismatched parameter lengths: "
                f"DGs={len(selected_dgs)}, "
                f"lambdas={len(lambda_values)}, "
                f"couplings={len(coupling_values)}"
            )

        # Determine the source description based on the selection
        if hasattr(self.computed_params, 'dg_values_dict'):
            for name, dg_values in self.computed_params.dg_values_dict.items():
                if dg_values == selected_dgs:
                    source_description = f"Cooperativity analysis: {name}"
                    break
            else:
                if self.computed_params.delta_g_values == selected_dgs:
                    source_description = "Standard DG calculator"
                else:
                    source_description = "Manually entered values"
        else:
            source_description = "Manually entered values"

        # Compute final rates with selected DGs
        print("\nCalculating final rates with selected DG values...")
        forward_rates, reverse_rates = self.rate_calculator.compute_marcus_rates(
            lambda_values=lambda_values,
            delta_g_values=selected_dgs,
            coupling_values=coupling_values
        )

        # Store results in computed_params for potential reuse
        self.computed_params.lambda_values = lambda_values
        self.computed_params.delta_g_values = selected_dgs
        self.computed_params.coupling_values = coupling_values
        self.computed_params.rates_forward = forward_rates
        self.computed_params.rates_backward = reverse_rates

        # Display final rates table and generate plots
        self._display_final_rates(
            sequence=sequence,
            dg_values=selected_dgs,
            lambda_values=lambda_values,
            coupling_values=coupling_values,
            forward_rates=forward_rates,
            reverse_rates=reverse_rates,
            source_description=source_description
        )

        return forward_rates, reverse_rates
