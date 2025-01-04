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
            
            # Calculate plane angle
            plane_angle = processor.calculate_plane_angle(
                atoms_dict[heme1], 
                atoms_dict[heme2]
            )
            
            # Classify stacking
            stacking_type = processor.classify_stacking(plane_angle)
            
            # Add row to table
            table.add_row(
                f"{heme1} → {heme2}", 
                f"{min_distance:.2f}", 
                f"{plane_angle:.2f}", 
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

        # Initial distance cutoff
        current_cutoff = 20.0
        processor.distance_cutoff = current_cutoff
      
        print("\nAutomatically detecting possible heme sequences. Please wait! ...")
        sys.stdout.flush()

        start_detect = time.time()
        while True:
            try:
                # Try linear topology first
                print("Detecting linear sequence...")
                sys.stdout.flush()
                linear_sequence = processor.detect_linear_sequence(atoms_dict)
                
                # Try branched topology
                print("Detecting branched sequences...")
                sys.stdout.flush()
                branched_sequences = []
                try:
                    branched_sequences = processor.detect_branched_sequence(atoms_dict)
                except ValueError:
                    pass
               
                print(f"Sequence detection complete. Total time: {time.time() - start_time:.2f} seconds")
                sys.stdout.flush()

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
                    processor.distance_cutoff = current_cutoff
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
                        min_distance = processor.calculate_min_distance(
                            atoms_dict[heme1], 
                            atoms_dict[heme2]
                        )
                        
                        # Calculate plane angle
                        plane_angle = processor.calculate_plane_angle(
                            atoms_dict[heme1], 
                            atoms_dict[heme2]
                        )
                        
                        # Classify stacking
                        stacking_type = processor.classify_stacking(plane_angle)
                        
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

    def run(self, pdb_file: str) -> EnergeticParameters:
        """
        Execute the energetic evaluation workflow with flexible calculator selection.
        """
        # Select heme sequence
        sequence = self.select_heme_sequence(pdb_file)
    
        # Use class instance's computed_params
        self.computed_params = EnergeticParameters()
        computed_params = self.computed_params  # Local reference for convenience
        dielectric_constants = None

        while True:
            # Clear cached calculator selection before showing menu
            self.interaction_manager.input_dict.pop('calculator_selection', None)
            self.interaction_manager.input_dict.pop('use_existing_interactions', None)

            # Define calculation options
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
            print("0) Compute All Quantities")
            
            # Get user selection one at a time
            valid_choices = [str(i) for i in range(len(calc_options) + 1)]
            selection = self.interaction_manager.prompt(
                "calculator_selection",
                "Enter calculator number:",
                choices=valid_choices
            )

            # Convert to int
            selected_num = int(selection)

            # Check for exit selection
            if selected_num == len(calc_options):  # If last option (Exit) is selected
                print("\nExiting program...")
                return computed_params

            # Handle "Compute All" selection
            if selected_num == 0:
                selected_nums = list(range(1, len(calc_options)))
            else:
                selected_nums = [selected_num]

            try:
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
                
                    print("\nCalculating interaction energies...")
                    computed_params.interaction_energies = self.compute_interactions(
                        sequence=sequence, 
                        dielectric_constants=computed_params.dielectric_constants  # Now using correct dielectric constants
                    )
                    print("Interaction energies calculation completed.")

                # Analyze Cooperativities
                if 4 in selected_nums:
                    print("\nLaunching cooperativity analysis...")
                    cooperativity_results = self.analyze_cooperativity(sequence=sequence)
                    
                    if cooperativity_results:
                        # Store full results in the instance's computed_params
                        self.computed_params.cooperativity_results = cooperativity_results
                        computed_params.cooperativity_results = cooperativity_results
                        
                        # Extract and store DG values for potential rate calculations
                        if 'sequential' in cooperativity_results:
                            seq_results = cooperativity_results['sequential']
                            if isinstance(seq_results, (list, tuple)):
                                # Multiple interaction levels
                                for i, level_result in enumerate(seq_results, 1):
                                    if isinstance(level_result, dict) and 'delta_G' in level_result:
                                        dg_values = [dg[2] for dg in level_result['delta_G']]
                                        self.computed_params.dg_values_dict[f'Sequential (Level {i})'] = dg_values
                                        computed_params.dg_values_dict[f'Sequential (Level {i})'] = dg_values
                            elif isinstance(seq_results, dict) and 'delta_G' in seq_results:
                                dg_values = [dg[2] for dg in seq_results['delta_G']]
                                self.computed_params.dg_values_dict['Sequential'] = dg_values
                                computed_params.dg_values_dict['Sequential'] = dg_values
                        
                        if 'geometric' in cooperativity_results:
                            geo_results = cooperativity_results['geometric']
                            if isinstance(geo_results, dict) and 'delta_G' in geo_results:
                                dg_values = [dg[2] for dg in geo_results['delta_G']]
                                self.computed_params.dg_values_dict['Geometric'] = dg_values
                                computed_params.dg_values_dict['Geometric'] = dg_values
                        
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
                # Calculate rates using this DG set
                forward_rates, reverse_rates = self.rate_calculator.compute_marcus_rates(
                    lambda_values=lambda_values,
                    delta_g_values=dg_set,
                    coupling_values=coupling_values
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
            from tabulate import tabulate
            
            # Create table data
            headers = ["Step", "DG Source", "ΔG (eV)", "Forward Rate", "Reverse Rate"]
            table_data = []
            
            for result in preview_results:
                desc = result['description']
                for step, (dg, kf, kr) in enumerate(zip(
                    result['dg_values'],
                    result['forward_rates'],
                    result['reverse_rates']
                )):
                    heme1, heme2 = sequence[step], sequence[step + 1]
                    table_data.append([
                        f"{heme1}→{heme2}",
                        desc,
                        f"{dg:.3f}",
                        f"{kf:.2E}",
                        f"{kr:.2E}"
                    ])
                # Add blank row between DG sets
                table_data.append(["", "", "", "", ""])
            
            print("\nRate Comparison:")
            print(tabulate(table_data, headers=headers, tablefmt="grid"))
            
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
        
        print(self.computed_params)

        # 1. Check standard DG calculator results
        if hasattr(self.computed_params, 'delta_g_values') and self.computed_params.delta_g_values:
            available_dgs.append(self.computed_params.delta_g_values)
            descriptions.append("Standard DG calculator")

        # 2. Check DG values from various sources (including cooperativity)
        if hasattr(self.computed_params, 'dg_values_dict'):
            for name, dg_values in self.computed_params.dg_values_dict.items():
                available_dgs.append(dg_values)
                descriptions.append(f"Cooperativity analysis: {name}")

        # Show available options
        while True:
            print("\nAvailable sources for ΔG values:")
            # Show numbered options for DG sets
            for i, desc in enumerate(descriptions, 1):
                dg_set = available_dgs[i-1]
                print(f"\n{i}) {desc}")
                print(f"   Number of ΔGs: {len(dg_set)}")
                print(f"   Values (eV): {', '.join(f'{dg:.3f}' for dg in dg_set)}")

            # Show additional options
            print("\nP) Preview rates for all DG sets")
            print("M) Enter DG values manually")
            
            # Get user selection
            choice = self.interaction_manager.prompt(
                "dg_source_selection",
                "\nSelect option (number, P, or M):",
                choices=[str(i) for i in range(1, len(descriptions) + 1)] + ['P', 'M']
            ).upper()
            
            if choice == 'P':
                # Show rate preview
                self._preview_rates(
                    available_dgs,
                    descriptions,
                    sequence,
                    lambda_values,
                    coupling_values
                )
                # Return to selection menu
                continue
                
            elif choice == 'M':
                print("\nEntering ΔG values manually...")
                num_steps = len(sequence) - 1
                manual_dgs = []
                
                print("\nEnter ΔG value for each electron transfer step (eV):")
                for i in range(num_steps):
                    while True:
                        dg = self.interaction_manager.prompt(
                            f"manual_dg_{i}",
                            f"Step {i+1} ({sequence[i]} → {sequence[i+1]}): ",
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
            
            else:
                # Return selected DG set
                choice_idx = int(choice) - 1
                return available_dgs[choice_idx]

    def compute_rates(self,
                    sequence: List[int],
                    lambda_values: List[float],
                    delta_g_values: List[float],
                    coupling_values: List[float]) -> Tuple[List[float], List[float]]:
        """
        Compute electron transfer rates.
        
        Args:
            sequence: List of heme IDs
            lambda_values: Reorganization energies
            delta_g_values: IGNORED - Will select DGs during execution
            coupling_values: Electronic coupling values
        
        Returns:
            Tuple of (forward_rates, reverse_rates)
        """
        # Get user-selected DG values
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
        
        # Compute final rates with selected DGs
        print("\nCalculating final rates with selected DG values...")
        return self.rate_calculator.compute_marcus_rates(
            lambda_values=lambda_values,
            delta_g_values=selected_dgs,
            coupling_values=coupling_values
        )



