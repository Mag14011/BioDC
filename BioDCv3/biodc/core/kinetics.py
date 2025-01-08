"""
Kinetics Evaluation module for BioDC.
Coordinates access to various kinetic calculation models.
"""
from pathlib import Path
from typing import List, Optional, Dict, Any
from dataclasses import dataclass, field
from rich.console import Console

from biodc.utils.interaction import InteractionManager
from biodc.core.kineval_modules.diffusion_calculator import DiffusionCalculator
from biodc.core.kineval_modules.flux_calculator import FluxCalculator
from biodc.core.kineval_modules.parameter_explorer import MonteCarloOptimizer

@dataclass
class KineticParameters:
    """Container for kinetic calculation parameters."""
    diffusion_results: Optional[Dict] = None
    flux_results: Optional[Dict] = None
    monte_carlo_results: Optional[Dict] = None
    
class KineticEvaluation:
    """
    Coordinates access to different kinetic calculation models.
    
    Provides menu-driven access to:
    - Analytical charge diffusion constants
    - Steady-state protein-limited electron flux
    - Monte Carlo parameter exploration
    """
    def __init__(self, 
                interaction_manager: InteractionManager,
                launch_dir: Path):
        """
        Initialize the kinetic evaluation workflow.
        
        Args:
            interaction_manager: Manages user interactions and input tracking
            launch_dir: Project launch directory
        """
        self.interaction_manager = interaction_manager
        self.launch_dir = launch_dir
        self.console = Console()
        
        # Set up kinetics evaluation directory
        self.ke_dir = launch_dir / "KE"
        self.ke_dir.mkdir(exist_ok=True)
        
        # Reference to energetics evaluation directory
        self.ee_dir = launch_dir / "EE"
        if not self.ee_dir.exists():
            raise FileNotFoundError(
                "EE directory not found. Please run energetics evaluation first."
            )
        
        # Validate environment and set calculator availability
        self._validate_environment()
            
        # Initialize calculators if available
        if self.all_calcs_available:
            self.diffusion_calculator = DiffusionCalculator(
                interaction_manager=interaction_manager,
                launch_dir=launch_dir,
                ee_dir=self.ee_dir
            )
            
            self.flux_calculator = FluxCalculator(
                interaction_manager=interaction_manager,
                launch_dir=launch_dir,
                ee_dir=self.ee_dir
            )
        
        # Monte Carlo optimizer always available if min.pdb exists
        self.monte_carlo_optimizer = MonteCarloOptimizer(
            interaction_manager=interaction_manager,
            launch_dir=launch_dir,
            ee_dir=self.ee_dir
        )

    def _validate_environment(self) -> None:
        """
        Validate available files and set calculator availability.
        Raises FileNotFoundError if min.pdb is missing.
        """
        # Check min.pdb (required for all calculations)
        if not (self.ee_dir / "min.pdb").exists():
            raise FileNotFoundError("Required file min.pdb not found in EE directory")
            
        # Check file combinations
        has_rates = (self.ee_dir / "rates.txt").exists()
        has_state_energies = (self.ee_dir / "StateEnergies.txt").exists()
        has_dg = (self.ee_dir / "DG.txt").exists()
        has_lambda = (self.ee_dir / "lambda.txt").exists()
        has_hda = (self.ee_dir / "Hda.txt").exists()
        
        # Set availability flag for diffusion and flux calculators
        self.all_calcs_available = has_rates or (has_lambda and has_hda and (has_dg or has_state_energies))
        
        # Print appropriate warnings
        if not self.all_calcs_available:
            self.console.print("\n[yellow]Warning: Only the Monte Carlo calculator will be available.")
            self.console.print("For full functionality, provide either:")
            self.console.print("  - rates.txt")
            self.console.print("  - or the combination of:")
            self.console.print("    * lambda.txt")
            self.console.print("    * Hda.txt")
            self.console.print("    * and either DG.txt or StateEnergies.txt")
        else:
            # Check for any missing files and warn specifically
            if not has_rates and (not all([has_lambda, has_hda]) or not (has_dg or has_state_energies)):
                self.console.print("\n[yellow]Warning: Some energetic files are missing.")
                if not has_lambda:
                    self.console.print("Missing: lambda.txt")
                if not has_hda:
                    self.console.print("Missing: Hda.txt")
                if not has_dg and not has_state_energies:
                    self.console.print("Missing: DG.txt or StateEnergies.txt")
                self.console.print("Some calculations may fail.")

    def run(self) -> KineticParameters:
        """
        Execute the kinetic evaluation workflow with calculator selection menu.
        
        Returns:
            KineticParameters containing results from selected calculations
        """
        # Initialize parameters container
        computed_params = KineticParameters()

        while True:
            # Clear cached calculator selection before showing menu
            self.interaction_manager.input_dict.pop('calculator_selection', None)

            # Define calculation options based on availability
            if self.all_calcs_available:
                calc_options = [
                    "Analytical Charge Diffusion Constant",
                    "Steady-state Electron Flux",
                    "Monte Carlo Parameter Exploration",
                    "Exit Program"
                ]
                print("\nSelect calculator to run:")
                for i, option in enumerate(calc_options, 1):
                    print(f"{i}) {option}")
                print("0) Run Diffusion and Flux Analysis")
                valid_choices = [str(i) for i in range(len(calc_options) + 1)]
            else:
                # Only show Monte Carlo option
                calc_options = [
                    "Monte Carlo Parameter Exploration",
                    "Exit Program"
                ]
                print("\nSelect calculator to run:")
                for i, option in enumerate(calc_options, 1):
                    print(f"{i}) {option}")
                valid_choices = [str(i) for i in range(1, len(calc_options) + 1)]
            
            # Get user selection
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

            try:
                if self.all_calcs_available:
                    # Handle combined diffusion and flux analysis
                    if selected_num == 0:
                        selected_nums = [1, 2]  # Only run diffusion and flux calculators
                    else:
                        selected_nums = [selected_num]

                    # Run diffusion calculation
                    if 1 in selected_nums:
                        print("\nCalculating analytical diffusion constant...")
                        results = self.diffusion_calculator.compute_diffusion_constant()
                        computed_params.diffusion_results = results
                        self.diffusion_calculator.display_results(results)
                        print("Diffusion constant calculation completed.")

                    # Run flux calculation
                    if 2 in selected_nums:
                        print("\nCalculating steady-state electron flux...")
                        results = self.flux_calculator.compute_steady_state_flux()
                        computed_params.flux_results = results
                        self.flux_calculator.save_results(results)
                        print("Flux calculation completed.")

                    # Run Monte Carlo optimization
                    if 3 in selected_nums:
                        print("\nLaunching Monte Carlo parameter exploration...")
                        results = self.monte_carlo_optimizer.explore_parameters()
                        computed_params.monte_carlo_results = results
                        if results:
                            print("\nMonte Carlo exploration complete!")
                            print(f"Results saved in: {self.ke_dir}/monte_carlo/")
                else:
                    # Only Monte Carlo is available
                    if selected_num == 1:
                        print("\nLaunching Monte Carlo parameter exploration...")
                        results = self.monte_carlo_optimizer.explore_parameters()
                        computed_params.monte_carlo_results = results
                        if results:
                            print("\nMonte Carlo exploration complete!")
                            print(f"Results saved in: {self.ke_dir}/monte_carlo/")

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