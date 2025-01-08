"""
Heme Cooperativity Analysis module for BioDC.
Integrates analyze_heme_cooperativity.py into the energetics evaluation workflow.
"""
import os
from pathlib import Path
from typing import Dict, List, Optional, Any
from dataclasses import dataclass
import numpy as np

@dataclass
class CooperativityParameters:
    """Parameters for cooperativity analysis."""
    energy_shift: float = 0.0
    interaction_scale: float = 1.0
    adjacent_only: bool = True
    make_periodic: bool = False
    model: str = "seq"  # "seq" or "geo"

class HemeCooperativityAnalyzer:
    """
    Analyzes heme cooperativity effects in electron transfer chains.
    Integrates with the main BioDC energetics evaluation system.
    """
    def __init__(self, interaction_manager, ee_dir: Path):
        """
        Initialize the cooperativity analyzer.
        
        Args:
            interaction_manager: BioDC's interaction manager for user prompts
            ee_dir: Path to the EE directory containing EnergyMatrix.txt
        """
        self.interaction_manager = interaction_manager
        self.ee_dir = ee_dir
        self.output_dir = ee_dir / "cooperativity_analysis"
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def analyze_cooperativity(self, sequence: List[int]) -> Optional[Dict]:
            """
            Launch cooperativity analysis using energy matrix from EE directory.
            
            Args:
                sequence: List of heme IDs
                
            Returns:
                Dictionary containing analysis results or None if analysis fails
            """
            try:
                # Get analysis parameters from user
                params = self._get_parameters_from_user()

                # Get redox plot preference
                plot_option = self.interaction_manager.prompt(
                    key="cooperativity_plot_option",
                    message="\nSelect redox fraction plots to generate:\n1) Oxidized fraction\n2) Reduced fraction\n3) Both\nEnter choice number",
                    choices=['1', '2', '3'],
                    input_type=str
                )
                plot_type = {
                    '1': 'ox',
                    '2': 'red',
                    '3': 'both'
                }[plot_option]

                # Import analysis functions
                from biodc.core.engeval_modules.analyze_heme_cooperativity import (
                    process_matrix, find_global_potential_range, calculate_fractions, plot_curves
                )

                matrix_file = self.ee_dir / "EnergyMatrix.txt"
                if not matrix_file.exists():
                    print(f"\nError: EnergyMatrix.txt not found in {self.ee_dir}")
                    return None

                print(f"\nUsing energy matrix from: {matrix_file}")
                
                results = {}
                all_potentials = []
                
                if params.model in ['seq', 'both']:
                    print("\nProcessing sequential model...")
                    seq_results = process_matrix(
                        name="sequential",
                        filepath=str(matrix_file),
                        model="seq",
                        energy_shift=params.energy_shift,
                        interaction_scale=params.interaction_scale,
                        adjacent_only=params.adjacent_only,
                        make_periodic=params.make_periodic,
                        output_dir=str(self.output_dir)
                    )
                    if seq_results:
                        results['sequential'] = seq_results
                        all_potentials.append(seq_results)

                if params.model in ['geo', 'both']:
                    print("\nProcessing geometric model...")
                    geo_results = process_matrix(
                        name="geometric",
                        filepath=str(matrix_file),
                        model="geo",
                        energy_shift=params.energy_shift,
                        interaction_scale=params.interaction_scale,
                        adjacent_only=params.adjacent_only,
                        make_periodic=params.make_periodic,
                        output_dir=str(self.output_dir)
                    )
                    if geo_results:
                        results['geometric'] = geo_results
                        all_potentials.append(geo_results)

                # Generate redox fraction plots if we have results
                # Generate redox fraction plots if we have results
                if results and all_potentials:
                    # Find global potential range for redox plots
                    lower_bound, upper_bound = find_global_potential_range(all_potentials)
                    E_range = np.linspace(lower_bound, upper_bound, 1600)

                    # Calculate fractions and create plots
                    plot_data = {'BioDC': {}}  # All our data is from BioDC
                
                    # For each model type, add its results to BioDC source
                    for name, result in results.items():
                        calculated_fractions = calculate_fractions(E_range, result, name)
                        plot_data['BioDC'].update(calculated_fractions)

                    # Create plots in output directory
                    output_plot_file = os.path.join(self.output_dir, f'redox_plot_{params.model}.png')
                    plot_curves(E_range, plot_data, output_plot_file, plot_type)
                    print(f"\nRedox fraction plots saved to: {output_plot_file}")

                return results if results else None

            except Exception as e:
                print(f"\nError during cooperativity analysis: {str(e)}")
                import traceback
                print(traceback.format_exc())
                return None

    def _get_parameters_from_user(self) -> CooperativityParameters:
        """Get analysis parameters from user interaction."""
        params = CooperativityParameters()

        # Get model type with numbered choices
        model_choice = self.interaction_manager.prompt(
            key="cooperativity_model",
            message="\nSelect analysis model:\n1) Sequential (thermodynamic pathway)\n2) Geometric (structural pathway)\n3) Both\nEnter choice number",
            choices=['1', '2', '3'],
            input_type=str
        )
        
        params.model = {
            '1': 'seq',
            '2': 'geo',
            '3': 'both'
        }[model_choice]

        # Get energy shift with default
        params.energy_shift = self.interaction_manager.prompt(
            key="cooperativity_energy_shift",
            message="\nEnter energy shift (eV)",
            input_type=float,
            default=0.0
        )

        # Get interaction scale with default
        params.interaction_scale = self.interaction_manager.prompt(
            key="cooperativity_interaction_scale",
            message="\nEnter interaction scaling factor",
            input_type=float,
            default=1.0
        )

        # Get adjacent only option with default True
        params.adjacent_only = self.interaction_manager.yes_no_prompt(
            key="cooperativity_adjacent_only",
            message="\nCalculate ΔG for adjacent hemes only?",
            default=True
        )

        # Get periodicity option if adjacent_only is True
        if params.adjacent_only:
            params.make_periodic = self.interaction_manager.yes_no_prompt(
                key="cooperativity_make_periodic",
                message="\nInclude step from last heme back to first?",
                default=True
            )

        return params

    def _generate_report(self, results: Dict[str, Any]) -> None:
        """Generate detailed analysis report."""
        report_file = self.output_dir / "cooperativity_report.txt"
        
        with open(report_file, 'w') as f:
            f.write("Heme Cooperativity Analysis Report\n")
            f.write("=================================\n\n")
            
            if 'sequential' in results:
                f.write("\nSequential (Thermodynamic) Model Results\n")
                f.write("-" * 40 + "\n")
                seq_results = results['sequential']
                
                # Write oxidation order
                f.write("\nOxidation Order:\n")
                oxidation_order = seq_results.get('oxidation_order', [])
                f.write(" → ".join(map(str, oxidation_order)) + "\n")
                
                # Write ΔG values
                f.write("\nΔG Values:\n")
                delta_g = seq_results.get('delta_G', [])
                for dg in delta_g:
                    f.write(f"Heme {dg[0]} → {dg[1]}: {dg[2]:.3f} eV\n")
            
            if 'geometric' in results:
                f.write("\nGeometric Model Results\n")
                f.write("-" * 40 + "\n")
                geo_results = results['geometric']
                
                # Write oxidation order
                f.write("\nOxidation Order:\n")
                oxidation_order = geo_results.get('oxidation_order', [])
                f.write(" → ".join(map(str, oxidation_order)) + "\n")
                
                # Write ΔG values
                f.write("\nΔG Values:\n")
                delta_g = geo_results.get('delta_G', [])
                for dg in delta_g:
                    f.write(f"Heme {dg[0]} → {dg[1]}: {dg[2]:.3f} eV\n")
