"""
Marcus Theory Rate Calculator for Electron Transfer in Multi-Heme Systems.

Computes electron transfer rates using Marcus theory for sequential 
heme-to-heme electron transfer steps.
"""

import math
from typing import List, Tuple
from pathlib import Path

from biodc.utils.interaction import InteractionManager

class RateCalculator:
    """
    Calculator for electron transfer rates using Marcus theory.
    
    Computes forward and reverse rates for electron transfer 
    between consecutive heme groups.
    """
    
    def __init__(self, 
                 interaction_manager: InteractionManager,
                 launch_dir: Path):
        """
        Initialize RateCalculator
        
        Args:
            interaction_manager: Manager for user interactions
            launch_dir: Project launch directory
        """
        self.interaction_manager = interaction_manager
        self.launch_dir = Path(launch_dir)
        self.ee_dir = self.launch_dir / "EE"
        
        # Physical constants
        self.constants = {
            'temperature': 300.0,  # Kelvin
            'pi': 3.141592654,
            'kB': 8.6173304E-5,    # eV/K
            'hbar': 6.582119514E-16  # eV·s
        }
    
    def compute_marcus_rates(
        self, 
        lambda_values: List[float], 
        delta_g_values: List[float], 
        coupling_values: List[float]
    ) -> Tuple[List[float], List[float]]:
        """
        Compute electron transfer rates using Marcus theory
        
        Args:
            lambda_values: Reorganization energies for each step (eV)
            delta_g_values: Reaction free energies for each step (eV)
            coupling_values: Electronic coupling values for each step (meV)
        
        Returns:
            Tuple of forward and reverse rate lists
        """
        # Prepare results storage
        num_steps = len(delta_g_values)
        prefactors = [0] * num_steps
        forward_activation = [0] * num_steps
        reverse_activation = [0] * num_steps
        forward_rates = [0] * num_steps
        reverse_rates = [0] * num_steps
        
        # Compute rates for each step
        for idx in range(num_steps):
            # Conversion: coupling from meV to eV
            hda = coupling_values[idx] / 1000.0
            
            # Prefactor calculation
            prefactors[idx] = (
                2 * self.constants['pi'] * (hda**2) / 
                (self.constants['hbar'] * 
                 math.sqrt(4 * self.constants['pi'] * 
                           lambda_values[idx] * 
                           self.constants['kB'] * 
                           self.constants['temperature']))
            )
            
            # Forward activation energy
            forward_activation[idx] = (
                (delta_g_values[idx] + lambda_values[idx])**2 / 
                (4 * lambda_values[idx])
            )
            
            # Forward rate
            forward_rates[idx] = prefactors[idx] * math.exp(
                (-1 * forward_activation[idx]) / 
                (self.constants['kB'] * self.constants['temperature'])
            )
            
            # Reverse activation energy
            reverse_activation[idx] = (
                ((-1 * delta_g_values[idx]) + lambda_values[idx])**2 / 
                (4 * lambda_values[idx])
            )
            
            # Reverse rate
            reverse_rates[idx] = prefactors[idx] * math.exp(
                (-1 * reverse_activation[idx]) / 
                (self.constants['kB'] * self.constants['temperature'])
            )
        
        # Write rates to file and display results
        self._write_rates_report(
            lambda_values,
            delta_g_values,
            coupling_values,
            forward_activation, 
            reverse_activation, 
            forward_rates, 
            reverse_rates
        )
        
        return forward_rates, reverse_rates
    
    def _write_rates_report(
        self, 
        lambda_values: List[float],
        delta_g_values: List[float],
        coupling_values: List[float],
        forward_activation: List[float],
        reverse_activation: List[float],
        forward_rates: List[float],
        reverse_rates: List[float]
    ):
        """
        Write detailed rates report with input parameters
        
        Args:
            lambda_values: Reorganization energies
            delta_g_values: Reaction free energies
            coupling_values: Electronic coupling values
            forward_activation: Forward activation energies
            reverse_activation: Reverse activation energies
            forward_rates: Forward electron transfer rates
            reverse_rates: Reverse electron transfer rates
        """
        with open(self.ee_dir / "rates.txt", 'w') as f:
            f.write("Marcus Theory Electron Transfer Rates Analysis\n")
            f.write("=" * 50 + "\n\n")
            
            # Input Parameters Section
            f.write("Input Parameters:\n")
            f.write("-" * 20 + "\n")
            f.write("Step\tλ (eV)\tΔG (eV)\tHda (meV)\n")
            f.write("-" * 45 + "\n")
            for idx, (lamb, dg, hda) in enumerate(
                zip(lambda_values, delta_g_values, coupling_values)
            ):
                f.write(
                    f"{idx+1}\t"
                    f"{lamb:.3f}\t"
                    f"{dg:.3f}\t"
                    f"{hda:.3f}\n"
                )
            f.write("\n")
            
            # Detailed Rates Analysis Section
            f.write("Electron Transfer Rates Analysis:\n")
            f.write("-" * 30 + "\n")
            f.write("Step\tForward Act.\tReverse Act.\tForward Rate\tReverse Rate\n")
            f.write("-" * 70 + "\n")
            
            for idx, (f_act, r_act, f_rate, r_rate) in enumerate(
                zip(forward_activation, reverse_activation, 
                    forward_rates, reverse_rates)
            ):
                f.write(
                    f"{idx+1}\t"
                    f"{f_act:.3E}\t"
                    f"{r_act:.3E}\t"
                    f"{f_rate:.3E}\t"
                    f"{r_rate:.3E}\n"
                )
        
        # Separate CSV file for easy parsing
        with open(self.ee_dir / "rates.csv", 'w') as f:
            f.write("step,lambda,delta_g,hda,forward_activation,reverse_activation,forward_rate,reverse_rate\n")
            
            for idx, (lamb, dg, hda, f_act, r_act, f_rate, r_rate) in enumerate(
                zip(lambda_values, delta_g_values, coupling_values,
                    forward_activation, reverse_activation, 
                    forward_rates, reverse_rates)
            ):
                f.write(
                    f"{idx+1},"
                    f"{lamb:.3E},"
                    f"{dg:.3E},"
                    f"{hda:.3E},"
                    f"{f_act:.3E},"
                    f"{r_act:.3E},"
                    f"{f_rate:.3E},"
                    f"{r_rate:.3E}\n"
                )
