"""
Marcus Theory Rate Calculator for Electron Transfer in Multi-Heme Systems.

Computes electron transfer rates using Marcus theory for sequential 
heme-to-heme electron transfer steps. Includes plotting functionality
for comprehensive analysis.
"""

import math
from typing import List, Tuple
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import seaborn as sns
from biodc.utils.interaction import InteractionManager
from biodc.utils.structure_analyzer import PDBProcessor

class RateCalculator:
    """
    Calculator for electron transfer rates using Marcus theory.
    
    Computes forward and reverse rates for electron transfer 
    between consecutive heme groups.
    """
    
    class RatePlotter:
        """Inner class for generating plots of electron transfer analysis."""
        
        def __init__(self):
            """Initialize the plotter with style settings."""
            plt.style.use('seaborn-darkgrid')
            self.colors = {
                'forward': '#2E86C1',  # Blue
                'reverse': '#E74C3C',  # Red
                'single': '#2ECC71'    # Green
            }
            
        def create_analysis_figure(
            self,
            lambda_values: List[float],
            delta_g_values: List[float],
            coupling_values: List[float],
            forward_activation: List[float],
            reverse_activation: List[float],
            forward_rates: List[float],
            reverse_rates: List[float],
            pdb_file: str,
            sequence_file: str = None,
            topology: str = 'linear',
            save_path: str = None,
            dpi: int = 300
        ) -> None:
            """
            Create a comprehensive figure with six panels analyzing electron transfer.
            
            Args:
                lambda_values: List of reorganization energies
                delta_g_values: List of reaction free energies
                coupling_values: List of electronic coupling values
                forward_activation: List of forward activation energies
                reverse_activation: List of reverse activation energies
                forward_rates: List of forward rates
                reverse_rates: List of reverse rates
                pdb_file: Path to PDB file for structural analysis
                sequence_file: Optional path to sequence file
                topology: Structure topology ('linear' or 'branched')
                save_path: Optional path to save figure
                dpi: Resolution for saved figure
            """
            # Create figure and GridSpec
            fig = plt.figure(figsize=(14, 9))
            gs = gridspec.GridSpec(2, 3, figure=fig)
            gs.update(wspace=0.3, hspace=0.3)
            
            # Generate sequence indices for x-axis
            x_seq = np.arange(1, len(lambda_values) + 1)
            
            # Panel A: Coupling vs sequence
            ax1 = fig.add_subplot(gs[0, 0])
            self._plot_coupling(ax1, x_seq, coupling_values)
            ax1.text(-0.1, 1.1, '(A)', transform=ax1.transAxes, fontsize=12, fontweight='bold')
            
            # Panel B: Lambda vs sequence
            ax2 = fig.add_subplot(gs[0, 1])
            self._plot_lambda(ax2, x_seq, lambda_values)
            ax2.text(-0.1, 1.1, '(B)', transform=ax2.transAxes, fontsize=12, fontweight='bold')
            
            # Panel C: ΔG vs sequence
            ax3 = fig.add_subplot(gs[0, 2])
            self._plot_delta_g(ax3, x_seq, delta_g_values)
            ax3.text(-0.1, 1.1, '(C)', transform=ax3.transAxes, fontsize=12, fontweight='bold')
            
            # Panel D: Activation energies
            ax4 = fig.add_subplot(gs[1, 0])
            self._plot_activation(ax4, x_seq, forward_activation, reverse_activation)
            ax4.text(-0.1, 1.1, '(D)', transform=ax4.transAxes, fontsize=12, fontweight='bold')
            
            # Panel E: Rates vs sequence
            ax5 = fig.add_subplot(gs[1, 1])
            self._plot_rates_sequence(ax5, x_seq, forward_rates, reverse_rates)
            ax5.text(-0.1, 1.1, '(E)', transform=ax5.transAxes, fontsize=12, fontweight='bold')
            
            # Panel F: Rate distributions by stacking
            ax6 = fig.add_subplot(gs[1, 2])
            self._plot_rate_distributions(
                ax6, 
                forward_rates, 
                reverse_rates, 
                pdb_file, 
                sequence_file, 
                topology
            )
            ax6.text(-0.1, 1.1, '(F)', transform=ax6.transAxes, fontsize=12, fontweight='bold')
            
            # Adjust layout and save/show
            plt.tight_layout()
            if save_path:
                plt.savefig(save_path, dpi=dpi, bbox_inches='tight')
            plt.show()
            
        def _plot_coupling(self, ax: plt.Axes, x_seq: np.ndarray, 
                          coupling_values: List[float]) -> None:
            """Plot electronic coupling values."""
            ax.plot(x_seq, coupling_values, 'o-', color=self.colors['single'])
            ax.set_xlabel('Sequence Position')
            ax.set_ylabel('Electronic Coupling (meV)')
            ax.set_title('Electronic Coupling vs. Sequence')
            
        def _plot_lambda(self, ax: plt.Axes, x_seq: np.ndarray, 
                        lambda_values: List[float]) -> None:
            """Plot reorganization energies."""
            ax.plot(x_seq, lambda_values, 'o-', color=self.colors['single'])
            ax.set_xlabel('Sequence Position')
            ax.set_ylabel('Reorganization Energy (eV)')
            ax.set_title('Reorganization Energy vs. Sequence')
            
        def _plot_delta_g(self, ax: plt.Axes, x_seq: np.ndarray, 
                         delta_g_values: List[float]) -> None:
            """Plot reaction free energies."""
            ax.plot(x_seq, delta_g_values, 'o-', 
                   label='Forward', color=self.colors['forward'])
            ax.plot(x_seq, [-dg for dg in delta_g_values], 'o-', 
                   label='Reverse', color=self.colors['reverse'])
            ax.set_xlabel('Sequence Position')
            ax.set_ylabel('ΔG (eV)')
            ax.set_title('Reaction Free Energy vs. Sequence')
            ax.legend()
            
        def _plot_activation(self, ax: plt.Axes, x_seq: np.ndarray,
                            forward_activation: List[float],
                            reverse_activation: List[float]) -> None:
            """Plot activation energies."""
            ax.plot(x_seq, forward_activation, 'o-', 
                   label='Forward', color=self.colors['forward'])
            ax.plot(x_seq, reverse_activation, 'o-', 
                   label='Reverse', color=self.colors['reverse'])
            ax.set_xlabel('Sequence Position')
            ax.set_ylabel('Activation Energy (eV)')
            ax.set_title('Activation Energy vs. Sequence')
            ax.legend()
            
        def _plot_rates_sequence(self, ax: plt.Axes, x_seq: np.ndarray,
                               forward_rates: List[float],
                               reverse_rates: List[float]) -> None:
            """Plot rates vs sequence."""
            ax.plot(x_seq, forward_rates, 'o-', 
                   label='Forward', color=self.colors['forward'])
            ax.plot(x_seq, reverse_rates, 'o-', 
                   label='Reverse', color=self.colors['reverse'])
            ax.set_xlabel('Sequence Position')
            ax.set_ylabel('Rate (s⁻¹)')
            ax.set_title('Electron Transfer Rates vs. Sequence')
            ax.set_yscale('log')
            ax.legend()
 
        def _plot_rate_distributions(self, ax: plt.Axes,
                                   forward_rates: List[float],
                                   reverse_rates: List[float],
                                   pdb_file: str,
                                   sequence_file: str = None,
                                   topology: str = 'linear') -> None:
            """Plot rate distributions by stacking type."""
            # Get structural classifications
            processor = PDBProcessor()
            geometry_data = processor.measure_heme_geometry(
                pdb_file, sequence_file, topology
            )
            
            # Print for debugging
            print("\nStacking types found:")
            for i, geom in enumerate(geometry_data):
                print(f"Pair {i+1}: {geom['stacking_type']}")
            
            # Organize rates by stacking type
            stacking_data = {
                'T-stacked': {'forward': [], 'reverse': []},
                'slip-stacked': {'forward': [], 'reverse': []},
                'co-planar': {'forward': [], 'reverse': []}
            }
            
            for i, geom in enumerate(geometry_data):
                stype = geom['stacking_type']  # Keep original case
                if stype not in stacking_data:
                    print(f"\nWarning: Unknown stacking type '{stype}'")
                    continue
                if i < len(forward_rates):
                    stacking_data[stype]['forward'].append(forward_rates[i])
                    stacking_data[stype]['reverse'].append(reverse_rates[i])
            
            # Create box plots
            positions = []
            data = []
            labels = []
            colors = []
            
            pos = 1
            for stype in ['slip-stacked', 'T-stacked', 'co-planar']:
                if stacking_data[stype]['forward']:
                    positions.extend([pos, pos + 0.5])
                    data.append(stacking_data[stype]['forward'])
                    data.append(stacking_data[stype]['reverse'])
                    labels.extend([f'{stype}\nforward', f'{stype}\nreverse'])
                    colors.extend([self.colors['forward'], self.colors['reverse']])
                pos += 2
                
            if data:  # Only create plot if we have data
                bp = ax.boxplot(data, positions=positions, patch_artist=True)
                
                # Style the box plots
                for patch, color in zip(bp['boxes'], colors):
                    patch.set_facecolor(color)
                    patch.set_alpha(0.6)
                
                # Add vertical lines between stacking types
                ymin, ymax = ax.get_ylim()
                ax.vlines([3, 5], ymin, ymax, linestyles='solid', colors='gray', alpha=0.5)
                
                ax.set_xticks(positions)
                ax.set_xticklabels(labels, rotation=45)
                ax.set_ylabel('Rate (s⁻¹)')
                ax.set_yscale('log')
                ax.set_title('Rate Distributions by Stacking Type')
                
            else:
                ax.text(0.5, 0.5, 'No structural data available',
                       ha='center', va='center')

    def __init__(self, 
                 interaction_manager: InteractionManager,
                 launch_dir: Path):
        """Initialize RateCalculator."""
        self.interaction_manager = interaction_manager
        self.launch_dir = Path(launch_dir)
        self.ee_dir = self.launch_dir / "EE"
        
        # Storage for computed values
        self.lambda_values = None
        self.delta_g_values = None
        self.coupling_values = None
        self.forward_activation = None
        self.reverse_activation = None
        self.forward_rates = None
        self.reverse_rates = None
        
        # Initialize plotter
        self.plotter = self.RatePlotter()
        
        # Physical constants
        self.constants = {
            'temperature': 300.0,  # Kelvin
            'pi': 3.141592654,
            'kB': 8.6173304E-5,    # eV/K
            'hbar': 6.582119514E-16  # eV·s
        }

    def _write_rates_report(
        self,
        lambda_values: List[float],
        delta_g_values: List[float],
        coupling_values: List[float],
        forward_activation: List[float],
        reverse_activation: List[float],
        forward_rates: List[float],
        reverse_rates: List[float]
    ) -> None:
        """
        Write detailed rates report with input parameters.
        
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

    def compute_marcus_rates(
        self, 
        lambda_values: List[float], 
        delta_g_values: List[float], 
        coupling_values: List[float],
        write_report: bool = True
    ) -> Tuple[List[float], List[float]]:
        """
        Compute electron transfer rates using Marcus theory.
        """
        # Store input values
        self.lambda_values = lambda_values
        self.delta_g_values = delta_g_values
        self.coupling_values = coupling_values
        
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
        
        # Store computed values
        self.forward_activation = forward_activation
        self.reverse_activation = reverse_activation
        self.forward_rates = forward_rates
        self.reverse_rates = reverse_rates
        
        # Write rates to file if requested
        if write_report:
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
        
    def plot_analysis(self, pdb_file: str, sequence_file: str = None, 
                     topology: str = 'linear', save_path: str = None,
                     dpi: int = 300) -> None:
        """Create comprehensive analysis plots."""
        if any(v is None for v in [
            self.lambda_values, self.delta_g_values, self.coupling_values,
            self.forward_activation, self.reverse_activation,
            self.forward_rates, self.reverse_rates
        ]):
            raise ValueError("No rate data available. Run compute_marcus_rates first.")
            
        self.plotter.create_analysis_figure(
            self.lambda_values,
            self.delta_g_values,
            self.coupling_values,
            self.forward_activation,
            self.reverse_activation,
            self.forward_rates,
            self.reverse_rates,
            pdb_file,
            sequence_file,
            topology,
            save_path,
            dpi
        )
