import os
import sys
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
from matplotlib.ticker import AutoMinorLocator
from dataclasses import dataclass
from typing import Optional, Dict, List, Tuple
################################################################################################################################################
import derrida
import blumberger
################################################################################################################################################

@dataclass
class PhysicalConstants:
    """Physical constants used in calculations"""
    r: float = 0.75E-7      # radius in cm
    c: float = 3.00E10      # speed of light in cm/s
    e: float = 1.602E-19    # elementary charge in C
    kb: float = 1.38E-23    # Boltzmann constant in J/K
    hbar: float = 1.05E-34  # reduced Planck constant in J·s
    deltabar: float = 3     # reorganization energy in cm^-1

@dataclass
class TransportProperties:
    """Container for transport properties"""
    conductivity: float    # S/cm
    mobility: float        # cm²/(V·s)
    diffusion_constant: float  # cm²/s
    hopping_rate: float   # s⁻¹
    conductance: float    # S

class PDBProcessor:
    def __init__(self):
        self.heme_atoms = [
            'FE', 'NA', 'C1A', 'C2A', 'C3A', 'C4A', 'CHB', 'C1B', 'NB', 
            'C2B', 'C3B', 'C4B', 'CHC', 'C1C', 'NC', 'C2C', 'C3C', 'C4C', 
            'CHD', 'C1D', 'ND', 'C2D', 'C3D', 'C4D', 'CHA'
        ]

    def read_linearized_sequence(self, filename: str) -> List[Tuple[int, int]]:
        """Read the linearized heme sequence file"""
        pairs = []
        with open(filename, 'r') as f:
            lines = f.readlines()
            for line in lines:
                try:
                    _, resid = map(int, line.strip().split())
                    pairs.append(resid)
                except ValueError:
                    continue
        return list(zip(pairs[:-1], pairs[1:]))

    def read_pdb_atoms(self, pdb_file: str) -> Dict[int, Dict[str, np.ndarray]]:
        """Read atom coordinates from PDB file"""
        atoms_dict = {}
        
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM  ') or line.startswith('HETATM'):
                    try:
                        atom_name = line[12:16].strip()
                        resid = int(line[22:26])
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        
                        if atom_name in self.heme_atoms:
                            if resid not in atoms_dict:
                                atoms_dict[resid] = {}
                            atoms_dict[resid][atom_name] = np.array([x, y, z])
                    except (ValueError, IndexError):
                        continue
        
        return atoms_dict

    def calculate_min_distance(self, coords1: Dict[str, np.ndarray], 
                             coords2: Dict[str, np.ndarray]) -> float:
        """Calculate minimum distance between two sets of coordinates"""
        min_dist = float('inf')
        
        for atom1 in self.heme_atoms:
            if atom1 not in coords1:
                continue
            for atom2 in self.heme_atoms:
                if atom2 not in coords2:
                    continue
                    
                dist = np.linalg.norm(coords1[atom1] - coords2[atom2])
                min_dist = min(min_dist, dist)
        
        return min_dist

    def measure_avg_dist(self, pdb_file: str, sequence_file: str) -> List[float]:
        """Calculate average distances between consecutive heme pairs"""
        residue_pairs = self.read_linearized_sequence(sequence_file)
        atoms_dict = self.read_pdb_atoms(pdb_file)
        
        distances = []
        for res1, res2 in residue_pairs:
            if res1 in atoms_dict and res2 in atoms_dict:
                min_dist = self.calculate_min_distance(atoms_dict[res1], atoms_dict[res2])
                distances.append(min_dist)
            else:
                print(f"Warning: Missing residue {res1} or {res2} in PDB file")
        
        return distances

    def measure_subunit_length(self, pdb_file: str, sequence_file: str) -> float:
        """Calculate distance between first and second-to-last FE atoms"""
        residue_pairs = self.read_linearized_sequence(sequence_file)
        atoms_dict = self.read_pdb_atoms(pdb_file)
        
        first_res = residue_pairs[0][0]
        second_to_last_res = residue_pairs[-1][0]
        
        if 'FE' in atoms_dict.get(first_res, {}) and 'FE' in atoms_dict.get(second_to_last_res, {}):
            return np.linalg.norm(
                atoms_dict[first_res]['FE'] - atoms_dict[second_to_last_res]['FE']
            )
        else:
            raise ValueError(f"Missing FE atoms in residues {first_res} or {second_to_last_res}")

class RateProcessor:
    def __init__(self, rates_file: str):
        self.rates_file = rates_file
        self.rates_data = None
        self._load_rates()

    def _load_rates(self):
        """Load rates from rates file"""
        try:
            self.rates_data = pd.read_csv(self.rates_file)
            self.ketf = self.rates_data['ketf'].tolist()
            self.ketb = self.rates_data['ketb'].tolist()
        except FileNotFoundError:
            print(f"Error: Rates file {self.rates_file} not found")
            raise
        except KeyError:
            print(f"Error: {self.rates_file} must contain 'ketf' and 'ketb' columns")
            raise

    def compute_flux(self) -> tuple:
        """Compute forward, backward, and average flux"""
        Jf, Jb = blumberger.flux(self.ketf, self.ketb)
        Javg = (Jf + Jb) / 2
        
        print(f"Forward Flux: {Jf:.2E}")
        print(f"Reverse Flux: {Jb:.2E}")
        print(f"Average forward/backward Flux: {Javg:.2E}")
        
        return Jf, Jb, Javg

    def compute_diffusion_coefficient(self, avg_heme_spacing: float) -> float:
        """Compute diffusion coefficient"""
        V, D = derrida.VD(self.ketf, self.ketb)
        D_final = D * (avg_heme_spacing**2)
        
        print(f"  Diffusion constant = {D_final:E} cm^2/S")
        with open('D.txt', 'w') as f:
            print(f"Diffusion constant = {D_final:E} (cm^2/S)", file=f)
        
        return D_final

class RedoxCurrentAnalyzer:
    def __init__(self, rates_file: str):
        self.constants = PhysicalConstants()
        self.pdb_processor = PDBProcessor()
        self.rate_processor = RateProcessor(rates_file)
        plt.style.use('default')

    def compute_properties_from_flux(self, Javg: float, T: float, crgden: float, 
                                   ahs: float, csa: float, lw: float) -> TransportProperties:
        """Compute transport properties from average charge flux"""
        current = Javg * self.constants.e
        V_typical = 0.5  
        conductance = current / V_typical
        conductivity = conductance * (lw/csa)
        mobility = conductivity / (crgden * self.constants.e)
        diffusion_constant = mobility * self.constants.kb * T / self.constants.e
        hopping_rate = diffusion_constant / (ahs * ahs)
        
        return TransportProperties(
            conductivity=conductivity,
            mobility=mobility,
            diffusion_constant=diffusion_constant,
            hopping_rate=hopping_rate,
            conductance=conductance
        )

    def compute_properties_from_diffusion(self, D: float, T: float, crgden: float, 
                                        ahs: float, csa: float, lw: float) -> TransportProperties:
        """Compute transport properties from diffusion coefficient"""
        mobility = (self.constants.e * D) / (self.constants.kb * T)
        conductivity = crgden * self.constants.e * mobility
        conductance = conductivity * csa / lw
        hopping_rate = D / (ahs * ahs)
        
        return TransportProperties(
            conductivity=conductivity,
            mobility=mobility,
            diffusion_constant=D,
            hopping_rate=hopping_rate,
            conductance=conductance
        )

    def _create_iv_plot(self, voltages: np.ndarray, currents_dict: Dict[str, list]) -> None:
        """Create current vs voltage plot with enhanced visibility"""
        fig, ax1 = plt.subplots(figsize=(3.3, 3.3))

        # Colors
        purple_color = '#800080'  # for flux-derived
        darker_orange = '#FF4500'  # DarkOrange for diffusion-derived
        light_green = '#90EE90'    # Light green for background region

        # Set up axes first
        ax1.set_xlabel('Voltage (V)')

        # Add light green background for the region -0.2 to +0.2 V
        ax1.axvspan(-0.2, 0.2, color=light_green, alpha=0.3, zorder=0)
    
        # Add "Bio-Relevant" annotation with relative positioning
        # Using transform=ax1.transAxes makes coordinates relative to axes (0-1)
        ax1.text(0.5, 0.80, 'Bio-Relevant', 
                color=light_green,
                path_effects=[path_effects.withStroke(linewidth=2, foreground='black')],
                ha='center',
                va='center',
                rotation=0,
                fontsize=10,
                zorder=1,
                transform=ax1.transAxes)  # This makes coordinates relative to the axes

        # Create second axis for flux data
        ax2 = ax1.twinx()

        # Plot flux-derived current first (background)
        if 'Flux' in currents_dict:
            ax1.plot(voltages, np.array(currents_dict['Flux']) * 1e3,
                    color=purple_color, linestyle='-', linewidth=1.5)
            ax1.set_ylabel('Flux-derived Current (fA)', color=purple_color)
            ax1.tick_params(axis='y', labelcolor=purple_color, direction='in')
            ax1.yaxis.set_minor_locator(AutoMinorLocator())
            ax1.tick_params(which='minor', direction='in')

        # Clear the current axis to ensure proper layering
        plt.sca(ax1)

        # Plot diffusion-derived current last (foreground)
        if 'Diffusion' in currents_dict:
            ax2.plot(voltages, np.array(currents_dict['Diffusion']) * 1e3,
                    color=darker_orange, linestyle='--', linewidth=2,
                    dashes=(5, 3))  # Longer dashes, shorter spaces
        
            # Rotate right y-axis label and add padding
            ax2.set_ylabel('Diffusion-derived Current (fA)', 
                        color=darker_orange,
                        rotation=270,
                        labelpad=15)  # Increase padding to prevent overlap
        
            ax2.tick_params(axis='y', labelcolor=darker_orange, direction='in')
            ax2.tick_params(axis='x', direction='in')
            ax2.xaxis.set_minor_locator(AutoMinorLocator())
            ax2.yaxis.set_minor_locator(AutoMinorLocator())
            ax2.tick_params(which='minor', direction='in')

        plt.tight_layout()
        plt.savefig('IV_curve.png', dpi=600, bbox_inches='tight')
        plt.close()


    def perform_voltage_sweep(self, properties_dict: Dict[str, TransportProperties], 
                            Gexp: Optional[float], csa: float, lw: float) -> None:
        """Perform voltage sweep using different conductances"""
        voltages = np.arange(-0.5, 0.55, 0.05)
        
        currents_dict = {}
        for source, props in properties_dict.items():
            currents_dict[source] = [props.conductance * V * 1E12 for V in voltages]
            
        if Gexp is not None:
            currents_dict['Experimental'] = [Gexp * V * 1E12 for V in voltages]
        
        self._create_iv_plot(voltages, currents_dict)
        
        header = " Voltage (V)"
        if Gexp is not None:
            header += " | Exp. Current (pA)"
        for source in properties_dict.keys():
            header += f" | {source} Current (pA)"
        
        with open("RedoxCurrentAnalysis.txt", 'w') as f:
            print(header, file=f)
            for i, V in enumerate(voltages):
                line = f" {V:11.3f}"
                if Gexp is not None:
                    line += f" | {currents_dict['Experimental'][i]:17.3f}"
                for source in properties_dict.keys():
                    line += f" | {currents_dict[source][i]:21.3E}"
                print(line, file=f)

    def compute_redox_current(self, pdb_file: str, sequence_file: str):
        """Main function to compute redox current and related properties"""
        print(" Please provide the following parameters: ")
        
        T = self._get_float_input("Temperature (K)? ")
        cps = self._get_int_input("Number of charges per subunit? ")
        
        distances = self.pdb_processor.measure_avg_dist(pdb_file, sequence_file)
        ahs = (round(sum(distances)/len(distances), 3)) * 1E-8
        print(f"  Average heme spacing is {ahs} cm")
        
        Ncnt = self._get_float_input("Fewest number of contacts at either protein-electrode interface? ")
        
        try:
            Jf, Jb, Javg = self.rate_processor.compute_flux()
            D = self.rate_processor.compute_diffusion_coefficient(ahs)
        except Exception as e:
            print(f"Error in rate calculations: {str(e)}")
            sys.exit(1)
        
        SubunitLength = self.pdb_processor.measure_subunit_length(pdb_file, sequence_file) * 1E-8
        print(f"\nThe subunit length measured is {SubunitLength:.2E} cm\n")
        
        lw = self._get_float_input("Length of wire (cm)? ")
        
        cpsul = cps / SubunitLength
        csa = math.pi * (self.constants.r ** 2)
        crgden = cpsul / csa

        # Calculate properties from both methods
        flux_props = self.compute_properties_from_flux(
            Javg=Javg, T=T, crgden=crgden, ahs=ahs, csa=csa, lw=lw
        )

        diffusion_props = self.compute_properties_from_diffusion(
            D=D, T=T, crgden=crgden, ahs=ahs, csa=csa, lw=lw
        )

        # Print results from both methods
        print("\nTransport properties from charge flux:")
        self._print_properties(flux_props)

        print("\nTransport properties from diffusion coefficient:")
        self._print_properties(diffusion_props)

        # Get experimental conductance if known
        Gexp = None
        if self._get_yes_no_input("Is the experimental conductance known (yes/no)? "):
            Gexp = self._get_float_input("Experimental Conductance (S)? ")

        # Perform voltage sweep with all available methods
        properties_dict = {
            'Flux': flux_props,
            'Diffusion': diffusion_props
        }

        self.perform_voltage_sweep(properties_dict, Gexp, csa, lw)

    def _print_properties(self, props: TransportProperties) -> None:
        """Helper method to print transport properties"""
        print(f"Conductance: {props.conductance:.2E} S")
        print(f"Conductivity: {props.conductivity:.2E} S/cm")
        print(f"Mobility: {props.mobility:.2E} cm²/(V·s)")
        print(f"Diffusion constant: {props.diffusion_constant:.2E} cm²/s")
        print(f"Hopping rate: {props.hopping_rate:.2E} s⁻¹")

    def _get_float_input(self, prompt: str) -> float:
        """Helper method to get float input with error handling"""
        while True:
            try:
                return float(input(f"  {prompt}"))
            except ValueError:
                print(" Your entry must be a floating-point number.")

    def _get_int_input(self, prompt: str) -> int:
        """Helper method to get integer input with error handling"""
        while True:
            try:
                return int(input(f"  {prompt}"))
            except ValueError:
                print(" Your entry must be an integer.")

    def _get_yes_no_input(self, prompt: str) -> bool:
        """Helper method to get yes/no input"""
        while True:
            response = input(prompt).lower()
            if response in ['yes', 'y']:
                return True
            elif response in ['no', 'n']:
                return False
            print(" Please answer 'yes' or 'no'.")

################################################################################################################################################
# Main execution
if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <pdb_file> <sequence_file> <rates_file>")
        sys.exit(1)

    pdb_file = sys.argv[1]
    sequence_file = sys.argv[2]
    rates_file = sys.argv[3]

    # Validate input files
    for file_path in [pdb_file, sequence_file, rates_file]:
        if not os.path.exists(file_path):
            print(f"Error: File {file_path} not found")
            sys.exit(1)

    # Run analysis
    analyzer = RedoxCurrentAnalyzer(rates_file)
    analyzer.compute_redox_current(pdb_file, sequence_file)
################################################################################################################################################
