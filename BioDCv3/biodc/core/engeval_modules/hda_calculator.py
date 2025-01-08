"""
Electronic Coupling Calculator for Multi-Heme Systems

Calculates electronic coupling based on geometric configuration 
of heme groups using plane angle and distance calculations.
"""

from typing import List, Dict, Optional
from pathlib import Path
import logging

from biodc.utils.structure_analyzer import PDBProcessor
from biodc.utils.interaction import InteractionManager

class CouplingCalculator:
    """
    Calculator for electronic coupling in multi-heme systems.
    
    Supports computation of coupling values based on geometric 
    configuration of heme groups.
    """
    
    def __init__(self, 
                 interaction_manager: InteractionManager,
                 launch_dir: Path):
        """
        Initialize CouplingCalculator
        
        Args:
            interaction_manager: Manager for user interactions
            launch_dir: Project launch directory
        """
        self.interaction_manager = interaction_manager
        self.launch_dir = Path(launch_dir)
        self.ee_dir = self.launch_dir / "EE"
        self.structure_analyzer = PDBProcessor()
        
        # Coupling values in meV
        self.coupling_values = {
            'strong': 8.0,    # Slip-stacked or co-planar
            'weak': 2.0,      # T-stacked
            'co-planar': 1.0, # Co-planar with special significance
            'negligible': 0.1 # Very weak coupling
        }
    
    def compute_electronic_coupling(self, 
                                    sequence: List[int], 
                                    pdb_file: str) -> List[float]:
        """
        Compute electronic coupling for a sequence of hemes.
        
        Args:
            sequence: Sequence of heme residue indices
            pdb_file: Path to PDB structure file
        
        Returns:
            List of coupling values in meV
        """
        # Check for existing Hda.txt
        hda_file = self.ee_dir / "Hda.txt"
        
        if hda_file.exists():
            # Prompt user about existing file
            choice = self.interaction_manager.prompt(
                "existing_hda",
                "\nExisting Hda.txt found. What would you like to do?\n"
                "1) Use existing values\n"
                "2) Recompute from geometry\n"
                "3) Enter values manually\n"
                "Choice: ",
                choices=['1', '2', '3']
            )
            
            if choice == '1':
                return self._read_existing_hda()
            elif choice == '3':
                return self._get_manual_couplings(sequence)
        
        # If no file or user chooses recompute
        return self._compute_from_geometry(sequence, pdb_file)
    
    def _compute_from_geometry(self, 
                               sequence: List[int], 
                               pdb_file: str) -> List[float]:
        """
        Compute couplings from geometric configuration
        
        Args:
            sequence: Sequence of heme residue indices
            pdb_file: Path to PDB structure file
        
        Returns:
            List of coupling values in meV
        """
        # Ask user if they want to compute from geometry
        compute = self.interaction_manager.yes_no_prompt(
            "compute_geometry",
            "\nShould we estimate the electronic coupling from the geometry?"
        )
        
        if not compute:
            return self._get_manual_couplings(sequence)
        
        # Read PDB atoms
        atoms_dict = self.structure_analyzer.read_pdb_atoms(pdb_file)
        
        # Prepare detailed coupling data
        couplings_data = []
        for i in range(len(sequence) - 1):
            # Get coordinates for current and next heme
            heme1_coords = atoms_dict[sequence[i]]
            heme2_coords = atoms_dict[sequence[i+1]]
            
            # Calculate plane angle
            angle = self.structure_analyzer.calculate_plane_angle(
                heme1_coords, heme2_coords
            )
            
            # Calculate minimum edge-to-edge distance
            min_distance = self.structure_analyzer.calculate_min_distance(
                heme1_coords, heme2_coords
            )
            
            # Classify stacking
            stacking_type = self.structure_analyzer.classify_stacking(angle)
            
            # Estimate coupling based on stacking type
            if stacking_type == "co-planar":
                coupling = self.coupling_values['co-planar']
            elif stacking_type == "slip-stacked":
                coupling = self.coupling_values['strong']
            elif stacking_type == "T-stacked":
                coupling = self.coupling_values['weak']
            else:
                coupling = self.coupling_values['negligible']
            
            # Store detailed information
            couplings_data.append({
                'donor_heme': sequence[i],
                'acceptor_heme': sequence[i+1],
                'angle': angle,
                'distance': min_distance,
                'stacking_type': stacking_type,
                'coupling': coupling
            })
        
        # Write comprehensive report
        self._write_coupling_report(couplings_data)
        
        # Return just the coupling values
        return [data['coupling'] for data in couplings_data]
    
    def _get_manual_couplings(self, sequence: List[int]) -> List[float]:
        """
        Get manually entered coupling values
        
        Args:
            sequence: Sequence of heme residue indices
        
        Returns:
            List of manually entered coupling values
        """
        couplings = []
        self.interaction_manager.prompt(
            "manual_coupling_entry",
            "\nEnter electronic coupling (Hda) for each transfer step:",
            choices=None,
            allow_empty=True
        )
        
        for i in range(len(sequence) - 1):
            value = self.interaction_manager.prompt(
                f"hda_manual_{i}",
                f"Step {i+1} (Heme {sequence[i]} → Heme {sequence[i+1]}) (meV): ",
                input_type=float
            )
            couplings.append(value)
        
        # Write manual values to file
        self._write_manual_hda(sequence, couplings)
        
        return couplings
    
    def _read_existing_hda(self) -> List[float]:
        """
        Read existing Hda values from file
        
        Supports reading from both computational and manual entry formats
        
        Returns:
            List of coupling values
        """
        couplings = []
        hda_file = self.ee_dir / "Hda.txt"
        
        with open(hda_file, 'r') as f:
            for line in f:
                # Check computational mode format
                if 'Hda coupling' in line:
                    coupling = float(line.split('=')[-1].strip().split()[0])
                    couplings.append(coupling)
                
                # Check manual mode format
                elif 'Hda coupling    =' in line:
                    coupling = float(line.split('=')[-1].strip())
                    couplings.append(coupling)
        
        # Prompt user to review values
        for i, coupling in enumerate(couplings):
            self.interaction_manager.prompt(
                f"existing_hda_{i}",
                f"Step {i+1} Hda: {coupling:.3f} meV",
                choices=None,
                allow_empty=True
            )
        
        return couplings
    
    def _write_coupling_report(self, couplings_data: List[Dict]):
        """
        Generate a detailed coupling report
        
        Args:
            couplings_data: List of detailed coupling calculations
        """
        with open(self.ee_dir / "Hda.txt", 'w') as f:
            f.write("Heme Plane Angle and Electronic Coupling Analysis\n")
            f.write("=" * 50 + "\n\n")
            
            for data in couplings_data:
                f.write(
                    f"Hda(HEM-{data['donor_heme']} <-> HEM-{data['acceptor_heme']}):\n"
                    f"  Plane angle     = {data['angle']:10.3f} deg\n"
                    f"  Edge-to-edge    = {data['distance']:10.3f} Å\n"
                    f"  Stacking type   = {data['stacking_type']}\n"
                    f"  Hda coupling    = {data['coupling']:6.3f} meV\n\n"
                )
    
    def _write_manual_hda(self, sequence: List[int], couplings: List[float]):
        """
        Write manually entered coupling values to file
        
        Args:
            sequence: Sequence of heme residue indices
            couplings: Manually entered coupling values
        """
        with open(self.ee_dir / "Hda.txt", 'w') as f:
            f.write("Manually Entered Electronic Coupling Values\n")
            f.write("=" * 45 + "\n\n")
            
            for i, (donor, acceptor, coupling) in enumerate(
                zip(sequence[:-1], sequence[1:], couplings)
            ):
                f.write(
                    f"Hda(HEM-{donor} <-> HEM-{acceptor}):\n"
                    f"  Hda coupling    = {coupling:6.3f} meV\n\n"
                )
