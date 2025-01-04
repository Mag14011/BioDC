# biodc/core/preparation/ligand_detection.py

"""
Module for detecting metal ligands in protein structures.
Currently focused on heme iron coordination but structured for potential expansion.
"""

from pathlib import Path
from typing import Dict, List, Set, Tuple
from dataclasses import dataclass
from Bio import PDB
import numpy as np

@dataclass
class MetalCoordination:
    """Container for metal coordination information."""
    metal_id: int
    metal_name: str  # e.g., 'FE' for heme iron
    residue_id: int
    residue_name: str
    atom_name: str
    distance: float

class LigandDetector:
    """Analyzes metal coordination in protein structures."""
    
    # Coordinating atoms to look for in residues
    COORDINATING_ATOMS = {
        'HIS': ['NE2', 'ND1'],  # His coordinates through nitrogens
        'MET': ['SD'],          # Met coordinates through sulfur
        'CYS': ['SG'],          # Cys coordinates through sulfur
        'TYR': ['OH'],          # Tyr coordinates through oxygen
        'ASP': ['OD1', 'OD2'],  # Asp coordinates through carboxylate oxygens
        'GLU': ['OE1', 'OE2'],  # Glu coordinates through carboxylate oxygens
        'ASN': ['OD1', 'ND2'],  # Asn can coordinate through O or N
        'GLN': ['OE1', 'NE2']   # Gln can coordinate through O or N
    }
    
    def __init__(self, pdb_path: str):
        """Initialize with PDB file path."""
        parser = PDB.PDBParser(QUIET=True)
        self.structure = parser.get_structure('protein', pdb_path)
        self.coordinated_residues: List[MetalCoordination] = []

    def find_thioether_bonds(self, heme_residues: List[str] = ['HEM', 'HEME', 'HEC']) -> Dict[int, List[int]]:
        """
        Find Cys residues making thioether bonds to heme vinyl groups.
        Searches incrementally from 2-4Å until an SG atom is found for each vinyl carbon.

        Args:
            heme_residues: List of residue names that could be hemes

        Returns:
            Dict mapping heme IDs to lists of Cys residue IDs that form thioether bonds
        """
        thioether_bonds = {}

        for model in self.structure:
            for chain in model:
                for heme in chain:
                    if heme.get_resname() in heme_residues:
                        heme_id = heme.get_id()[1]
                        thioether_bonds[heme_id] = []

                        # Get vinyl carbons
                        cab = None
                        cac = None
                        for atom in heme:
                            if atom.get_name() == 'CAB':
                                cab = atom
                            elif atom.get_name() == 'CAC':
                                cac = atom

                        if not (cab and cac):
                            continue

                        # Search for each vinyl carbon separately
                        for vinyl_carbon in [cab, cac]:
                            found_cys = False
                            # Increment search distance until we find a Cys
                            for distance in np.arange(2.0, 4.0, 0.1):
                                if found_cys:
                                    break

                                for res in chain:
                                    if res.get_resname() == 'CYS':
                                        if 'SG' in res:
                                            sg_atom = res['SG']
                                            if sg_atom - vinyl_carbon <= distance:
                                                thioether_bonds[heme_id].append(res.get_id()[1])
                                                found_cys = True
                                                break

        return thioether_bonds

    def find_metal_coordinations(self, metal_name: str = 'FE',
                            residue_names: List[str] = ['HEM', 'HEME', 'HEC']) -> List[MetalCoordination]:
        """
        Find residues coordinating to Fe centers.
        Searches incrementally from 2-4Å until ligands are found.
        """
        # Clear existing coordinations
        self.coordinated_residues = []

        # Find all metal atoms
        metal_atoms = []
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    if residue.get_resname() in residue_names:
                        for atom in residue:
                            if atom.get_name() == metal_name:
                                metal_atoms.append((residue.get_id()[1], atom))

        # For each metal center, find its coordinating residues
        for metal_id, metal_atom in metal_atoms:
            ligands_found = 0

            # Search with increasing distance until we find at least two ligands
            for distance in np.arange(2.0, 4.0, 0.1):
                if ligands_found >= 2:
                    break

                # Check all potential coordinating residues
                for model in self.structure:
                    for chain in model:
                        for residue in chain:
                            res_name = residue.get_resname()
                            if res_name in self.COORDINATING_ATOMS:
                                res_id = residue.get_id()[1]

                                # Skip if we already found this residue
                                if any(coord.residue_id == res_id
                                    for coord in self.coordinated_residues):
                                    continue

                                # Check coordinating atoms
                                for atom_name in self.COORDINATING_ATOMS[res_name]:
                                    if atom_name in residue:
                                        atom = residue[atom_name]
                                        coord_distance = atom - metal_atom

                                        if coord_distance <= distance:
                                            self.coordinated_residues.append(
                                                MetalCoordination(
                                                    metal_id=metal_id,
                                                    metal_name=metal_name,
                                                    residue_id=res_id,
                                                    residue_name=res_name,
                                                    atom_name=atom_name,
                                                    distance=coord_distance
                                                )
                                            )
                                            ligands_found += 1
                                            break  # Move to next residue

        return self.coordinated_residues

    def print_coordination_info(self):
        """Print information about found coordinations."""
        print("\nMetal coordination analysis:")
        for coord in self.coordinated_residues:
            print(f" {coord.metal_name} {coord.metal_id}: {coord.residue_name} {coord.residue_id} "
                  f"coordinates through {coord.atom_name} "
                  f"(distance: {coord.distance:.2f}Å)")
    
    def get_excluded_residues(self) -> Dict[str, Set[int]]:
        """
        Get sets of residue IDs that should be excluded from titration
        due to metal coordination.
        """
        excluded = {
            'HIS': set(),  # Titratable
            'ASP': set(),  # Titratable
            'GLU': set(),  # Titratable
            'TYR': set()   # Titratable
        }
        
        for coord in self.coordinated_residues:
            if coord.residue_name in excluded:
                excluded[coord.residue_name].add(coord.residue_id)
        
        return excluded

def ligand_detection(pdb_path: str) -> Dict[str, set]:
    """
    Detect metal-coordinating residues in PDB structure.

    Args:
        pdb_path: Path to PDB file

    Returns:
        Dictionary of residue types to sets of residue IDs that coordinate metals
    """
    detector = LigandDetector(pdb_path)
    detector.find_metal_coordinations()
    return detector.get_excluded_residues()

