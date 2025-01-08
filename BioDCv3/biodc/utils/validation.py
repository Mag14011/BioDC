# biodc/utils/validation.py

"""
Validation utilities for BioDC package.
Provides functions for validating input parameters, 
file integrity, and system requirements.
"""

import os
import sys
from pathlib import Path
from typing import Dict, Any, Optional, List

import subprocess
from Bio import PDB

class BioDCValidator:
    """
    Comprehensive validation class for BioDC package requirements.
    """
    
    @staticmethod
    def check_system_dependencies() -> Dict[str, bool]:
        """
        Check for required external dependencies.
        
        Returns:
            Dictionary of dependency status
        """
        dependencies = {
            'VMD': ['vmd', '--version'],
            'AmberTools': ['tleap', '-h'],
            'CPPTRAJ': ['cpptraj', '-h'],
            'SANDER': ['sander', '-h'],
            'Python': [sys.executable, '--version']
        }
        
        results = {}
        for name, cmd in dependencies.items():
            try:
                subprocess.run(
                    cmd, 
                    stdout=subprocess.PIPE, 
                    stderr=subprocess.PIPE, 
                    text=True
                )
                results[name] = True
            except (FileNotFoundError, subprocess.CalledProcessError):
                results[name] = False
        
        return results

    @staticmethod
    def validate_pdb(pdb_path: str) -> Dict[str, Any]:
        """
        Validate PDB file integrity and characteristics.
        
        Args:
            pdb_path: Path to PDB file
        
        Returns:
            Dictionary of validation results
        """
        results = {
            'file_exists': False,
            'is_valid_pdb': False,
            'total_atoms': 0,
            'total_residues': 0,
            'consecutive_numbering': False,
            'contains_heme': False
        }
        
        # Check file existence
        if not os.path.exists(pdb_path):
            return results
        results['file_exists'] = True
        
        try:
            # Use BioPython to parse PDB
            parser = PDB.PDBParser(QUIET=True)
            structure = parser.get_structure('protein', pdb_path)
            
            # Count atoms and residues
            atoms = list(structure.get_atoms())
            residues = list(structure.get_residues())
            
            results['total_atoms'] = len(atoms)
            results['total_residues'] = len(residues)
            results['is_valid_pdb'] = True
            
            # Check for consecutive residue numbering
            res_ids = [res.id[1] for res in residues]
            results['consecutive_numbering'] = (
                len(res_ids) == len(set(res_ids)) and 
                res_ids == list(range(min(res_ids), max(res_ids) + 1))
            )
            
            # Check for heme
            results['contains_heme'] = any(
                res.resname in ['HEM', 'HEME'] for res in residues
            )
            
        except Exception as e:
            print(f"Error parsing PDB: {e}")
            results['is_valid_pdb'] = False
        
        return results

    @staticmethod
    def validate_configuration(config: Dict[str, Any]) -> List[str]:
        """
        Validate configuration dictionary for required keys.
        
        Args:
            config: Configuration dictionary
        
        Returns:
            List of missing or invalid configuration keys
        """
        required_keys = [
            'pdb_file',
            'forcefield_dir',
            'solvent_type'
        ]
        
        missing_keys = [
            key for key in required_keys 
            if key not in config or not config[key]
        ]
        
        # Additional validation for specific keys
        if 'solvent_type' in config and config['solvent_type'] not in ['explicit', 'implicit']:
            missing_keys.append('invalid_solvent_type')
        
        return missing_keys

def validate_biodc_setup(pdb_path: Optional[str] = None, 
                        config: Optional[Dict[str, Any]] = None) -> bool:
    """
    Comprehensive validation of BioDC setup.
    
    Args:
        pdb_path: Optional path to PDB file
        config: Optional configuration dictionary
    
    Returns:
        Boolean indicating whether setup is valid
    """
    validator = BioDCValidator()
    
    # Check system dependencies
    dependencies = validator.check_system_dependencies()
    print("\nSystem Dependencies:")
    for name, status in dependencies.items():
        print(f"{name}: {'✓' if status else '✗'}")
    
    # Validate PDB if path provided
    if pdb_path:
        pdb_validation = validator.validate_pdb(pdb_path)
        print("\nPDB Validation:")
        for key, value in pdb_validation.items():
            print(f"{key}: {value}")
    
    # Validate configuration if provided
    if config:
        config_errors = validator.validate_configuration(config)
        print("\nConfiguration Validation:")
        if config_errors:
            print("Missing or invalid keys:", config_errors)
            return False
        else:
            print("Configuration: Valid")
    
    # Return overall status
    return all(dependencies.values()) and \
           (not pdb_path or pdb_validation['is_valid_pdb']) and \
           (not config or not config_errors)
