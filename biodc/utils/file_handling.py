# biodc/utils/file_handling.py
"""
File handling utilities for BioDC package.
Provides functions for reading configuration and input files.
"""

import os
from pathlib import Path
from typing import Dict, Any, Optional

def read_input(config_path: Optional[str] = None) -> Dict[str, str]:
    """
    Read input configuration from a key-value formatted file.
    
    Args:
        config_path: Path to configuration file
    
    Returns:
        Dictionary of configuration parameters
    """
    if not config_path:
        return {}
    
    # Ensure config path exists
    if not os.path.exists(config_path):
        return {}
    
    input_dict = {}
    try:
        with open(config_path, 'r') as f:
            for line in f:
                # Strip whitespace and skip empty lines
                line = line.strip()
                if not line or '=' not in line:
                    continue
                
                # Split on first '=' and strip whitespace
                key, value = line.split('=', 1)
                input_dict[key.strip()] = value.strip()
    except (IOError, PermissionError) as e:
        print(f"Error reading configuration file: {e}")
    
    return input_dict

def write_configuration(config: Dict[str, Any], 
                        output_path: str) -> None:
    """
    Write configuration to a key-value formatted file.
    
    Args:
        config: Configuration dictionary
        output_path: Path to write the configuration file
    """
    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Write configuration in key-value format
    with open(output_path, 'w') as f:
        for key, value in config.items():
            f.write(f"{key} = {value}\n")

def validate_file_path(file_path: str, 
                       must_exist: bool = True, 
                       is_directory: bool = False) -> Path:
    """
    Validate and normalize file or directory path.
    
    Args:
        file_path: Path to validate
        must_exist: Whether the path must exist
        is_directory: Whether the path should be a directory
    
    Returns:
        Normalized Path object
    
    Raises:
        FileNotFoundError: If file/directory does not exist when must_exist is True
        NotADirectoryError: If is_directory is True but path is not a directory
        ValueError: If path is invalid
    """
    # Normalize path
    path = Path(file_path).expanduser().resolve()
    
    # Check existence if required
    if must_exist and not path.exists():
        raise FileNotFoundError(f"Path does not exist: {path}")
    
    # Check if it's a directory when required
    if is_directory and not path.is_dir():
        raise NotADirectoryError(f"Path is not a directory: {path}")
    
    return path