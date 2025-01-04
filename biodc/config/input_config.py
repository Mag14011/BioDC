# biodc/config/input_config.py
"""
Configuration schema and validation for BioDC input parameters.
"""

from typing import Dict, Any, Optional
from pydantic import BaseModel, Field, validator

class BioDCInputConfig(BaseModel):
    """
    Comprehensive input configuration for BioDC.
    Validates and provides type checking for input parameters.
    """
    
    # Structure Preparation Parameters
    pdb_file: Optional[str] = Field(
        None, 
        description="Path to input PDB file"
    )
    
    # Disulfide Selection
    select_disulfides: Optional[bool] = Field(
        False, 
        description="Whether to select disulfide linkages"
    )
    disulfide_pairs: Optional[list] = Field(
        None, 
        description="List of disulfide pairs to modify"
    )
    
    # Mutation Selection
    select_mutations: Optional[bool] = Field(
        False, 
        description="Whether to perform mutations"
    )
    mutation_details: Optional[Dict[str, str]] = Field(
        None, 
        description="Specific mutation details"
    )
    
    # pH Active Sites
    constant_ph_dynamics: Optional[bool] = Field(
        False, 
        description="Run molecular dynamics with titratable residues"
    )
    ph_active_residues: Optional[Dict[str, list]] = Field(
        None, 
        description="Residue IDs for pH-active sites"
    )
    
    # Workflow Selection
    workflow_division: Optional[str] = Field(
        None, 
        description="Selected workflow division (0/1/2/3)",
        pattern=r'^[0-3]$'
    )
    
    @validator('workflow_division')
    def validate_workflow_division(cls, v):
        """Validate workflow division input"""
        if v not in ['0', '1', '2', '3']:
            raise ValueError("Workflow division must be 0, 1, 2, or 3")
        return v
    
    @validator('disulfide_pairs', always=True)
    def validate_disulfide_pairs(cls, v, values):
        """
        Validate disulfide pairs when disulfide selection is True
        """
        if values.get('select_disulfides') and not v:
            raise ValueError("Disulfide pairs must be specified when select_disulfides is True")
        return v

def validate_input_config(input_dict: Dict[str, Any]) -> BioDCInputConfig:
    """
    Validate and create a BioDC input configuration.
    
    Args:
        input_dict: Dictionary of input parameters
    
    Returns:
        Validated BioDCInputConfig instance
    """
    try:
        return BioDCInputConfig(**input_dict)
    except Exception as e:
        print(f"Input configuration validation error: {e}")
        raise