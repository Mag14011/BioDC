"""
Configuration management for BioDC package.
"""

from typing import Dict, Optional

class Config:
    """
    Configuration handler for BioDC.
    
    Manages input parameters and provides validation.
    """
    
    def __init__(self, config: Optional[Dict] = None):
        """
        Initialize configuration.
        
        Args:
            config: Optional dictionary of configuration parameters
        """
        # Default configuration
        self._config = config or {}
    
    def get(self, key: str, default=None):
        """
        Retrieve a configuration value.
        
        Args:
            key: Configuration parameter key
            default: Default value if key not found
        
        Returns:
            Configuration value or default
        """
        return self._config.get(key, default)
    
    def __getitem__(self, key: str):
        """
        Get configuration value using dictionary-like access.
        
        Raises:
            KeyError if key not found
        """
        return self._config[key]
