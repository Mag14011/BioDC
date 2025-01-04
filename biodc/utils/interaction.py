# biodc/utils/interaction.py

import os
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, Optional, Union
import logging

logger = logging.getLogger("biodc")

class InteractionManager:
    def __init__(self,
                launch_dir: Optional[Union[str, Path]] = None,
                input_dict: Optional[Dict] = None):
        self.launch_dir = Path(launch_dir) if launch_dir else Path.cwd()
        self.launch_dir.mkdir(parents=True, exist_ok=True)

        self.input_file = self.launch_dir / "input.txt"
        self.interactive_log = self.launch_dir / "InteractiveInput.txt"
        self.prior_file = self.launch_dir / "PriorInteractiveInput.txt"

        # Initialize input sources with priority
        self.input_dict = {}
        if self.input_file.exists():
            self.input_dict.update(self._read_config_file(self.input_file))
        if input_dict:
            self.input_dict.update(input_dict)

        self.recorded_interactions = {}
        self._setup_logging()

    def _setup_logging(self):
        # Ensure the interactive log file exists and is writable
        self.interactive_log.touch()
        self._log_file = self.interactive_log.open('a')

    def _read_config_file(self, file_path: Path) -> Dict[str, str]:
        config = {}
        try:
            with open(file_path, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line and '=' in line:
                        key, value = line.split('=', 1)
                        config[key.strip()] = value.strip()
            logger.info(f"Loaded configuration from {file_path}")
        except Exception as e:
            logger.error(f"Error reading config file {file_path}: {e}")
        return config

    def _record_interaction(self, key: str, value: Any):
        str_value = str(value)
        self.recorded_interactions[key] = str_value
        self.input_dict[key] = str_value

        try:
            # Use Path's write method for atomic writing
            current_content = self.interactive_log.read_text() if self.interactive_log.exists() else ""
            updated_content = current_content + f"{key} = {str_value}\n"
            self.interactive_log.write_text(updated_content)

            print(f"Successfully logged: {key} = {str_value}")
        except Exception as e:
            print(f"Logging error details: {e}")
            logger.error(f"Error recording interaction: {e}")

    def prompt(self,
            key: str,
            message: str,
            choices: Optional[list] = None,
            allow_empty: bool = False,
            input_type: type = str,
            default: Optional[Any] = None) -> Any:
        
        if key in self.input_dict:
            value = self.input_dict[key]
            try:
                typed_value = input_type(value)
                if choices and typed_value not in choices:
                    raise ValueError
                return typed_value
            except ValueError:
                logger.warning(f"Invalid predefined input for {key}: {value}")

        while True:
            try:
                import click
                if choices:
                    value = click.prompt(message, type=click.Choice(choices),
                                        show_choices=True,
                                        default=default)
                else:
                    value = click.prompt(message, type=input_type,
                                        default=default)
                self._record_interaction(key, value)
                return value
                
            except ImportError:
                # If we have a default value, show it in the prompt
                if default is not None:
                    prompt_text = f"{message} [default: {default}]: "
                else:
                    prompt_text = message if message.endswith(': ') else f"{message}: "

                raw_value = input(prompt_text).strip()

                # Handle empty input (just pressing Enter)
                if not raw_value:
                    if default is not None:
                        self._record_interaction(key, default)
                        return default
                    elif allow_empty:
                        self._record_interaction(key, '')
                        return None
                    else:
                        print("Input cannot be empty. Please try again.")
                        continue

                # Try to convert non-empty input to the required type
                try:
                    value = input_type(raw_value)
                except ValueError:
                    print(f"Please enter a valid {input_type.__name__}")
                    continue

                if choices and value not in choices:
                    print(f"Invalid choice. Please choose from {choices}")
                    continue
                    
                self._record_interaction(key, value)
                return value

    def yes_no_prompt(self, key: str, message: str, default: Optional[bool] = None) -> bool:
        if key in self.input_dict:
            value = str(self.input_dict[key]).lower()
            result = value in ['yes', 'y', 'true', '1']
            return result

        while True:
            try:
                import click
                result = click.confirm(message, default=default)
                self._record_interaction(key, 'yes' if result else 'no')
                return result
            except ImportError:
                response = input(f"{message} (yes/no) ").lower().strip()
                if response in ['yes', 'y']:
                    self._record_interaction(key, 'yes')
                    return True
                elif response in ['no', 'n']:
                    self._record_interaction(key, 'no')
                    return False
                print("Please respond with 'yes' or 'no'")

    def get_recorded_interactions(self) -> Dict[str, str]:
        return self.recorded_interactions.copy()

    def get_input_dict(self) -> Dict[str, str]:
        return self.input_dict.copy()

    def __del__(self):
        if hasattr(self, '_log_file'):
            self._log_file.close()
