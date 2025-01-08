# Standard library imports
import base64
import importlib.resources
import importlib.util
import os
import queue
import subprocess
import sys
import threading
import time
from io import StringIO
from pathlib import Path
from typing import Dict, List, Optional

# Third-party imports
try:
    import pymol
    from pymol import cmd
    PYMOL_AVAILABLE = True
except ImportError:
    PYMOL_AVAILABLE = False

from Bio import PDB
from PIL import Image, ImageTk

# GUI-related imports
import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext, ttk
from tkinterweb import HtmlFrame

# Import the specific functions we need
from biodc.core.prep_modules.initialization import select_pdb, verify_programs
from biodc.core.preparation import PreparedStructure, generate_cpin as prep_generate_cpin
from biodc.core.prep_modules.select_disulfides import DisulfideFinder, select_disulfides
from biodc.core.prep_modules.select_mutate import select_mutate, StructureMutator
from biodc.core.prep_modules.ligand_detection import ligand_detection
from biodc.core.prep_modules.select_ph_active_sites import PHActiveSites, LigandDetector, ResidueSelector
from biodc.core.prep_modules.create_res_indexing import ResidueIndexing
from biodc.core.prep_modules.process_residues import PDBProcessor, HemeEnvironment
from biodc.core.prep_modules import generate_tleap
from biodc.core.prep_modules.struct_relax import struct_relax

def get_biodc_forcefield_dir():
    """
    Locate the forcefields directory within the biodc package.
    
    Returns:
        Path to the forcefields directory
    """
    try:
        # Try using importlib.resources first
        biodc_path = Path(importlib.resources.files('biodc'))
        forcefield_dir = biodc_path / 'data' / 'forcefield'
        
        # Verify the directory exists
        if forcefield_dir.is_dir():
            return forcefield_dir
        
        # Fallback method using package location
        spec = importlib.util.find_spec('biodc')
        if spec and spec.submodule_search_locations:
            pkg_path = Path(spec.submodule_search_locations[0])
            forcefield_dir = pkg_path / 'data' / 'forcefields'
            
            if forcefield_dir.is_dir():
                return forcefield_dir
        
        # Absolute fallback
        raise FileNotFoundError("Could not locate biodc forcefields directory")
    
    except Exception as e:
        print(f"Error locating forcefields directory: {e}")
        raise

class ConsoleRedirector(StringIO):
    def __init__(self, text_widget):
        super().__init__()
        self.text_widget = text_widget
        self.queue = queue.Queue()
        self.update_thread = threading.Thread(target=self._update_widget, daemon=True)
        self.update_thread.start()

    def write(self, string):
        self.queue.put(string)

    def flush(self):
        pass

    def _update_widget(self):
        while True:
            try:
                while True:
                    string = self.queue.get_nowait()
                    self.text_widget.configure(state='normal')
                    self.text_widget.insert('end', string)
                    self.text_widget.see('end')
                    self.text_widget.configure(state='disabled')
                    self.queue.task_done()
            except queue.Empty:
                self.text_widget.after(100, self._update_widget)
                break

DARK_BLUE = '#2E5CA6'    # Very dark navy
DARK_GREEN = '#156B3F'   # Very dark forest
DARK_RED = '#A31212'     # Very dark red

class ColorButton(ttk.Button):
    def __init__(self, master=None, **kwargs):
        super().__init__(master, **kwargs)
        self.configure(style=kwargs.get('style', 'Custom.TButton'))
        
        # Set up hover behavior
        self.bind('<Enter>', self._on_enter)
        self.bind('<Leave>', self._on_leave)
        self.bind('<Map>', self._on_map)
    
    def _on_enter(self, event):
        self.state(['active'])
        
    def _on_leave(self, event):
        self.state(['!active'])
        
    def _on_map(self, event):
        # Force color refresh
        style = ttk.Style()
        style.configure(self['style'], background=self._get_color())
        
    def _get_color(self):
        return self._style_colors.get(self['style'], '#1a4b8c')
        
    def _set_button_color(self, color):
        """Set button color and active state color"""
        style = ttk.Style()
        style.configure(self['style'], background=color)
        style.map(self['style'], background=[('active', color)])
        
    _style_colors = {
        'Custom.TButton': DARK_BLUE,
        'Nav.TButton': DARK_GREEN,
        'Action.TButton': DARK_RED
    }

class CustomButton(ColorButton):
    def __init__(self, master=None, **kwargs):
        super().__init__(master, style='Custom.TButton', **kwargs)
        if sys.platform == 'darwin':
            self.after(10, lambda: self._set_button_color(DARK_BLUE))

class NavButton(ColorButton):
    def __init__(self, master=None, **kwargs):
        super().__init__(master, style='Nav.TButton', **kwargs)
        if sys.platform == 'darwin':
            self.after(10, lambda: self._set_button_color(DARK_GREEN))

class ActionButton(ColorButton):
    def __init__(self, master=None, **kwargs):
        super().__init__(master, style='Action.TButton', **kwargs)
        if sys.platform == 'darwin':
            self.after(10, lambda: self._set_button_color(DARK_RED))

class InitializationFrame(ttk.LabelFrame):

    def __init__(self, parent, console):
        super().__init__(parent, text="Initialize Structure", padding=10)
        self.console = console
        self.input_dict = {}
        self.launch_dir = Path.cwd()
        self.current_pdb = None
        self.pdb_downloader = PDB.PDBList()
        self.create_widgets()

    def create_widgets(self):
        # Program verification section (keep existing code)
        verify_frame = ttk.LabelFrame(self, text="Program Verification", padding=5)
        verify_frame.pack(fill='x', pady=5)
        
        # Add an informative label for program requirements
        requirements_label = ttk.Label(
            verify_frame, 
            text="Required programs: tleap, cpptraj, and either sander or pmemd (from Amber/AmberTools). "
                "\n Press the 'Verify' button to make sure they are properly installed. "
                "If needed, they can be downloaded free of charge for academcis at https://ambermd.org/GetAmber.php",
            wraplength=600
        )
        requirements_label.pack(fill='x', pady=(0,5))

        self.verify_status = {}
        for program in ['tleap', 'cpptraj', 'sander', 'pmemd']:
            var = tk.StringVar(value="⚪")  # Unchecked
            self.verify_status[program] = var
            ttk.Label(verify_frame, text=f"{program}:").pack(side='left', padx=5)
            ttk.Label(verify_frame, textvariable=var).pack(side='left', padx=2)

        CustomButton(verify_frame, text="Verify", 
                command=self.verify_programs).pack(side='right', padx=5)

        # PDB selection section with both local and download options
        pdb_frame = ttk.LabelFrame(self, text="PDB Selection", padding=5)
        pdb_frame.pack(fill='x', pady=10)

        # Local PDB file selection
        local_frame = ttk.Frame(pdb_frame)
        local_frame.pack(fill='x', pady=5)
        
        ttk.Label(local_frame, text="Select Local PDB File:").pack(anchor='w')
        
        self.pdb_path = tk.StringVar()
        ttk.Entry(local_frame, textvariable=self.pdb_path, 
                width=30, style='Beige.TEntry').pack(side='left', padx=5, pady=5)
        ttk.Button(local_frame, text="Browse", 
                command=self.browse_pdb).pack(side='left', padx=5)

        # PDB download section
        download_frame = ttk.Frame(pdb_frame)
        download_frame.pack(fill='x', pady=5)
        
        ttk.Label(download_frame, text="Download from RCSB:").pack(anchor='w')
        
        download_entry_frame = ttk.Frame(download_frame)
        download_entry_frame.pack(fill='x')
        
        self.pdb_id = tk.StringVar()
        ttk.Entry(download_entry_frame, textvariable=self.pdb_id, 
                width=10, style='Beige.TEntry').pack(side='left', padx=5, pady=5)
        ttk.Button(download_entry_frame, text="Download", 
                command=self.download_pdb).pack(side='left', padx=5)
                
        # PDB info section
        self.pdb_info = ttk.Label(self, text="", justify='left')
        self.pdb_info.pack(fill='x', pady=5)

    def download_pdb(self):
        """Download PDB from RCSB."""
        pdb_id = self.pdb_id.get().strip().upper()
        if not pdb_id:
            messagebox.showerror("Error", "Please enter a PDB ID")
            return
            
        try:
            print(f"\nDownloading PDB {pdb_id}...")
            
            # Create SPR directory if it doesn't exist
            spr_dir = self.launch_dir / "SPR"
            spr_dir.mkdir(exist_ok=True)
            
            # Download the PDB file to the launch directory first
            pdb_file = self.pdb_downloader.retrieve_pdb_file(
                pdb_id,
                pdir=str(self.launch_dir),  # Download to launch dir first
                file_format="pdb"
            )
            
            if pdb_file and Path(pdb_file).exists():
                print(f"Successfully downloaded {pdb_id}")
                
                # Update input dictionary with PDB info
                self.input_dict['OriginalPDB'] = pdb_id
                
                # Process the downloaded PDB through select_pdb
                try:
                    result = select_pdb(self.launch_dir, self.input_dict)
                    if result:
                        print(f"\nPDB initialization successful: {result}")
                        # Update current PDB to the processed version
                        if (self.launch_dir / "SPR" / f"{result}.pdb").exists():
                            self.current_pdb = str(self.launch_dir / "SPR" / f"{result}.pdb")
                            self.pdb_path.set(self.current_pdb)
                        
                        self.pdb_info.config(
                            text=f"Downloaded and processed PDB: {pdb_id}\nFile ready for use."
                        )
                        return result
                        
                except Exception as e:
                    raise Exception(f"Error processing downloaded PDB: {str(e)}")
                    
            else:
                raise Exception(f"Failed to download PDB {pdb_id}")
                
        except Exception as e:
            print(f"\nError downloading/processing PDB: {str(e)}")
            messagebox.showerror("Error", f"Failed to download/process PDB: {str(e)}")

    def verify_programs(self):
        """Verify required programs using verify_programs()."""
        print("\nVerifying required programs...")
                
        try:
            # Existing verification code follows...
            missing_programs = []
            program_checks = {
                'tleap': ['tleap', '-h'],
                'cpptraj': ['cpptraj', '-h'],
                'sander': ['sander', '-h'],
                'pmemd': ['pmemd', '-h']
            }

            for program, cmd in program_checks.items():
                try:
                    subprocess.run(
                        cmd, 
                        stdout=subprocess.PIPE, 
                        stderr=subprocess.PIPE, 
                        timeout=5
                    )
                    self.verify_status[program].set("✓")
                except Exception as e:
                    print(f"Error checking {program}: {e}")
                    missing_programs.append(program)
                    self.verify_status[program].set("✗")

            if missing_programs:
                print(f"\nMissing programs: {', '.join(missing_programs)}")
                return False
            else:
                print("\nAll required programs verified successfully!")
                return True
        except Exception as e:
            print(f"Error during program verification: {str(e)}")
            for program in self.verify_status:
                self.verify_status[program].set("?")  # Error
            return False
    
    def browse_pdb(self):
        """Handle PDB file selection."""
        filename = filedialog.askopenfilename(
            title="Select PDB file",
            filetypes=[("PDB files", "*.pdb"), ("All files", "*.*")]
        )
        if filename:
            self.pdb_path.set(filename)
            print(f"\nSelected PDB file: {filename}")
            
            # Update input dictionary with PDB info
            self.input_dict['OriginalPDB'] = Path(filename).stem
            
            try:
                # Call select_pdb function directly
                result = select_pdb(self.launch_dir, self.input_dict)
                if result:
                    print(f"PDB initialization successful: {result}")
                    # Update current PDB to the renumbered version if it exists
                    if (self.launch_dir / "SPR" / f"{result}.pdb").exists():
                        self.current_pdb = str(self.launch_dir / "SPR" / f"{result}.pdb")
                    else:
                        self.current_pdb = filename
                    
                    self.pdb_info.config(
                        text=f"Initialized PDB: {result}\nFile ready for processing."
                    )
                    return result
            except Exception as e:
                print(f"Error during PDB initialization: {str(e)}")
                self.pdb_info.config(
                    text=f"Error during PDB initialization: {str(e)}"
                )
                messagebox.showerror("Error", f"Failed to initialize PDB: {str(e)}")

    def validate(self) -> bool:
        """Validate initialization step."""
        if not self.pdb_path.get():
            messagebox.showerror("Error", "Please select a PDB file.")
            return False
            
        if not Path(self.pdb_path.get()).exists():
            messagebox.showerror("Error", "Selected PDB file does not exist.")
            return False
            
        # Verify programs if not already verified
        if any(var.get() == "⚪" for var in self.verify_status.values()):
            if not self.verify_programs():
                return False
                
        return True

    def get_data(self) -> Dict:
        """Get collected data from this step."""
        # Create SPR directory if it doesn't exist
        spr_dir = self.launch_dir / "SPR"
        spr_dir.mkdir(exist_ok=True)
        
        # Use current_pdb if available, otherwise use original pdb_path
        pdb_path = self.current_pdb if self.current_pdb else self.pdb_path.get()
        
        # If this isn't already a file in SPR, copy it to SPR
        if Path(pdb_path).parent != spr_dir:
            import shutil
            dest_path = str(spr_dir / Path(pdb_path).name)
            shutil.copy2(pdb_path, dest_path)
            pdb_path = dest_path
            
        # Change working directory to SPR
        os.chdir(spr_dir)
        
        return {
            'pdb_file': pdb_path,
            'input_dict': self.input_dict
        }

    def update_pdb_info(self, pdb_path):
        """Update PDB info after mutation."""
        self.current_pdb = pdb_path  # Update current PDB tracking
        self.pdb_path.set(pdb_path)
        self.input_dict['OriginalPDB'] = Path(pdb_path).stem
        self.pdb_info.config(
            text=f"Current PDB: {Path(pdb_path).name}\nFile ready for processing."
        )

class SPRFrame(ttk.Frame):
    def __init__(self, parent):
        super().__init__(parent)
        self.parent = parent

        style = ttk.Style()
        style.theme_use('default')

        style.configure('Beige.TEntry', fieldbackground='beige', background='beige')

        # Set dark wood background color
        wood_color = '#4A3728'  # Dark wood brown

        # Configure styles for all ttk widgets
        style.configure('TFrame', background=wood_color)
        style.configure('TLabelframe', background=wood_color)
        style.configure('TLabelframe.Label', background=wood_color, foreground='white')
        style.configure('TPanedwindow', background=wood_color)

        # Much darker button colors for Mac
        dark_blue = '#2E5CA6'    # Very dark navy
        dark_green = '#156B3F'   # Very dark forest
        dark_red = '#A31212'     # Very dark red

        # Configure button styles with explicit element layouts
        for btn_style, bg_color, active_color in [
            ('Custom.TButton', dark_blue, '#030F23'),
            ('Nav.TButton', dark_green, '#021509'),
            ('Action.TButton', dark_red, '#A31212')
        ]:
            # Base style configuration
            style.configure(btn_style,
                background=bg_color,
                foreground='white',
                padding=10,
                font=('Helvetica', 12, 'bold'),
                relief='solid',
                borderwidth=1
            )
            
            # Map configuration for active state
            style.map(btn_style,
                foreground=[('active', 'white')],
                background=[('active', active_color)]
            )
            
            # Explicit button element layout
            style.layout(btn_style, [
                ('Button.border', {'children': [
                    ('Button.padding', {'children': [
                        ('Button.label', {'sticky': 'nswe'})
                    ], 'sticky': 'nswe'})
                ], 'sticky': 'nswe', 'border': '1'})
            ])

        # Configure other widgets to match the dark theme
        style.configure('TLabel', background=wood_color, foreground='white')
        style.configure('Treeview', 
            background=wood_color,
            fieldbackground=wood_color,
            foreground='white'
        )
        
        style.configure('Beige.TEntry', 
                        fieldbackground='beige', 
                        background='beige', 
                        insertwidth=2,  # Cursor width
                        insertcolor='black',  # Cursor color
                        insertbackground='black'  # Ensures cursor is visible
        )

        # If you want the root window to be wood colored, do this:
#       self.parent.configure(bg=wood_color)

        # Create split pane layout
        self.paned = ttk.PanedWindow(self, orient='horizontal')
        self.paned.pack(expand=True, fill='both')
        
        # Left frame for steps
        self.main_container = ttk.Frame(self.paned)
        self.paned.add(self.main_container, weight=7)
        
        # Right frame for console output
        self.console_frame = ttk.LabelFrame(self.paned, text="Console Output")
        self.paned.add(self.console_frame, weight=3)
        
        # Create console output widget
        self.create_console()
        
        # Create frames for each step
        self.create_all_frames()
        
        # Create navigation
        self.create_navigation()
        
        # Show first step
        self.show_step(0)
        
        # Redirect stdout and stderr to our console
        sys.stdout = ConsoleRedirector(self.console)
        sys.stderr = ConsoleRedirector(self.console)
        
        print("SPR GUI initialized")
        print("Please verify programs and select a PDB file to begin...")

    def create_console(self):
        """Create console output widget."""
        console_container = ttk.Frame(self.console_frame)
        console_container.pack(fill='both', expand=True)
        
        # Create text widget with scrollbar
        self.console = scrolledtext.ScrolledText(
            console_container,
            wrap=tk.WORD,
            height=10,
            background='black',
            foreground='light green',
            font=('Courier', 12)
        )
        self.console.pack(expand=True, fill='both')
        self.console.configure(state='disabled')

        # Add console controls
        control_frame = ttk.Frame(self.console_frame)
        control_frame.pack(fill='x', padx=5, pady=2)
        
        ActionButton(
            control_frame, 
            text="Clear Console",
            command=self.clear_console
        ).pack(side='left', padx=2)
        
        NavButton(
            control_frame,
            text="Save Output",
            command=self.save_console
        ).pack(side='left', padx=2)

    def clear_console(self):
        """Clear the console output."""
        self.console.configure(state='normal')
        self.console.delete(1.0, tk.END)
        self.console.configure(state='disabled')
        print("Console cleared.")

    def save_console(self):
        """Save console output to a file."""
        filename = filedialog.asksaveasfilename(
            defaultextension=".log",
            filetypes=[("Log files", "*.log"), ("All files", "*.*")]
        )
        if filename:
            with open(filename, 'w') as f:
                f.write(self.console.get(1.0, tk.END))
            print(f"\nConsole output saved to: {filename}")

    def create_all_frames(self):
        """Create all step frames."""
        self.frames = {}
        
        # Initialize structure frame
        self.frames['Initialize'] = InitializationFrame(self.main_container, self.console)
        
        # Create Disulfide frame
        self.frames['Disulfide'] = DisulfideFrame(self.main_container, self.console)

        # Create Mutations frame
        self.frames['Mutations'] = MutationsFrame(self.main_container, self.console)

        # Create pH-active site selection frame
        self.frames['pH'] = PHFrame(self.main_container, self.console)

        # Create Processing frame
        self.frames['Processing'] = ProcessingFrame(self.main_container, self.console)

        # Create visualization frame
        self.frames['Visualization'] = VisualizationFrame(self.main_container, self.console)

    def create_navigation(self):
        """Create navigation buttons."""
        nav_frame = ttk.Frame(self.main_container)
        nav_frame.pack(fill='x', pady=10)
        
        self.prev_button = ActionButton(nav_frame, text="Previous", command=self.prev_step)
        self.prev_button.pack(side='left', padx=5)
        
        self.next_button = NavButton(nav_frame, text="Next", command=self.next_step)
        self.next_button.pack(side='right', padx=5)
        
        self.progress_label = ttk.Label(nav_frame, text=f"Step 1/{len(self.frames)}")
        self.progress_label.pack()

    def show_step(self, step_index):
        """Show the specified step and update navigation."""
        steps = list(self.frames.keys())
        
        # Hide all frames
        for frame in self.frames.values():
            frame.pack_forget()
        
        # Show current frame
        self.frames[steps[step_index]].pack(fill='both', expand=True)
        
        # Update navigation
        self.prev_button['state'] = 'normal' if step_index > 0 else 'disabled'
        self.next_button['state'] = 'normal' if step_index < len(steps) - 1 else 'disabled'
        self.progress_label['text'] = f"Step {step_index + 1}/{len(steps)}"
        
        self.current_step = step_index
        print(f"\nMoving to step: {steps[step_index]}")

    def next_step(self):
        """Move to next step with validation."""
        current_frame = self.frames[list(self.frames.keys())[self.current_step]]
        if hasattr(current_frame, 'validate') and not current_frame.validate():
            return
            
        self.show_step(min(self.current_step + 1, len(self.frames) - 1))

    def prev_step(self):
        """Move to previous step."""
        self.show_step(max(self.current_step - 1, 0))

class DisulfideFrame(ttk.LabelFrame):
    def __init__(self, parent, console):
        super().__init__(parent, text="Disulfide Bond Selection", padding=10)
        self.console = console
        self.disulf_res_list = " "
        self.disulf_pair_id = []
        self.potential_pairs = []
        self.create_widgets()
        
    def create_widgets(self):
        # Auto-detection section
        detect_frame = ttk.LabelFrame(self, text="Automatic Detection", padding=5)
        detect_frame.pack(fill='x', pady=5)
        
        NavButton(detect_frame, text="Detect Potential Disulfide Bonds", 
                  command=self.detect_disulfides).pack(pady=5)
        
        # Display potential pairs
        self.pairs_frame = ttk.LabelFrame(self, text="Potential Disulfide Pairs", padding=5)
        self.pairs_frame.pack(fill='x', pady=5)
        
        # Create treeview for pairs
        columns = ("Pair", "Residue 1", "Residue 2", "Distance", "Status")
        self.pairs_tree = ttk.Treeview(self.pairs_frame, columns=columns, show="headings")
        
        for col in columns:
            self.pairs_tree.heading(col, text=col)
            self.pairs_tree.column(col, width=100)
            
        self.pairs_tree.pack(fill='x', pady=5)
        
        # Add buttons for pair selection
        btn_frame = ttk.Frame(self.pairs_frame)
        btn_frame.pack(fill='x', pady=5)
        
        NavButton(btn_frame, text="Accept Selected", 
                  command=self.accept_selected).pack(side='left', padx=5)
        ActionButton(btn_frame, text="Reject Selected", 
                  command=self.reject_selected).pack(side='left', padx=5)
        
        # Manual entry section
        manual_frame = ttk.LabelFrame(self, text="Manual Entry", padding=5)
        manual_frame.pack(fill='x', pady=5)
        
        # List of manual pairs
        self.manual_tree = ttk.Treeview(manual_frame, columns=("Res1", "Res2"), show="headings", height=4)
        self.manual_tree.heading("Res1", text="Residue 1")
        self.manual_tree.heading("Res2", text="Residue 2")
        self.manual_tree.pack(fill='x', pady=5)
        
        entry_frame = ttk.Frame(manual_frame)
        entry_frame.pack(fill='x', pady=5)
        
        ttk.Label(entry_frame, text="Residue 1:").pack(side='left', padx=5)
        self.res1_var = tk.StringVar()
        ttk.Entry(entry_frame, textvariable=self.res1_var, width=10, style='Beige.TEntry').pack(side='left', padx=5)
        
        ttk.Label(entry_frame, text="Residue 2:").pack(side='left', padx=5)
        self.res2_var = tk.StringVar()
        ttk.Entry(entry_frame, textvariable=self.res2_var, width=10, style='Beige.TEntry').pack(side='left', padx=5)
        
        btn_frame = ttk.Frame(manual_frame)
        btn_frame.pack(fill='x', pady=5)
        
        NavButton(btn_frame, text="Add Pair", 
                  command=self.add_manual_pair).pack(side='left', padx=5)
        ActionButton(btn_frame, text="Remove Selected", 
                  command=self.remove_manual_pair).pack(side='left', padx=5)

    def detect_disulfides(self):
        """Detect potential disulfide bonds in the structure."""
        try:
            # Get the SPRFrame instance by traversing up the widget hierarchy
            spr_frame = self.master.master.master
            pdb_path = spr_frame.frames['Initialize'].get_data()['pdb_file']  # This will now get the current/renumbered PDB
            finder = DisulfideFinder(pdb_path)
            self.potential_pairs = finder.find_potential_disulfides()
            
            # Clear existing items
            for item in self.pairs_tree.get_children():
                self.pairs_tree.delete(item)
            
            # Add detected pairs to treeview
            for i, (res1, res2, distance) in enumerate(self.potential_pairs, 1):
                self.pairs_tree.insert('', 'end', values=(
                    f"Pair {i}",
                    f"CYS {res1}",
                    f"CYS {res2}",
                    f"{distance:.2f} Å",
                    "Pending"
                ))
                
            print(f"\nDetected {len(self.potential_pairs)} potential disulfide pairs")
            
        except Exception as e:
            print(f"\nError detecting disulfide bonds: {str(e)}")
            messagebox.showerror("Error", f"Failed to detect disulfide bonds: {str(e)}")

    def accept_selected(self):
        """Accept selected pairs from the treeview."""
        selected = self.pairs_tree.selection()
        for item in selected:
            values = self.pairs_tree.item(item)['values']
            pair_idx = int(values[0].split()[1]) - 1
            res1, res2, _ = self.potential_pairs[pair_idx]
            
            # Add to disulfide lists if not already present
            if [res1, res2] not in self.disulf_pair_id:
                self.disulf_pair_id.append([res1, res2])
                self.disulf_res_list += f"{res1} {res2} "
                
            # Update status in treeview
            self.pairs_tree.set(item, "Status", "Accepted")
            print(f"\nAccepted disulfide pair: CYS {res1} - CYS {res2}")
            
    def reject_selected(self):
        """Reject selected pairs from the treeview."""
        selected = self.pairs_tree.selection()
        for item in selected:
            values = self.pairs_tree.item(item)['values']
            pair_idx = int(values[0].split()[1]) - 1
            res1, res2, _ = self.potential_pairs[pair_idx]
            
            # Remove from disulfide lists if present
            if [res1, res2] in self.disulf_pair_id:
                self.disulf_pair_id.remove([res1, res2])
                self.disulf_res_list = self.disulf_res_list.replace(f"{res1} {res2} ", " ")
                
            # Update status in treeview
            self.pairs_tree.set(item, "Status", "Rejected")
            print(f"\nRejected disulfide pair: CYS {res1} - CYS {res2}")
            
    def validate_cysteine(self, residue_number):
        """Validate if a residue number corresponds to a Cysteine."""
        try:
            spr_frame = self.master.master.master
            pdb_path = spr_frame.frames['Initialize'].get_data()['pdb_file']  # Get current PDB
            structure = PDB.PDBParser(QUIET=True).get_structure('protein', pdb_path)
            
            # Check all models, chains, and residues
            for model in structure:
                for chain in model:
                    for residue in chain:
                        if residue.get_id()[1] == residue_number:
                            return residue.get_resname() == "CYS"
            return False
        except Exception as e:
            print(f"\nError validating residue: {str(e)}")
            return False

    def add_manual_pair(self):
        """Add manually entered disulfide pair with validation."""
        try:
            res1 = int(self.res1_var.get())
            res2 = int(self.res2_var.get())
            
            # Validate both residues are Cysteines
            if not self.validate_cysteine(res1):
                messagebox.showerror("Error", f"Residue {res1} is not a Cysteine")
                return
            if not self.validate_cysteine(res2):
                messagebox.showerror("Error", f"Residue {res2} is not a Cysteine")
                return
            
            if [res1, res2] not in self.disulf_pair_id:
                self.disulf_pair_id.append([res1, res2])
                self.disulf_res_list += f"{res1} {res2} "
                print(f"\nManually added disulfide pair: CYS {res1} - CYS {res2}")
                
                # Add to treeview
                self.manual_tree.insert('', 'end', values=(f"CYS {res1}", f"CYS {res2}"))
                
                # Clear entry fields
                self.res1_var.set("")
                self.res2_var.set("")
            else:
                print("\nThis pair is already in the list")
                
        except ValueError:
            print("\nPlease enter valid residue numbers")
            messagebox.showerror("Error", "Please enter valid residue numbers")
            
    def remove_manual_pair(self):
        """Remove selected manual pair."""
        selected = self.manual_tree.selection()
        for item in selected:
            values = self.manual_tree.item(item)['values']
            res1 = int(values[0].split()[1])
            res2 = int(values[1].split()[1])
            
            if [res1, res2] in self.disulf_pair_id:
                self.disulf_pair_id.remove([res1, res2])
                self.disulf_res_list = self.disulf_res_list.replace(f"{res1} {res2} ", " ")
                self.manual_tree.delete(item)
                print(f"\nRemoved disulfide pair: CYS {res1} - CYS {res2}")
            
    def validate(self) -> bool:
        """Validate disulfide selection step."""
        # Write disulfide definitions if pairs exist
        if self.disulf_pair_id:
            spr_frame = self.master.master.master
            launch_dir = spr_frame.frames['Initialize'].launch_dir
            with open(launch_dir / 'DisulfideDefinitions.txt', 'w') as f:
                f.write(self.disulf_res_list)
        return True

    def get_data(self) -> Dict:
        """Get collected data from this step."""
        return {
            'disulf_res_list': self.disulf_res_list,
            'disulf_pair_id': self.disulf_pair_id
        }

class MutationsFrame(ttk.LabelFrame):
    def __init__(self, parent, console):
        super().__init__(parent, text="Mutations", padding=10)
        self.console = console
        self.mutations = []
        self.residue_data = []
        self.create_widgets()
        
    def create_widgets(self):
        # Residue list
        list_frame = ttk.LabelFrame(self, text="Structure Residues", padding=5)
        list_frame.pack(fill='both', expand=True, pady=5)
        
        # Create treeview with scrollbar
        tree_frame = ttk.Frame(list_frame)
        tree_frame.pack(fill='both', expand=True)
        
        self.tree_scroll = ttk.Scrollbar(tree_frame)
        self.tree_scroll.pack(side='right', fill='y')
        
        columns = ("Chain", "ResID", "Original", "New", "Status")
        self.residue_tree = ttk.Treeview(tree_frame, columns=columns, show="headings", 
                                       height=15, yscrollcommand=self.tree_scroll.set)
        
        for col in columns:
            self.residue_tree.heading(col, text=col)
            self.residue_tree.column(col, width=80)
            
        self.residue_tree.pack(side='left', fill='both', expand=True)
        self.tree_scroll.config(command=self.residue_tree.yview)
        
        # Bind double-click to edit
        self.residue_tree.bind('<Double-1>', self.on_double_click)
        
        # Control buttons
        btn_frame = ttk.Frame(self)
        btn_frame.pack(fill='x', pady=5)
        
        CustomButton(btn_frame, text="Load Structure", 
                  command=self.load_structure).pack(side='left', padx=5)
        NavButton(btn_frame, text="Apply Mutations", 
                  command=self.apply_mutations).pack(side='left', padx=5)
        ActionButton(btn_frame, text="Reset Selected", 
                  command=self.reset_mutations).pack(side='left', padx=5)

    def load_structure(self):
        try:
            spr_frame = self.master.master.master
            init_data = spr_frame.frames['Initialize'].get_data()
            pdb_path = init_data['pdb_file']
            
            print(f"\nAttempting to load structure from: {pdb_path}")
            structure = PDB.PDBParser(QUIET=True).get_structure('protein', pdb_path)
            
            # Clear existing items
            for item in self.residue_tree.get_children():
                self.residue_tree.delete(item)
            
            self.residue_data = []  # Reset residue data
            
            # Add all residues to treeview
            for model in structure:
                for chain in model:
                    for residue in chain:
                        res_id = residue.get_id()
                        
                        # Include all residues, not just standard ones
                        res_name = residue.get_resname()
                        chain_id = chain.get_id()
                        res_num = str(res_id[1])
                        
                        self.residue_data.append({
                            'chain': chain_id,
                            'resid': res_num,
                            'original': res_name,
                            'new': res_name
                        })
                        
                        self.residue_tree.insert('', 'end', values=(
                            chain_id, res_num, res_name, res_name, "Original"
                        ))
            
            print(f"\nLoaded {len(self.residue_data)} residues from structure")
            
        except Exception as e:
            print(f"\nError loading structure: {str(e)}")
            messagebox.showerror("Error", f"Failed to load structure: {str(e)}")
      
    def on_double_click(self, event):
        item = self.residue_tree.selection()[0]
        column = self.residue_tree.identify_column(event.x)
        
        # Only allow editing the "New" column (index 3)
        if column == '#4':  # New residue column
            x, y, w, h = self.residue_tree.bbox(item, column)
            
            # Create entry widget
            entry = ttk.Entry(self.residue_tree, width=10, style='Beige.TEntry')
            entry.place(x=x, y=y, width=w, height=h)
            
            # Get current value
            value = self.residue_tree.item(item)['values'][3]
            entry.insert(0, value)
            entry.select_range(0, tk.END)
            entry.focus()
            
            def on_enter(event):
                new_value = entry.get().upper()
                if len(new_value) == 3:  # Validate three-letter code
                    values = list(self.residue_tree.item(item)['values'])
                    if new_value != values[2]:  # If different from original
                        values[3] = new_value  # Update new residue
                        values[4] = "Pending"  # Update status
                        self.residue_tree.item(item, values=values)
                    else:
                        values[3] = values[2]  # Reset to original
                        values[4] = "Original"  # Reset status
                        self.residue_tree.item(item, values=values)
                entry.destroy()
                
            def on_escape(event):
                entry.destroy()
                
            entry.bind('<Return>', on_enter)
            entry.bind('<Escape>', on_escape)
            entry.bind('<FocusOut>', lambda e: entry.destroy())
            
    def reset_mutations(self):
        selected = self.residue_tree.selection()
        for item in selected:
            values = list(self.residue_tree.item(item)['values'])
            values[3] = values[2]  # Reset new to original
            values[4] = "Original"  # Reset status
            self.residue_tree.item(item, values=values)

    def apply_mutations(self):
        try:
            mutations = []
            for item in self.residue_tree.get_children():
                values = self.residue_tree.item(item)['values']
                if values[4] == "Pending":
                    mutations.append((values[2], values[1], values[3]))  # orig, pos, new
            
            if not mutations:
                print("\nNo mutations to apply")
                return
                
            spr_frame = self.master.master.master
            init_data = spr_frame.frames['Initialize'].get_data()
            pdb_path = init_data['pdb_file']  # This will now be the current (possibly renumbered) PDB
            launch_dir = spr_frame.frames['Initialize'].launch_dir
            
            # Prepare input dictionary with mutation information
            input_dict = {
                "NumMut": len(mutations),
            }
            
            # Add mutation details to input dictionary
            for i, (orig, pos, new) in enumerate(mutations):
                input_dict[f"Mutation_{i+1}"] = f"{orig} {pos} {new}"
            
            # Get base PDB name without extension
            pdb_base = Path(pdb_path).stem
            
            # Apply mutations using the prepared input dictionary
            result_pdb = select_mutate(pdb_base, launch_dir, input_dict)
            
            # Update the PDB path in Initialize frame if mutations were applied
            if result_pdb == "mutated":
                new_pdb_path = str(launch_dir / "SPR" / "mutated.pdb")
                spr_frame.frames['Initialize'].update_pdb_info(new_pdb_path)
                print(f"\nUpdated PDB path to: {new_pdb_path}")
            
            print("\nMutations applied successfully")
            
            # Update status in tree
            for item in self.residue_tree.get_children():
                values = self.residue_tree.item(item)['values']
                if values[4] == "Pending":
                    values = list(values)
                    values[4] = "Applied"
                    self.residue_tree.item(item, values=values)
                        
        except Exception as e:
            print(f"\nError applying mutations: {str(e)}")
            messagebox.showerror("Error", f"Failed to apply mutations: {str(e)}")

    def validate(self) -> bool:
        return True
    
    def get_data(self) -> Dict:
        mutations = []
        for item in self.residue_tree.get_children():
            values = self.residue_tree.item(item)['values']
            if values[4] in ["Pending", "Applied"]:
                mutations.append((values[2], values[1], values[3]))
        return {
            'mutations': mutations
        }

class PHFrame(ttk.LabelFrame):
     
    def __init__(self, parent, console):
        super().__init__(parent, text="pH-Active Sites Selection", padding=10)
        self.console = console
        self.sites = PHActiveSites()
        self.create_widgets()

    def create_widgets(self):
        # Metal coordination info
        info_frame = ttk.LabelFrame(self, text="Metal Coordination Info", padding=5)
        info_frame.pack(fill='x', pady=5)
        self.coord_text = scrolledtext.ScrolledText(info_frame, height=5)
        self.coord_text.pack(fill='x', padx=5, pady=5)
        
        # Residue selection frame
        select_frame = ttk.LabelFrame(self, text="Residue Selection", padding=5)
        select_frame.pack(fill='both', expand=True, pady=5)
        
        # Create residue selection tables
        self.residue_tables = {}
        residue_types = ["ASP", "GLU", "HIS", "TYR", "LYS", "PRN"]
        
        for i, res_type in enumerate(residue_types):
            frame = ttk.Frame(select_frame)
            frame.grid(row=i//3, column=i%3, padx=5, pady=5, sticky='nsew')
            
            ttk.Label(frame, text=f"{res_type} Residues").pack()
            
            # Create treeview
            tree = ttk.Treeview(frame, columns=("ID", "Status"), show="headings", height=6)
            tree.heading("ID", text="ID")
            tree.heading("Status", text="Status")
            tree.column("ID", width=50)
            tree.column("Status", width=80)
            tree.pack(fill='x')
            
            self.residue_tables[res_type] = tree
            
        # Configure grid
        select_frame.columnconfigure(0, weight=1)
        select_frame.columnconfigure(1, weight=1)
        select_frame.columnconfigure(2, weight=1)
        
        # Control buttons
        btn_frame = ttk.Frame(self)
        btn_frame.pack(fill='x', pady=5)
        
        NavButton(btn_frame, text="Analyze Structure", 
                  command=self.analyze_structure).pack(side='left', padx=5)
                       
        CustomButton(btn_frame, text="Select/Deselect All", 
                  command=self.toggle_all).pack(side='left', padx=5)
                  
        # Bind double-click for all trees
        for tree in self.residue_tables.values():
            tree.bind('<Double-1>', self.toggle_selection)

    def toggle_selection(self, event):
        """Toggle residue selection status on double-click."""
        tree = event.widget
        item = tree.selection()[0]
        values = list(tree.item(item)['values'])
        values[1] = "Selected" if values[1] == "Available" else "Available"
        tree.item(item, values=values)

    def toggle_all(self):
        """Toggle all residues between Selected and Available."""
        # Check first item's status to determine new status
        first_tree = list(self.residue_tables.values())[0]
        if first_tree.get_children():
            first_status = first_tree.item(first_tree.get_children()[0])['values'][1]
            new_status = "Available" if first_status == "Selected" else "Selected"
            
            # Apply to all items
            for tree in self.residue_tables.values():
                for item in tree.get_children():
                    values = list(tree.item(item)['values'])
                    values[1] = new_status
                    tree.item(item, values=values)
    
    def analyze_structure(self):
        try:
            spr_frame = self.master.master.master
            init_data = spr_frame.frames['Initialize'].get_data()
            launch_dir = spr_frame.frames['Initialize'].launch_dir
            pdb_path = init_data['pdb_file']  # This will get the current/renumbered PDB
            pdb_name = Path(pdb_path).stem
            
            # Clear existing entries
            for tree in self.residue_tables.values():
                for item in tree.get_children():
                    tree.delete(item)
            
            # Analyze metal coordination
            detector = LigandDetector(f"{pdb_name}.pdb")
            coordinated = detector.find_metal_coordinations()
            
            # Update coordination info text
            self.coord_text.delete('1.0', tk.END)
            detector.print_coordination_info()
            excluded = detector.get_excluded_residues()
            
            if any(excluded.values()):
                self.coord_text.insert(tk.END, "\nExcluded metal-coordinated residues:\n")
                for res_type, ids in excluded.items():
                    if ids:
                        self.coord_text.insert(tk.END, f"{res_type}: {', '.join(map(str, sorted(ids)))}\n")
            
            # Find titratable residues
            selector = ResidueSelector(pdb_name, launch_dir)
            
            for res_type in self.residue_tables:
                residlist = selector.find_titratable_residues(res_type)
                if res_type in excluded:
                    residlist = [r for r in residlist if int(r) not in excluded[res_type]]
                    
                tree = self.residue_tables[res_type]
                for res_id in residlist:
                    tree.insert('', 'end', values=(res_id, "Available"))
                    
                # Store selections
                if res_type == "ASP": self.sites.asp_ids = " ".join(residlist) if residlist else ""
                elif res_type == "GLU": self.sites.glu_ids = " ".join(residlist) if residlist else ""
                elif res_type == "HIS": self.sites.his_ids = " ".join(residlist) if residlist else ""
                elif res_type == "LYS": self.sites.lys_ids = " ".join(residlist) if residlist else ""
                elif res_type == "TYR": self.sites.tyr_ids = " ".join(residlist) if residlist else ""
                elif res_type == "PRN": self.sites.prn_ids = " ".join(residlist) if residlist else ""
                    
            print("\nStructure analysis complete")
            
        except Exception as e:
            print(f"\nError analyzing structure: {str(e)}")
            messagebox.showerror("Error", f"Failed to analyze structure: {str(e)}")

    def validate(self) -> bool:
        # Collect selected residues
        for res_type, tree in self.residue_tables.items():
            selected = []
            for item in tree.get_children():
                values = tree.item(item)['values']
                if values[1] == "Selected":
                    # Convert to string to ensure consistency
                    selected.append(str(values[0]))
                    
            selection = " ".join(selected)
            if res_type == "ASP": self.sites.asp_ids = selection
            elif res_type == "GLU": self.sites.glu_ids = selection
            elif res_type == "HIS": self.sites.his_ids = selection
            elif res_type == "LYS": self.sites.lys_ids = selection
            elif res_type == "TYR": self.sites.tyr_ids = selection
            elif res_type == "PRN": self.sites.prn_ids = selection
        
        return True
    
    def get_data(self) -> Dict:
        return {
            'asp_ids': self.sites.asp_ids,
            'glu_ids': self.sites.glu_ids,
            'his_ids': self.sites.his_ids,
            'lys_ids': self.sites.lys_ids,
            'tyr_ids': self.sites.tyr_ids,
            'prn_ids': self.sites.prn_ids
        }

class RedoxSelectionFrame(ttk.LabelFrame):
    """Frame for selecting redox states of hemes."""
    def __init__(self, parent):
        super().__init__(parent, text="Heme Redox States", padding=5)
        self.redox_states = {}  # Store selections: {heme_id: 'ox'/'red'}
        self.create_widgets()
        
    def create_widgets(self):
        self.tree = ttk.Treeview(self, columns=("Heme", "Type", "Ligation", "Redox"), 
                                show="headings", height=6)
        
        self.tree.heading("Heme", text="Heme ID")
        self.tree.heading("Type", text="Type")
        self.tree.heading("Ligation", text="Ligation")
        self.tree.heading("Redox", text="Redox")
        
        # Adjust column widths
        self.tree.column("Heme", width=80)
        self.tree.column("Type", width=80)
        self.tree.column("Ligation", width=80)
        self.tree.column("Redox", width=100)
        
        self.tree.pack(fill='x', pady=5)
        
        # Bind double-click to toggle redox state
        self.tree.bind('<Double-1>', self.toggle_redox_state)
        
        # Add help text
        ttk.Label(self, text="Double-click a row to toggle between oxidized and reduced states").pack()
        
    def load_hemes(self, indexing_file):
        """Load hemes from ResIndexing.txt."""
        self.redox_states.clear()
        for item in self.tree.get_children():
            self.tree.delete(item)
            
        try:
            with open(indexing_file) as f:
                for line in f:
                    parts = line.strip().split()
                    if not parts:
                        continue
                        
                    if len(parts) == 7:  # c-type: CysB CysC HisP Distal Heme c HX
                        heme_id = parts[4]
                        heme_type = "c-type"
                        ligation = parts[6]
                    elif len(parts) == 5:  # b-type: HisP Distal Heme b HX
                        heme_id = parts[2]
                        heme_type = "b-type"
                        ligation = parts[4]
                    else:
                        continue
                    
                    # Default to oxidized state
                    self.redox_states[heme_id] = 'ox'
                    self.tree.insert('', 'end', values=(heme_id, heme_type, ligation, "Oxidized"))
                    
        except Exception as e:
            print(f"\nError loading hemes from ResIndexing.txt: {str(e)}")

    def toggle_redox_state(self, event):
        """Toggle redox state of selected heme."""
        region = self.tree.identify("region", event.x, event.y)
        if region != "cell":
            return
            
        item = self.tree.identify('item', event.x, event.y)
        if not item:
            return
            
        values = self.tree.item(item)['values']
        if not values:
            return
            
        heme_id = str(values[0])  # Convert to string to match dictionary keys
        
        # Toggle state
        current = self.redox_states.get(heme_id, 'ox')
        new_state = 'red' if current == 'ox' else 'ox'
        self.redox_states[heme_id] = new_state
        
        # Update display
        values = list(values)
        values[3] = "Reduced" if new_state == 'red' else "Oxidized"
        self.tree.item(item, values=values)
        print(f"Toggled heme {heme_id} to {new_state}")  # Debug print        

    def get_redox_states(self):
        """Return list of redox states in order of hemes."""
        states = []
        for item in self.tree.get_children():
            values = self.tree.item(item)['values']
            heme_id = str(values[0])  # Convert to string to match dictionary keys
            state = 'R' if self.redox_states.get(heme_id, 'ox') == 'red' else 'O'
            states.append(state)
        print(f"Returning redox states: {states}")  # Debug print
        return states

class ProcessingFrame(ttk.LabelFrame):

    def __init__(self, parent, console):
        super().__init__(parent, text="Structure Processing", padding=10)
        self.console = console
        self.steps = [
            "Create Residue Indexing",
            "Process Residue Indexing",
            "Generate TLeap Input",
            "Reorder Structure",
            "Validate Titratable Sites",
            "Generate CPIN File",
            "Structure Relaxation"
        ]
        self.current_step = 0
        
        # Initialize ResIndex variables
        self.indexing_var = tk.StringVar(value='auto')
        
        # Initialize TLeap variables
        self.prefix_var = tk.StringVar()
        self.solvent_var = tk.StringVar(value='exp')
        self.box_type_var = tk.StringVar(value='rec')
        self.buffer_var = tk.StringVar(value='10.0')
        self.na_var = tk.StringVar(value='0')
        self.cl_var = tk.StringVar(value='0')
        
        self.create_widgets()

    def create_widgets(self):
        # Main horizontal container
        main_container = ttk.Frame(self)
        main_container.pack(expand=True, fill='both')
        
        # Left: Checklist
        checklist_frame = ttk.LabelFrame(main_container, text="Processing Steps", padding=5)
        checklist_frame.pack(side='left', fill='y', padx=5, pady=5)
        
        self.status_vars = []
        self.step_labels = []
        for i, step in enumerate(self.steps):
            var = tk.StringVar(value="⚪")  # Pending
            self.status_vars.append(var)
            
            frame = ttk.Frame(checklist_frame)
            frame.pack(fill='x', pady=2)
            
            ttk.Label(frame, textvariable=var).pack(side='left', padx=2)
            label = ttk.Label(frame, text=step)
            label.pack(side='left', padx=5)
            self.step_labels.append(label)
        
        # Right: Instructions and Controls
        center_frame = ttk.Frame(main_container)
        center_frame.pack(side='left', fill='both', expand=True, padx=5, pady=5)
        
        # Instructions box at top
        self.instructions_frame = ttk.LabelFrame(center_frame, text="Instructions", padding=5)
        self.instructions_frame.pack(fill='x', pady=5)
        
        self.instructions_text = ttk.Label(self.instructions_frame, wraplength=400, justify='left')
        self.instructions_text.pack(fill='x', padx=5, pady=5)
        
        # Controls container below instructions
        self.controls_container = ttk.Frame(center_frame)
        self.controls_container.pack(fill='both', expand=True, pady=5)
        
        # Create all the control frames but don't pack them yet
        self.create_control_frames(self.controls_container)
        
        # Navigation buttons at bottom
        nav_frame = ttk.Frame(center_frame)
        nav_frame.pack(fill='x', pady=5)
        
        self.start_button = NavButton(nav_frame, text="Start Processing", 
                                    command=self.start_current_step)
        self.start_button.pack(side='left', padx=5)
        
        self.next_button = NavButton(nav_frame, text="Next Step", 
                                    command=self.next_step, state='disabled')
        self.next_button.pack(side='left', padx=5)

        # Add the progress label
        self.progress_label = ttk.Label(nav_frame, text="Step 1/7")
        self.progress_label.pack(side='right', padx=5)

        # Add this at the very end of create_widgets:
        print("DEBUG: End of create_widgets, showing initial step")
        self.show_step(0)  # Force show initial step

    def create_control_frames(self, container):
        """Create all the step-specific control frames."""

        # 1. Create Residue Indexing controls
        self.resindex_frame = ttk.Frame(container)
        
        # Radio buttons for automatic/manual selection
        selection_frame = ttk.Frame(self.resindex_frame)
        selection_frame.pack(fill='x', pady=5)
        print("DEBUG: Created selection_frame with auto/manual controls")

        self.indexing_var = tk.StringVar(value='auto')
        ttk.Radiobutton(selection_frame, text="Automatic Detection", 
                    variable=self.indexing_var, value='auto',
                    command=self.toggle_manual_entry).pack(side='left', padx=5)
        ttk.Radiobutton(selection_frame, text="Manual Entry", 
                    variable=self.indexing_var, value='manual',
                    command=self.toggle_manual_entry).pack(side='left', padx=5)
        
        # Manual entry text area and save button
        self.manual_frame = ttk.Frame(self.resindex_frame)
        self.manual_text = scrolledtext.ScrolledText(self.manual_frame, height=10, width=50)
        self.manual_text.pack(fill='both', expand=True, pady=5)
        
        # Add default template text
        self.manual_text.insert('1.0', """
# Template for ResIndexing.txt:
# For c-type hemes (7 fields):
# CysB CysC HisP Distal Heme c HX
#
# For b-type hemes (5 fields):
# HisP Distal Heme b HX
#
# Where X is the ligand code:
# H = His, M = Met, C = Cys, Y = Tyr
# D = Asp, E = Glu, N = Asn, Q = Gln, K = Lys
""")
        # Save button
        NavButton(self.manual_frame, text="Save ResIndexing.txt",
                command=self.save_resindexing).pack(pady=5)
        
        # 2. Process Residue Indexing controls (includes Redox State Selection)
        self.process_frame = ttk.Frame(container)
        # Redox state table will be created dynamically when needed
        self.redox_frame = None  # Will be created when needed

        # 3. TLeap Input controls
        self.tleap_frame = ttk.Frame(container)
        
        # Output prefix
        prefix_frame = ttk.Frame(self.tleap_frame)
        prefix_frame.pack(fill='x', pady=2)
        ttk.Label(prefix_frame, text="Output Prefix:").pack(side='left', padx=5)
        ttk.Entry(prefix_frame, textvariable=self.prefix_var, style='Beige.TEntry').pack(side='left', padx=5)

        # Solvent type
        solv_frame = ttk.Frame(self.tleap_frame)
        solv_frame.pack(fill='x', pady=2)
        ttk.Label(solv_frame, text="Solvent Type:").pack(side='left', padx=5)
        ttk.Radiobutton(solv_frame, text="Explicit", 
                    variable=self.solvent_var, value='exp',
                    command=self.toggle_solvent_options).pack(side='left', padx=5)
        ttk.Radiobutton(solv_frame, text="Implicit", 
                    variable=self.solvent_var, value='imp',
                    command=self.toggle_solvent_options).pack(side='left', padx=5)

        # Box settings (for explicit solvent)
        self.box_frame = ttk.LabelFrame(self.tleap_frame, text="Box Settings", padding=5)
        box_type_frame = ttk.Frame(self.box_frame)
        box_type_frame.pack(fill='x', pady=2)
        ttk.Label(box_type_frame, text="Box Type:").pack(side='left', padx=5)
        ttk.Radiobutton(box_type_frame, text="Rectangular", 
                    variable=self.box_type_var, value='rec').pack(side='left', padx=5)
        ttk.Radiobutton(box_type_frame, text="Octahedral", 
                    variable=self.box_type_var, value='octahed').pack(side='left', padx=5)

        buffer_frame = ttk.Frame(self.box_frame)
        buffer_frame.pack(fill='x', pady=2)
        ttk.Label(buffer_frame, text="Buffer Size (Å):").pack(side='left', padx=5)
        ttk.Entry(buffer_frame, textvariable=self.buffer_var, width=8, style='Beige.TEntry').pack(side='left', padx=5)

        ion_frame = ttk.Frame(self.box_frame)
        ion_frame.pack(fill='x', pady=2)
        ttk.Label(ion_frame, text="Na+ count:").pack(side='left', padx=5)
        ttk.Entry(ion_frame, textvariable=self.na_var, width=5, style='Beige.TEntry').pack(side='left', padx=5)
        ttk.Label(ion_frame, text="Cl- count:").pack(side='left', padx=5)
        ttk.Entry(ion_frame, textvariable=self.cl_var, width=5, style='Beige.TEntry').pack(side='left', padx=5)

        # 4. Reorder Structure controls
        self.reorder_frame = ttk.Frame(container)
        # Simple frame - just needs space for the buttons which are handled by the main frame

        # 5. Validate Titratable Sites controls
        self.validate_sites_frame = ttk.Frame(container)
        # Selection interface will be created dynamically when needed
        self.validation_frame = None  # Will be created when needed

        # 6. Generate CPIN controls
        self.cpin_frame = ttk.Frame(container)
        # Simple frame - just needs space for the buttons which are handled by the main frame

        # 7. Structure Relaxation controls
        self.relax_frame = ttk.Frame(container)
        # Minimization method selection
        method_frame = ttk.Frame(self.relax_frame)
        method_frame.pack(fill='x', pady=5)
        
        ttk.Label(method_frame, text="Minimization Method:").pack(side='left', padx=5)
        self.min_method_var = tk.StringVar(value='sander')
        ttk.Radiobutton(method_frame, text="SANDER", 
                        variable=self.min_method_var, 
                        value='sander',
                        command=self.toggle_processor_entry).pack(side='left', padx=5)
        ttk.Radiobutton(method_frame, text="PMEMD", 
                        variable=self.min_method_var, 
                        value='pmemd',
                        command=self.toggle_processor_entry).pack(side='left', padx=5)
        
        # Processor count entry
        self.processor_frame = ttk.Frame(self.relax_frame)
        self.processor_frame.pack(fill='x', pady=5)
        ttk.Label(self.processor_frame, text="Number of Processors:").pack(side='left', padx=5)
        self.processor_var = tk.StringVar(value='1')
        vcmd = (self.register(self.validate_processor_count), '%P')
        self.processor_entry = ttk.Entry(self.processor_frame, 
                                    textvariable=self.processor_var,
                                    width=5,
                                    style='Beige.TEntry',
                                    validate='key',
                                    validatecommand=vcmd)
        self.processor_entry.pack(side='left', padx=5)

    def show_step(self, step_idx):
        """Show the specified step and update navigation."""
#       print(f"DEBUG: Entering show_step with step_idx = {step_idx}")
#       print(f"DEBUG: Current step = {self.current_step}")

        self.current_step = step_idx
        
        # Reset buttons
        self.start_button['state'] = 'normal'
        self.next_button['state'] = 'disabled'
        
        # Hide all control frames
        all_frames = [
            self.resindex_frame,
            self.process_frame,
            self.tleap_frame,
            self.reorder_frame,
            self.validate_sites_frame,
            self.cpin_frame,
            self.relax_frame
        ]
        for frame in all_frames:
            frame.pack_forget()
        
        # Hide dynamic frames if they exist
        if hasattr(self, 'redox_frame') and self.redox_frame:
            self.redox_frame.pack_forget()
        if hasattr(self, 'validation_frame') and self.validation_frame:
            self.validation_frame.pack_forget()
        if hasattr(self, 'box_frame'):
            self.box_frame.pack_forget()
        
        # Update step-specific widgets and info           
        if step_idx == 0:  # Create Residue Indexing
            print("DEBUG: Inside step 0 section")
            self.start_button['text'] = "Start Processing"
            
            # Pack the resindex frame and selection frame
#           print("DEBUG: About to pack resindex_frame")
            self.resindex_frame.pack(in_=self.controls_container, fill='both', expand=True)
#           print("DEBUG: resindex_frame packed")
    
            # Make sure the auto/manual selection frame is packed
#           print("DEBUG: Children of resindex_frame:", self.resindex_frame.winfo_children())
            selection_frame = [child for child in self.resindex_frame.winfo_children() 
                            if isinstance(child, ttk.Frame)][0]
            selection_frame.pack(fill='x', pady=5)
            
            # Set initial instruction text
            self.instructions_text['text'] = """
We need to create a file (ResIndexing.txt) that identifies the IDs of the proximal His and distal ligands bound to the heme group. Currently BioDC supports, in principle, b- and c-type hemes with His-Dist ligation, where Dist = His, Met, Cys, Tyr, Asp, Glu, Asn, Gln, and Lys. However, forcefield parameters are currently only available for His-His and His-Met ligation. Stay tuned! 

BioDC will auto-"magically" detect whether your system has b- or c-type hemes (they can be mixed together) and the type of ligation by scanning distances 2–4 Å away from the Fe, CAB and CAC atoms in increments of 0.1 Å. The program stops at the minimum distance needed to classify a heme. This distance approach may fail in some cases. For these, you can manually tell BioDC the residue IDs of the Fe ligands and (if present) thioehter-linked Cys by switching to Manual mode. 

For automatic detection, press 'Start Processing.' For manual entry, use the toggle switch, fill out the template file, and press save. Either way, you can press 'Next Step' at the bottom to proceed."""

            # Hide manual frame initially
            if hasattr(self, 'manual_frame'):
                self.manual_frame.pack_forget()

        elif step_idx == 1:  # Process Residue Indexing
            self.start_button['text'] = "Process Structure"
            self.process_frame.pack(in_=self.controls_container, fill='both', expand=True)
            self.instructions_text['text'] = """
Your ResIndexing file will be used to re-label the residues according to the appropriate parameterizations of b- and c-type His-[Dist] ligated hemes.

Press 'Process Structure' to start the process. You will then be asked to set the redox state of each heme by double-clicking on table entries. 
Once you're satisfied with the selections, press 'Apply Redox States' and then 'Next Step' to continue."""
            
        elif step_idx == 2:  # Generate TLeap Input
            self.start_button['text'] = "Generate and Run TLeap"
            self.tleap_frame.pack(in_=self.controls_container, fill='x', pady=5)
            if self.solvent_var.get() == 'exp':
                self.box_frame.pack(in_=self.tleap_frame, fill='x', pady=5)
            self.instructions_text['text'] = """
Now we will submit the processed PDB to TLEaP to generate topology and coordinate files.

Please enter an output file name (without extension) and set the details for the solvent and simulation box."""
            
        elif step_idx == 3:  # Reorder Structure
            self.start_button['text'] = "Reorder Structure"
            self.reorder_frame.pack(in_=self.controls_container, fill='x', pady=5)
            self.instructions_text['text'] = """
We will use CPPTRAJ to reorder the residues by connectivity for consistency with other AMBER software.

TLEaP, in the prior step, required all the heme groups to come after all the protein residues, which does not follow connectivity for multi-chain protein systems. The ordering will be corrected once you press 'Reorder Structure.' Press 'Next Step' after that to continue."""
            
        elif step_idx == 4:  # Validate Titratable Sites
            self.start_button['text'] = "Validate Titratable Sites"
            self.validate_sites_frame.pack(in_=self.controls_container, fill='both', expand=True)
            
            # Check if we have any titratable residues
            spr_frame = self.master.master.master
            ph_data = spr_frame.frames['pH'].get_data()
            has_titratable = any(ph_data.values())
            
            if has_titratable:
                self.instructions_text['text'] = """
Remember that you selected titratable residues by residue ID earlier. Well, we need to check that those residue IDs are still true since we reordered the structure.

We asked earlier about residues that should be treated as titratable because the residue names for selected ASP, GLU, and HIS need to be changed to AS4, GL4, and HIP during the structure processing. It was also necessary during structure processing to put all the protein residues before any of the heme residues to oblige TLEaP in getting AMBER topology and coordinate files. But now that we've reordered the structure, we need to check, and if necessary, re-select titratable residues. Also, during processing, the propionic acid groups of the hemes were split off as their own residues, and you can now select them to be titratable if you wish. 
Simply press 'Validate Titratable Sites', make your selections by double-clicking and pressing 'Confirm Selections', and then 'Next Step' to continue."""
            else:
                self.instructions_text['text'] = """
No titratable residues were selected, skipping validation."""

                self.next_button['state'] = 'normal'
                self.start_button['state'] = 'disabled'
                self.status_vars[4].set("✓")
            
        elif step_idx == 5:  # Generate CPIN
            self.start_button['text'] = "Generate CPIN"
            self.cpin_frame.pack(in_=self.controls_container, fill='x', pady=5)
            
            # Check if we have any titratable residues
            spr_frame = self.master.master.master
            ph_data = spr_frame.frames['pH'].get_data()
            has_titratable = any(ph_data.values())
            
            if has_titratable:
                self.instructions_text['text'] = """
Because you selected some residues to be titratable, we will generate the CPIN file for constant pH molecular dynamics.
Press 'Generate CPIN' and then 'Next Step' to continue."""
            else:
                self.instructions_text['text'] = """
No titratable residues selected, skipping CPIN generation."""

                self.next_button['state'] = 'normal'
                self.start_button['state'] = 'disabled'
                self.status_vars[5].set("✓")
            
        elif step_idx == 6:  # Structure Relaxation
            self.start_button['text'] = "Start Relaxation"
            self.relax_frame.pack(in_=self.controls_container, fill='x', pady=5)
            self.instructions_text['text'] = """
Finally, we will perform energy minimization on the structure.

Please select the minimization engine and number of processors. Note that the use of SANDER here is only for a single processor."""
        
            # Hide the Next Step button completely
            self.next_button.pack_forget()
        
        # Update navigation
        self.progress_label['text'] = f"Step {step_idx + 1}/{len(self.steps)}"
        print(f"\nMoving to step: {self.steps[step_idx]}")

    def validate_indexing_method(self, method):
        """Validate indexing method input."""
        return method in ['auto', 'man']

    def start_current_step(self):
        """Start processing current step."""
        if self.current_step == 0:  # Create Residue Indexing
            self.create_residue_indexing()
        elif self.current_step == 1:  # Process Residue Indexing
            self.process_residues()
        elif self.current_step == 2:  # Generate TLeap Input
            self.run_tleap()
        elif self.current_step == 3:  # Reorder Structure
            self.reorder_structure()
        elif self.current_step == 4:  # Validate Titratable Sites
            self.validate_titratable_sites()
        elif self.current_step == 5:  # Generate CPIN
            self.generate_cpin()
        elif self.current_step == 6:  # Structure Relaxation
            self.relax_structure()

    def create_residue_indexing(self):
        """Handle residue indexing creation."""
        try:
            spr_frame = self.master.master.master
            init_data = spr_frame.frames['Initialize'].get_data()
            launch_dir = spr_frame.frames['Initialize'].launch_dir
            pdb_path = Path(init_data['pdb_file'])
            
            # Ensure we're in SPR directory
            spr_dir = launch_dir / "SPR"
            spr_dir.mkdir(exist_ok=True)
            os.chdir(spr_dir)
            
            if self.indexing_var.get() == 'auto':
                indexing = ResidueIndexing(pdb_path, launch_dir)
                indexing.analyze_heme_environments()
                
                if indexing.validate_environments():
                    indexing.write_indexing_files()
                    
                    stdout = sys.stdout
                    string_io = StringIO()
                    sys.stdout = string_io
                    
                    indexing.print_environment_summary()
                    
                    sys.stdout = stdout
                    print(string_io.getvalue())
                    
                    print("""
Please verify the identified residues in the ResIndexing.txt file.
If corrections are needed, save them in CorrectedResIndexing.txt.

Note: The program identifies ligands based on distance criteria.
Please check that the assigned ligands are correct and make any necessary adjustments in CorrectedResIndexing.txt.""")
                    
                    self.status_vars[0].set("✓")
                    self.next_button['state'] = 'normal'
                    self.start_button['state'] = 'disabled'
                    
                else:
                    raise ValueError("Failed to validate heme environments")
                    
            else:  # Manual
                print("""
    To create ResIndexing.txt by hand:
    > Create a txt file with an editor of your choosing (e.g. 
        vi ResIndexing.txt). 
        
    > For c-type hemes, each line needs 7 space-separated fields:
        CysB CysC HisP Distal Heme c HX
        
    > For b-type hemes, each line needs 5 space-separated fields:
        HisP Distal Heme b HX
        
    Where X is the code for the distal ligand:
        H = His, M = Met, C = Cys, Y = Tyr, D = Asp,
        E = Glu, N = Asn, Q = Gln, K = Lys
        
    See documentation for detailed format explanation.

    Click 'Next Step' once you have created the file in the SPR directory.""")
                self.next_button['state'] = 'normal'
                self.status_vars[0].set("✓")
                
        except Exception as e:
            print(f"\nError during residue indexing: {str(e)}")
            self.status_vars[0].set("✗")
            messagebox.showerror("Error", f"Failed to create residue indexing: {str(e)}")

    def toggle_manual_entry(self):
        """Show/hide manual entry area based on selection."""
        if self.indexing_var.get() == 'manual':
            self.manual_frame.pack(fill='both', expand=True)
        else:
            self.manual_frame.pack_forget()

    def save_resindexing(self):
        """Save manual entry text to ResIndexing.txt."""
        try:
            # Get text content
            content = self.manual_text.get('1.0', tk.END).strip()
            
            # Remove template comments if they're still there
            if content.startswith("# Template"):
                # Find the first non-comment line
                lines = content.split('\n')
                content_lines = [line for line in lines if not line.strip().startswith('#')]
                content = '\n'.join(content_lines)
            
            if not content.strip():
                messagebox.showerror("Error", "Please enter the residue indexing content.")
                return
                
            # Ensure SPR directory exists
            spr_frame = self.master.master.master
            launch_dir = spr_frame.frames['Initialize'].launch_dir
            spr_dir = launch_dir / "SPR"
            spr_dir.mkdir(exist_ok=True)
            
            # Save file
            with open(spr_dir / "ResIndexing.txt", 'w') as f:
                f.write(content)
                
            print("\nResIndexing.txt saved successfully in SPR directory.")
            
            # Mark step as complete
            self.status_vars[0].set("✓")
            self.next_button['state'] = 'normal'
            
        except Exception as e:
            print(f"\nError saving ResIndexing.txt: {str(e)}")
            messagebox.showerror("Error", f"Failed to save ResIndexing.txt: {str(e)}")

    def process_residues(self):
        """Handle residue processing step."""
        try:
            spr_frame = self.master.master.master
            init_data = spr_frame.frames['Initialize'].get_data()
            launch_dir = spr_frame.frames['Initialize'].launch_dir
            pdb_path = init_data['pdb_file']  # This will get the current active PDB
            
            # Update file naming to use the current PDB name
            pdb_name = Path(pdb_path).stem
            print(f"\nStarting residue processing using PDB: {pdb_name}")
            disulf_data = spr_frame.frames['Disulfide'].get_data()
            ph_data = spr_frame.frames['pH'].get_data()
            
            # Create new redox selection frame if needed
            if not self.redox_frame:
                self.redox_frame = RedoxSelectionFrame(self.process_frame)
            
                # Create the process button inside the redox_frame
                button_frame = ttk.Frame(self.redox_frame)
                button_frame.pack(fill='x', pady=5)
                
                self.process_button = NavButton(
                    button_frame, 
                    text="Apply Redox states",
                    command=self._process_with_redox
                )
                self.process_button.pack(pady=5)

            # Load hemes from ResIndexing.txt
            self.redox_frame.pack(in_=self.process_frame, fill='x', pady=5)
            self.redox_frame.load_hemes(Path.cwd() / "ResIndexing.txt")
                                    
            # Disable the next/start buttons until processing is complete
            self.next_button['state'] = 'disabled'
            self.start_button['state'] = 'disabled'
            
        except Exception as e:
            print(f"\nError during redox state selection: {str(e)}")
            self.status_vars[1].set("✗")
            messagebox.showerror("Error", f"Failed to setup redox selection: {str(e)}")

    def _process_with_redox(self):
        """Continue processing after redox states are selected."""
        # Store original stdout at the start
        stdout = sys.stdout
        string_io = StringIO()
        
        try:
            spr_frame = self.master.master.master
            init_data = spr_frame.frames['Initialize'].get_data()
            launch_dir = spr_frame.frames['Initialize'].launch_dir
            pdb_path = init_data['pdb_file']  # This will get the current active PDB
            pdb_name = Path(pdb_path).stem
            disulf_data = spr_frame.frames['Disulfide'].get_data()
            ph_data = spr_frame.frames['pH'].get_data()
            
            # Get redox states
            redox_states = self.redox_frame.get_redox_states()
            print(f"Processing with redox states: {redox_states}")
            
            # Redirect stdout to capture process output
            sys.stdout = string_io
            
            # Process the residues
            processor = PDBProcessor(Path(pdb_path))  # Use full path to current PDB
            
            # Add redox states to input_dict
            if not hasattr(self, 'input_dict'):
                self.input_dict = {}
            self.input_dict['RedoxState'] = redox_states
            
            # Process disulfides if any
            if disulf_data['disulf_res_list']:
                processor.process_disulfides(disulf_data['disulf_res_list'])
            
            # Process titratable residues
            processor.process_titratable_residues(
                ph_data['asp_ids'],
                ph_data['glu_ids'],
                ph_data['his_ids']
            )
            
            # Calculate heme shifts and move hemes
            indexing_file = Path.cwd() / "ResIndexing.txt"
            processor.heme_shifts = processor.calculate_heme_shifts(indexing_file)
            processor.shift_hemes()
            
            # Process heme environments and collect environment data
            heme_environments = []
            indexing_data = []
            
            with open(indexing_file) as f:
                for line in f:
                    parts = line.strip().split()
                    if not parts or len(parts) not in [5, 7]:
                        continue
                        
                    try:
                        is_c_type = len(parts) == 7
                        heme_id = int(parts[4 if is_c_type else 2])
                        shifted_id = processor.heme_shifts[heme_id]
                        
                        # Get redox state
                        redox_state = 'ox' if redox_states[len(indexing_data)] == 'O' else 'red'
                        
                        # Create environment for processor
                        env = HemeEnvironment(
                            heme_id=heme_id,
                            shifted_id=shifted_id,
                            new_id=shifted_id,
                            his_p=int(parts[2 if is_c_type else 0]),
                            distal_ligand=int(parts[3 if is_c_type else 1]),
                            distal_type='HIS',  # Will be mapped based on ligand code
                            is_c_type=is_c_type,
                            cys_b=int(parts[0]) if is_c_type else None,
                            cys_c=int(parts[1]) if is_c_type else None,
                            redox_state=redox_state
                        )
                        processor.process_heme_environment(env)
                        heme_environments.append(env)
                        
                        # Create indexing data entry
                        data = {
                            'heme_id': heme_id,
                            'shifted_id': shifted_id,
                            'new_id': shifted_id,
                            'his_p': int(parts[2 if is_c_type else 0]),
                            'distal_ligand': int(parts[3 if is_c_type else 1]),
                            'distal_type': parts[-1],  # Keep original code (HH, HM, etc)
                            'is_c_type': is_c_type,
                            'cys_b': int(parts[0]) if is_c_type else None,
                            'cys_c': int(parts[1]) if is_c_type else None,
                            'redox_state': redox_state
                        }
                        indexing_data.append(data)
                        
                    except Exception as e:
                        print(f"Error processing line {line}: {e}")
                        continue
            
            # Save final structure
            output_path = Path.cwd() / "processed.pdb"
            processor.save_structure(output_path)
            
            # Store structure info for tleap
            structure_counts = processor.count_heme_types()
            
            # Transform structure info into format needed by generate_tleap
            self.structure_info = {}
            for heme_type, states in structure_counts.items():
                for redox_state, count in states.items():
                    if count > 0:  # Only include types that exist
                        self.structure_info[(heme_type, redox_state)] = count
                        
            # Add disulfide bond data to indexing data
            if disulf_data['disulf_pair_id']:
                indexing_data.append({
                    'type': 'disulfide',
                    'pairs': disulf_data['disulf_pair_id']
                })

            # Save indexing data
            self.indexing_data = indexing_data
            
            # Validate
            if processor.validate_structure(indexing_file):
                self.status_vars[1].set("✓")  # Set Process Residue Indexing as complete
                print("\nStructure processing completed successfully")
            else:
                raise ValueError("Structure validation failed")
                
            # Restore stdout and update output
            sys.stdout = stdout
            print(string_io.getvalue())
            
            # Clean up UI
            self.redox_frame.pack_forget()
            self.next_button['state'] = 'normal'
            self.start_button['state'] = 'disabled'
            
        except Exception as e:
            # Make sure we restore stdout before handling the error
            sys.stdout = stdout
            print(f"\nError during residue processing: {str(e)}")
            self.status_vars[1].set("✗")
            messagebox.showerror("Error", f"Failed to process residues: {str(e)}")
            return
            
        finally:
            # Always restore stdout
            sys.stdout = stdout

    def next_step(self):
        """Move to next processing step."""
        if self.current_step < len(self.steps) - 1:
            self.show_step(self.current_step + 1)
    
    def validate(self) -> bool:
        return all(var.get() == "✓" for var in self.status_vars)
    
    def get_data(self) -> Dict:
        return {}

    def toggle_solvent_options(self):
        """Show/hide box settings based on solvent type."""
        if self.solvent_var.get() == 'explicit':
            self.box_frame.pack(fill='x', pady=5)
        else:
            self.box_frame.pack_forget()

    def validate_tleap_settings(self) -> bool:
        """Validate TLeap input settings."""
        try:
            if not self.prefix_var.get().strip():
                raise ValueError("Output prefix is required")

            if self.solvent_var.get() == 'explicit':
                buffer_size = float(self.buffer_var.get())
                if buffer_size <= 0:
                    raise ValueError("Buffer size must be positive")

                na_count = int(self.na_var.get())
                cl_count = int(self.cl_var.get())
                if na_count < 0 or cl_count < 0:
                    raise ValueError("Ion counts must be non-negative")

            return True
        except ValueError as e:
            messagebox.showerror("Validation Error", str(e))
            return False

    def run_tleap(self):
        """Handle TLeap input generation and processing."""
        try:
            # Validate TLeap settings first
            if not self.validate_tleap_settings():
                return
                
            # Get necessary data from other frames
            spr_frame = self.master.master.master
            init_data = spr_frame.frames['Initialize'].get_data()
            launch_dir = spr_frame.frames['Initialize'].launch_dir
            
            # Paths
            pdb_path = Path.cwd() / "processed.pdb"
            ff_path = get_biodc_forcefield_dir()
            
            # Verify required data is available
            if not hasattr(self, 'structure_info') or not hasattr(self, 'indexing_data'):
                raise ValueError("Structure processing step must be completed first")

            # Disable buttons and show progress
            self.start_button['state'] = 'disabled'
            self.next_button['state'] = 'disabled'
            print("Generating TLeap input and running...\n")

            def run_tleap_thread():
                try:
                    # Prepare input dictionary for TLeap
                    tleap_input = {
                        'prefix': self.prefix_var.get().strip(),
                        'solvent_type': self.solvent_var.get(),
                        'box_type': self.box_type_var.get() if self.solvent_var.get() == 'explicit' else None,
                        'buffer_size': float(self.buffer_var.get()) if self.solvent_var.get() == 'explicit' else None,
                        'na_count': int(self.na_var.get()) if self.solvent_var.get() == 'explicit' else 0,
                        'cl_count': int(self.cl_var.get()) if self.solvent_var.get() == 'explicit' else 0
                    }

                    output_prefix, solvent_type = generate_tleap.generate_tleap_input(
                        pdb_path.name,  # Just the filename
                        ff_path, 
                        self.structure_info,
                        self.indexing_data,
                        tleap_input,
                        launch_dir
                    )
                    
                    # Store output information for later steps
                    self.output_prefix = output_prefix
                    self.solvent_type = solvent_type
                    
                    # Update UI in main thread
                    self.after(0, self.tleap_success, output_prefix, solvent_type)
                    
                except Exception as e:
                    # Update UI in main thread
                    self.after(0, self.tleap_error, str(e))

            # Start the thread
            thread = threading.Thread(target=run_tleap_thread, daemon=True)
            thread.start()
            
        except Exception as e:
            print(f"\nError during TLeap input generation: {str(e)}")
            self.status_vars[2].set("✗")
            messagebox.showerror("Error", f"Failed to generate TLeap input: {str(e)}")

    def tleap_success(self, output_prefix, solvent_type):
        """Handle successful TLeap execution."""
        self.status_vars[2].set("✓")
        print(  
            f"TLeap input generated successfully.\n\n"
            f"Output prefix: {output_prefix}\n"
            f"Solvent type: {solvent_type}\n\n"
            f"Generated files:\n"
            f"- {output_prefix}.prmtop\n"
            f"- {output_prefix}.rst7"
        )
        
        self.next_button['state'] = 'normal'
        self.start_button['state'] = 'disabled'

    def run_tleap(self):
        """Handle TLeap input generation and processing."""
        try:
            # Validate TLeap settings first
            if not self.validate_tleap_settings():
                return
                
            # Get necessary data from other frames
            spr_frame = self.master.master.master
            init_data = spr_frame.frames['Initialize'].get_data()
            launch_dir = spr_frame.frames['Initialize'].launch_dir
            
            # Paths
            pdb_path = Path.cwd() / "processed.pdb"
            ff_path = get_biodc_forcefield_dir()
            
            # Verify required data is available
            if not hasattr(self, 'structure_info') or not hasattr(self, 'indexing_data'):
                raise ValueError("Structure processing step must be completed first")

            # Disable buttons and show progress
            self.start_button['state'] = 'disabled'
            self.next_button['state'] = 'disabled'
            print("Generating TLeap input and running...\n")

            def run_tleap_thread():
                try:
                    # Prepare complete input dictionary for TLeap with all required fields
                    prefix = self.prefix_var.get().strip()
                    solvent_type = self.solvent_var.get()
                    
                    tleap_input = {
                        'OutPrefix': prefix,  # Match the key expected by InteractionManager
                        'SolvEnv': solvent_type,  # Match the key expected by InteractionManager
                    }
                    
                    # Add explicit solvent settings if needed
                    if solvent_type == 'exp':
                        tleap_input.update({
                            'BoxShape': self.box_type_var.get(),
                            'BufferSize': self.buffer_var.get(),
                            'NaCount': self.na_var.get(),
                            'ClCount': self.cl_var.get()
                        })

                    # Generate and run tleap
                    output_prefix, solvent_type = generate_tleap.generate_tleap_input(
                        pdb_path.name,  # Just the filename
                        ff_path, 
                        self.structure_info,
                        self.indexing_data,
                        tleap_input,  # Pass our complete input dictionary
                        launch_dir
                    )
                    
                    # Store output information for later steps
                    self.output_prefix = output_prefix
                    self.solvent_type = solvent_type
                    
                    # Update UI in main thread
                    self.after(0, self.tleap_success, output_prefix, solvent_type)
                    
                except Exception as e:
                    # Update UI in main thread
                    self.after(0, self.tleap_error, str(e))

            # Start the thread
            thread = threading.Thread(target=run_tleap_thread, daemon=True)
            thread.start()
            
        except Exception as e:
            print(f"\nError during TLeap input generation: {str(e)}")
            self.status_vars[2].set("✗")
            messagebox.showerror("Error", f"Failed to generate TLeap input: {str(e)}")

    def tleap_success(self, output_prefix, solvent_type):
        """Handle successful TLeap execution."""
        self.status_vars[2].set("✓")
        print(  
            f"TLeap input generated successfully.\n\n"
            f"Output prefix: {output_prefix}\n"
            f"Solvent type: {solvent_type}\n\n"
            f"Generated files:\n"
            f"- {output_prefix}.prmtop\n"
            f"- {output_prefix}.rst7"
        )
        
        self.next_button['state'] = 'normal'
        self.start_button['state'] = 'disabled'

    def tleap_error(self, error_msg):
        """Handle TLeap execution error."""
        print(f"\nError during TLeap input generation: {error_msg}")
        self.status_vars[2].set("✗")
        self.start_button['state'] = 'normal'
        messagebox.showerror("Error", f"Failed to generate TLeap input: {error_msg}")

    def reorder_structure(self):
        """Handle structure reordering step."""
        try:
            if not hasattr(self, 'output_prefix'):
                raise ValueError("TLeap output prefix not found. Please complete TLeap step first.")
                
            # Disable button and update status
            self.start_button['state'] = 'disabled'
            print("Starting structure reordering...\n")
            
            def reorder_thread():
                try:
                    # Create cpptraj input
                    reorder_input = f"""parm {self.output_prefix}.prmtop
    trajin {self.output_prefix}.rst7
    fixatomorder parmout {self.output_prefix}_reord.prmtop
    trajout {self.output_prefix}_orig.pdb
    trajout {self.output_prefix}_reord.pdb topresnum
    trajout {self.output_prefix}_reord.rst7 topresnum
    run
    quit"""
                    
                    # Write cpptraj input file
                    with open("ReorderRes.in", "w") as f:
                        f.write(reorder_input)
                    
                    # Run cpptraj
                    result = subprocess.run(
                        "cpptraj -i ReorderRes.in",
                        shell=True,
                        capture_output=True,
                        text=True,
                        check=True
                    )
                    
                    # Store output file paths
                    self.reordered_files = (
                        f"{self.output_prefix}_reord.prmtop",
                        f"{self.output_prefix}_orig.pdb",
                        f"{self.output_prefix}_reord.pdb"
                    )
                    
                    # Store reordered prmtop as final prmtop
                    self.final_prmtop = self.reordered_files[0]
                    
                    # Update UI in main thread
                    self.after(0, self.reorder_success, result.stdout)
                    
                except Exception as e:
                    self.after(0, self.reorder_error, str(e))
            
            # Start the thread
            thread = threading.Thread(target=reorder_thread, daemon=True)
            thread.start()
            
        except Exception as e:
            self.reorder_error(str(e))

    def reorder_success(self, output):
        """Handle successful reordering."""
        print(f"Structure reordering completed successfully.\n\n")
        print(f"Generated files:\n")
        print(f"- {self.reordered_files[0]} (Reordered topology)\n")
        print(f"- {self.reordered_files[1]} (Original PDB)\n")
        print(f"- {self.reordered_files[2]} (Reordered PDB)\n\n")
        print("CPPTRAJ Output:\n")
        print(output)
        
        # Update status and enable next step
        self.status_vars[3].set("✓")
        self.next_button['state'] = 'normal'

    def reorder_error(self, error_msg):
        """Handle reordering error."""
        print(f"Error during structure reordering:\n{error_msg}")
        self.start_button['state'] = 'normal'
        self.status_vars[3].set("✗")
        messagebox.showerror("Error", f"Failed to reorder structure:\n{error_msg}")

    def create_validation_frame(self):
        """Create the validation interface frame."""
        self.validation_frame = ttk.LabelFrame(self, text="Residue Validation", padding=5)
        
        # Status message
        self.validation_message = ttk.Label(self.validation_frame, 
                                        wraplength=400, justify='left')
        self.validation_message.pack(fill='x', pady=5)
        
        # Create residue selection tables
        self.validation_tables = {}
        table_frame = ttk.Frame(self.validation_frame)
        table_frame.pack(fill='both', expand=True)
        
        # Will be populated only when needed
        self.validation_tables = {}
        
        # Confirm button
        self.confirm_button = NavButton(self.validation_frame, 
                                    text="Confirm Selections",
                                    command=self.confirm_validation)
        self.confirm_button.pack(pady=5)

    def create_residue_table(self, parent, res_type: str) -> ttk.Treeview:
        """Create a treeview for residue selection."""
        frame = ttk.Frame(parent)
        frame.pack(side='left', padx=5, pady=5, fill='both', expand=True)
        
        ttk.Label(frame, text=f"{res_type} Residues").pack()
        
        tree = ttk.Treeview(frame, columns=("ID", "Status"), show="headings", height=6)
        tree.heading("ID", text="ID")
        tree.heading("Status", text="Status")
        tree.column("ID", width=50)
        tree.column("Status", width=80)
        tree.pack(fill='x')
        
        # Bind double-click for selection toggle
        tree.bind('<Double-1>', lambda e: self.toggle_residue_selection(e))
        
        return tree

    def toggle_residue_selection(self, event):
        """Toggle residue selection status on double-click."""
        tree = event.widget
        if not tree.selection():
            return
        item = tree.selection()[0]
        values = list(tree.item(item)['values'])
        values[1] = "Selected" if values[1] == "Available" else "Available"
        tree.item(item, values=values)

    def validate_titratable_sites(self):
        """Handle validation of titratable residues after reordering."""
        try:
            if not hasattr(self, 'reordered_files'):
                raise ValueError("Structure must be reordered before validation.")
                
            # Get required data
            spr_frame = self.master.master.master
            launch_dir = spr_frame.frames['Initialize'].launch_dir
            ph_frame = spr_frame.frames['pH']
            ph_data = ph_frame.get_data()
            
            # Print debugging info
            print("\nStarting titratable residue validation...")
            print(f"Launch dir: {launch_dir}")
            print("Current pH selections:", ph_data)
            print(f"Reordered PDB: {self.reordered_files[2]}")
            
            # Check if we have any titratable residues
            has_titratable = any(ph_data.values())
            print(f"Has titratable residues: {has_titratable}")
            
            if not has_titratable:
                print("No titratable residues selected, skipping validation.\n")
                self.status_vars[4].set("✓")
                self.next_button['state'] = 'normal'
                return
            
            # Disable buttons and show progress
            self.start_button['state'] = 'disabled'
            print("Checking titratable residue validity...\n")
            
            def validation_thread():
                try:
                    print("\nStarting validation thread...")
                    
                    # Create ResidueSelector for reordered structure
                    pdb_base = self.reordered_files[2].removesuffix('.pdb')
                    print(f"Creating ResidueSelector for: {pdb_base}")
                    selector = ResidueSelector(pdb_base, launch_dir)
                    
                    # Check each residue type
                    validation_needed = False
                    validation_data = {}
                    
                    residue_types = {
                        "ASP": ("AS4", ph_data['asp_ids']),
                        "GLU": ("GL4", ph_data['glu_ids']),
                        "HIS": ("HIP", ph_data['his_ids']),
                        "LYS": ("LYS", ph_data['lys_ids']),
                        "TYR": ("TYR", ph_data['tyr_ids'])
                    }
                    
                    print("\nChecking each residue type:")
                    for res_type, (proc_name, selected_ids) in residue_types.items():
                        print(f"\n{res_type}:")
                        print(f"  Process name: {proc_name}")
                        print(f"  Selected IDs: {selected_ids}")
                        
                        if selected_ids:
                            # Get currently valid residues
                            print(f"  Finding titratable {proc_name} residues...")
                            current_ids = selector.find_titratable_residues(proc_name)
                            print(f"  Current valid IDs: {current_ids}")
                            
                            selected_set = set(selected_ids.split())
                            print(f"  Selected set: {selected_set}")
                            
                            # Check for invalid selections
                            invalid = [rid for rid in selected_set if rid not in current_ids]
                            print(f"  Invalid selections: {invalid}")
                            
                            if invalid:
                                validation_needed = True
                                validation_data[res_type] = {
                                    'invalid': invalid,
                                    'available': current_ids,
                                    'selected': selected_set - set(invalid)
                                }
                                print(f"  Validation data added for {res_type}")
                    
                    # Check for PRN residues
                    print("\nChecking for PRN residues...")
                    prn_residues = selector.find_titratable_residues("PRN")
                    print(f"PRN residues found: {prn_residues}")
                    if prn_residues:
                        validation_needed = True
                        validation_data['PRN'] = {
                            'invalid': [],
                            'available': prn_residues,
                            'selected': set()
                        }
                        print("PRN validation data added")
                    
                    print(f"\nValidation needed: {validation_needed}")
                    if validation_needed:
                        print("Validation data:", validation_data)
                    
                    # Update UI in main thread
                    print("\nUpdating UI...")
                    self.after(0, self.show_validation_interface 
                            if validation_needed 
                            else self.validation_success, 
                            validation_data)
                    
                except Exception as e:
                    print(f"\nError in validation thread: {str(e)}")
                    import traceback
                    traceback.print_exc()
                    self.after(0, self.validate_error, str(e))
            
            # Start the thread
            print("\nStarting validation thread...")
            thread = threading.Thread(target=validation_thread, daemon=True)
            thread.start()
            
        except Exception as e:
            print(f"\nError in validate_titratable_sites: {str(e)}")
            import traceback
            traceback.print_exc()
            self.validate_error(str(e))

    def show_validation_interface(self, validation_data):
        """Show interface for reselecting invalid residues."""
        print("\nShowing validation interface...")
        try:
            # Create validation frame if needed
            if not self.validation_frame:
                self.validation_frame = ttk.LabelFrame(self.validate_sites_frame, text="Residue Validation", padding=5)
                
                # Status message
                self.validation_message = ttk.Label(self.validation_frame, 
                                                wraplength=400, justify='left')
                self.validation_message.pack(fill='x', pady=5)
                
                # Create confirm button
                self.confirm_button = NavButton(self.validation_frame, 
                                            text="Confirm Selections",
                                            command=self.confirm_validation)
                self.confirm_button.pack(side='bottom', pady=5)
                
                # Dictionary to store validation tables
                self.validation_tables = {}
            
            # Show validation frame
            print("Packing validation frame...")
            self.validation_frame.pack(in_=self.validate_sites_frame, fill='both', expand=True, pady=10)
            
            # Update message
            print("Updating message...")
            message = "Some previously selected residues need to be validated:\n\n"
            for res_type, data in validation_data.items():
                if data['invalid']:
                    message += f"{res_type}: Invalid residues {', '.join(data['invalid'])}\n"
                elif res_type == 'PRN':
                    message += f"PRN residues are now available for selection.\n"
            message += "\nPlease review and update selections below."
            self.validation_message.config(text=message)
            
            print("Creating/updating tables...")
            # Create/update tables for each residue type needing validation
            for res_type, data in validation_data.items():
                print(f"\nProcessing {res_type}:")
                if res_type not in self.validation_tables:
                    print(f"Creating new table for {res_type}")
                    self.validation_tables[res_type] = self.create_residue_table(
                        self.validation_frame, res_type
                    )
                
                # Clear existing items
                tree = self.validation_tables[res_type]
                for item in tree.get_children():
                    tree.delete(item)
                
                # Add available residues
                print(f"Available residues: {data['available']}")
                print(f"Selected residues: {data['selected']}")
                for res_id in sorted(data['available'], key=int):
                    status = "Selected" if res_id in data['selected'] else "Available"
                    tree.insert('', 'end', values=(res_id, status))
            
            # Store validation data for reference
            self.validation_data = validation_data
            print("\nValidation interface setup complete")
            
        except Exception as e:
            print(f"\nError in show_validation_interface: {str(e)}")
            import traceback
            traceback.print_exc()
            self.validate_error(str(e))

    def confirm_validation(self):
        """Handle validation confirmation."""
        print("\n=== CONFIRM VALIDATION CALLED ===")
        
        # Collect new selections
        new_selections = {}
        for res_type, tree in self.validation_tables.items():
            print(f"\nProcessing {res_type} tree:")
            selected = []
            for item in tree.get_children():
                values = tree.item(item)['values']
                print(f"  Residue {values[0]}: Status = {values[1]}")
                if values[1] == "Selected":
                    selected.append(str(values[0]))
            
            if selected:
                new_selections[res_type] = " ".join(selected)
                print(f"  Selected {res_type} residues: {new_selections[res_type]}")
        
        print("\nFinal New Selections:")
        for res_type, selection in new_selections.items():
            print(f"  {res_type}: {selection}")
        
        # If no selections were made
        if not new_selections:
            print("NO RESIDUES SELECTED")
            self.validation_frame.pack_forget()
            self.validation_success({})
            return

        # Update pH frame with new selections
        spr_frame = self.master.master.master
        ph_frame = spr_frame.frames['pH']
        
        # Map residue types to pH frame variables
        type_map = {
            "ASP": "asp_ids",
            "GLU": "glu_ids",
            "HIS": "his_ids",
            "LYS": "lys_ids",
            "TYR": "tyr_ids",
            "PRN": "prn_ids"
        }
        
        # Update values
        for res_type, selection in new_selections.items():
            print(f"\nProcessing {res_type}")
            if res_type in type_map:
                attribute_name = type_map[res_type]
                current_value = getattr(ph_frame.sites, attribute_name)
                print(f"  Current {attribute_name}: '{current_value}'")
                setattr(ph_frame.sites, attribute_name, selection)
                updated_value = getattr(ph_frame.sites, attribute_name)
                print(f"  Updated {attribute_name}: '{updated_value}'")
        
        # Hide validation frame
        self.validation_frame.pack_forget()
        
        # Show success message
        self.validation_success(new_selections)

    def validation_success(self, validation_data):
        """Handle successful validation."""
        print("Titratable residue validation completed successfully.\n\n")
        
        # Add selection summary
        for res_type, selection in validation_data.items():
            if isinstance(selection, str) and selection:  # If it's a selection string
                print(f"Selected {res_type} residues: {selection}\n")
            elif isinstance(selection, dict) and selection.get('selected'):  # If it's validation data
                selected = " ".join(sorted(selection['selected']))
                print(f"Selected {res_type} residues: {selected}\n")
        
        # Update status and enable next step
        self.status_vars[4].set("✓")
        self.next_button['state'] = 'normal'

    def validate_error(self, error_msg):
        """Handle validation error."""
        print(f"Error during titratable residue validation:\n{error_msg}")
        self.start_button['state'] = 'normal'
        self.status_vars[4].set("✗")
        messagebox.showerror("Error", f"Failed to validate titratable residues:\n{error_msg}")

    def validate_processor_count(self, value):
        """Validate processor count entry."""
        if not value:  # Allow empty value
            return True
        try:
            count = int(value)
            return count > 0
        except ValueError:
            return False

    def toggle_processor_entry(self):
        """Enable/disable processor entry based on minimization method."""
        if self.min_method_var.get() == 'sander':
            self.processor_entry['state'] = 'disabled'
            self.processor_var.set('1')
        else:
            self.processor_entry['state'] = 'normal'

    def generate_cpin(self):
        """Handle CPIN file generation."""
        print("\n === GENERATE CPIN METHOD CALLED ===")
        
        try:
            print("Checking for reordered files...")
            if not hasattr(self, 'reordered_files'):
                print("ERROR: No reordered files found")
                raise ValueError("Structure must be reordered before generating CPIN.")
            
            print(f"Reordered files: {self.reordered_files}")

            if not hasattr(self, 'reordered_files'):
                raise ValueError("Structure must be reordered before generating CPIN.")
                
            # Get required data
            spr_frame = self.master.master.master
            ph_frame = spr_frame.frames['pH']
            ph_data = ph_frame.get_data()
            
            # Check if we have any titratable residues
            has_titratable = any(ph_data.values())
            if not has_titratable:
                self.status_vars[4].set("✓")
                self.next_button['state'] = 'normal'
                return
            
            # Create PreparedStructure object
            prep = PreparedStructure(
                pdb=self.reordered_files[2],  # reordered PDB
                sel_asp_ids=ph_data['asp_ids'],
                sel_glu_ids=ph_data['glu_ids'],
                sel_his_ids=ph_data['his_ids'],
                sel_lys_ids=ph_data['lys_ids'],
                sel_tyr_ids=ph_data['tyr_ids'],
                sel_prn_ids=ph_data['prn_ids']
            )
            
            # Disable button and update status
            self.start_button['state'] = 'disabled'
            print("Generating CPIN file...\n")
            
            def cpin_thread():
                try:
                    # Generate CPIN file
                    cpin_file = prep_generate_cpin(
                        self.output_prefix,
                        self.reordered_files[0],  # reordered prmtop
                        prep
                    )
                    
                    if cpin_file:
                        # Update final prmtop
                        self.final_prmtop = f"{self.output_prefix}_new.prmtop"
                        self.after(0, self.cpin_success, cpin_file)
                    else:
                        self.after(0, self.cpin_error, "No titratable residues to process")
                    
                except Exception as e:
                    self.after(0, self.cpin_error, str(e))
            
            # Start the thread
            thread = threading.Thread(target=cpin_thread, daemon=True)
            thread.start()
            
        except Exception as e:
            self.cpin_error(str(e))

    def cpin_success(self, cpin_file):
        """Handle successful CPIN generation."""
        print(f"CPIN file generation completed successfully.\n\n")
        print(f"Generated files:\n")
        print(f"- {cpin_file}\n")
        print(f"- {self.final_prmtop}\n")
        print(f"- {self.output_prefix}_new.rst7\n")
        
        # Update status and enable next step
        self.status_vars[5].set("✓")
        self.next_button['state'] = 'normal'

    def cpin_error(self, error_msg):
        """Handle CPIN generation error."""
        print(f"Error during CPIN generation:\n{error_msg}")
        self.start_button['state'] = 'normal'
        self.status_vars[4].set("✗")
        messagebox.showerror("Error", f"Failed to generate CPIN file:\n{error_msg}")

    def relax_structure(self):
        """Handle structure relaxation with real-time min.out tracking."""
        try:
            # Get required data
            spr_frame = self.master.master.master
            launch_dir = spr_frame.frames['Initialize'].launch_dir
            
            if not hasattr(self, 'output_prefix'):
                raise ValueError("Output prefix not found. Please complete previous steps first.")
            
            # Validate processor count if PMEMD selected
            if self.min_method_var.get() == 'pmemd':
                try:
                    nproc = int(self.processor_var.get())
                    if nproc < 1:
                        raise ValueError("Number of processors must be at least 1")
                except ValueError:
                    messagebox.showerror("Error", "Please enter a valid number of processors")
                    return
                
            # Disable button and update status
            self.start_button['state'] = 'disabled'
            print("Starting structure relaxation...\n")
            
            def monitor_minout(filepath, stop_event, console_widget):
                """Monitor and print min.out contents in real-time, like tail -f.
                Handles both new files and file overwrites."""
                try:
                    print(f"Monitoring file: {filepath}")
                    
                    # Wait for file to exist
                    while not os.path.exists(filepath):
                        if stop_event.is_set():
                            return
                        time.sleep(0.5)
                    
                    # Get initial file size
                    initial_size = os.path.getsize(filepath)
                    last_size = initial_size
                    
                    with open(filepath, 'r') as f:
                        while not stop_event.is_set():
                            current_size = os.path.getsize(filepath)
                            
                            # File was truncated/overwritten
                            if current_size < last_size:
                                f.seek(0)
                                last_size = 0
                            # New content was added
                            elif current_size > last_size:
                                f.seek(last_size)
                                
                            # Read new content
                            new_content = f.read()
                            if new_content:
                                def update_console(text):
                                    console_widget.configure(state='normal')
                                    console_widget.insert(tk.END, text)
                                    console_widget.see(tk.END)
                                    console_widget.configure(state='disabled')
                                console_widget.after(0, update_console, new_content)
                                last_size = f.tell()
                            
                            time.sleep(0.1)
                                
                except Exception as e:
                    print(f"Error monitoring min.out: {e}")
                    import traceback
                    traceback.print_exc()

            def relax_thread():
                try:
                    # Create structure directory
                    struc_dir = launch_dir / "SPR"
                    struc_dir.mkdir(exist_ok=True)
                    
                    # Prepare to monitor min.out
                    minout_path = struc_dir / "min.out"
                    stop_monitoring = threading.Event()
                    
                    # Start min.out monitoring thread
                    monitor_thread = threading.Thread(
                        target=monitor_minout, 
                        args=(str(minout_path), stop_monitoring, self.console),
                        daemon=True
                    )
                    monitor_thread.start()
                    
                    # Prepare input dict with relaxation settings
                    input_dict = {
                        'StructRelaxCompChoice': 'S' if self.min_method_var.get() == 'sander' else 'P'
                    }
                    
                    if self.min_method_var.get() == 'pmemd':
                        input_dict['NProc'] = self.processor_var.get()
                    
                    try:
                        # Run relaxation
                        struct_relax(
                            launch_dir,
                            struc_dir,
                            self.output_prefix,
                            self.solvent_type,
                            input_dict
                        )
                        
                        # Stop monitoring thread
                        stop_monitoring.set()
                        monitor_thread.join(timeout=2)
                        
                        # Update UI in main thread
                        self.after(0, self.relax_success)
                        
                    except Exception as e:
                        # Stop monitoring thread
                        stop_monitoring.set()
                        monitor_thread.join(timeout=2)
                        
                        # Update UI in main thread
                        self.after(0, self.relax_error, str(e))
                    
                except Exception as e:
                    # Update UI in main thread
                    self.after(0, self.relax_error, str(e))
            
            # Start the thread
            thread = threading.Thread(target=relax_thread, daemon=True)
            thread.start()
            
        except Exception as e:
            self.relax_error(str(e))

    def relax_success(self):
        """Handle successful structure relaxation."""
        print("Structure relaxation completed successfully.\n\n")
        print("The structure preparation is now complete!\n")

        # Generate and load final PDB
        final_pdb = self.generate_final_pdb()
        if final_pdb:
            spr_frame = self.master.master.master
            vis_frame = spr_frame.frames['Visualization']
            vis_frame.current_pdb = str(final_pdb)
            vis_frame.load_file(str(final_pdb))
        
        # Switch to visualization frame
        spr_frame.show_step(5)  # Assuming Visualization is step 6 but 0-based
    
        # Update status and enable next step
        self.status_vars[6].set("✓")
        self.next_button['state'] = 'normal'

    def relax_error(self, error_msg):
        """Handle relaxation error."""
        print(f"Error during structure relaxation:\n{error_msg}")
        self.start_button['state'] = 'normal'
        self.status_vars[5].set("✗")
        messagebox.showerror("Error", f"Failed to relax structure:\n{error_msg}")

    def generate_final_pdb(self):
        """Convert minimized structure to PDB using cpptraj"""
        try:
            cpptraj_input = f"""parm {self.output_prefix}_reord.prmtop
    trajin min.rst7
    trajout min.pdb
    run"""
            
            with open("convert_min.in", "w") as f:
                f.write(cpptraj_input)
                
            result = subprocess.run(
                "cpptraj -i convert_min.in",
                shell=True,
                capture_output=True,
                text=True,
                check=True
            )
            
            return Path.cwd() / "min.pdb"
            
        except Exception as e:
            print(f"\nError converting minimized structure: {str(e)}")
            return None

class VisualizationFrame(ttk.LabelFrame):
    def __init__(self, parent, console):
        super().__init__(parent, text="Structure Visualization", padding=10)
        self.console = console
        self.current_pdb = None
        
        # Initialize PyMOL with proper settings
        if PYMOL_AVAILABLE:
            # Initialize PyMOL with proper settings
            pymol.finish_launching(['pymol', '-qc'])
            cmd.set('ray_trace_mode', 1)
            cmd.bg_color('white')
        else:
            print("\nPyMOL is not installed. Structure visualization will be limited.")
            print("To install PyMOL, use one of the following methods:")
            print("  - conda install -c conda-forge pymol-open-source")
            print("  - sudo apt-get install pymol (Ubuntu/Debian)")
            print("  - brew install pymol (MacOS)")
        
        self.create_widgets()
        self.setup_mouse_interactions()
   
    def create_widgets(self):
        # Control frame for visualization options
        control_frame = ttk.Frame(self)
        control_frame.pack(fill='x', pady=5)

        if PYMOL_AVAILABLE:
            # Representation buttons
            # Representation buttons
            rep_frame = ttk.LabelFrame(control_frame, text="Representations", padding=5)
            rep_frame.pack(fill='x', pady=5)
        
            representations = [
                ("Ribbon", self.show_ribbon),
                ("Cartoon", self.show_cartoon),
                ("Surface", self.show_surface),
                ("Highlight Hemes", self.highlight_hemes)
            ]
        
            for text, command in representations:
                ttk.Button(rep_frame, text=text, command=command).pack(side='left', padx=5)
        
            # File load button
            CustomButton(rep_frame, text="Load File", command=self.load_file).pack(side='left', padx=5)
        
            # View control frame
            view_frame = ttk.LabelFrame(control_frame, text="View Controls", padding=5)
            view_frame.pack(fill='x', pady=5)
        
            # Rotation controls
            rot_frame = ttk.Frame(view_frame)
            rot_frame.pack(side='left', padx=10)
        
            ttk.Button(rot_frame, text="↺", command=lambda: self.rotate('y', -30)).pack(side='left', padx=2)
            ttk.Button(rot_frame, text="↻", command=lambda: self.rotate('y', 30)).pack(side='left', padx=2)
            ttk.Button(rot_frame, text="⟲", command=lambda: self.rotate('x', -30)).pack(side='left', padx=2)
            ttk.Button(rot_frame, text="⟳", command=lambda: self.rotate('x', 30)).pack(side='left', padx=2)
        
            # Zoom controls
            zoom_frame = ttk.Frame(view_frame)
            zoom_frame.pack(side='left', padx=10)
        
            ttk.Button(zoom_frame, text="Zoom In", command=lambda: self.zoom(0.8)).pack(side='left', padx=2)
            ttk.Button(zoom_frame, text="Zoom Out", command=lambda: self.zoom(1.2)).pack(side='left', padx=2)
            ttk.Button(zoom_frame, text="Reset View", command=self.reset_view).pack(side='left', padx=2)
        
            # Image display area
            self.image_label = ttk.Label(self)
            self.image_label.pack(fill='both', expand=True)

        else:
            # Show message about PyMOL being required
            msg_frame = ttk.LabelFrame(self, text="PyMOL Required", padding=10)
            msg_frame.pack(fill='both', expand=True, pady=20)

            ttk.Label(msg_frame,
                    text="PyMOL is required for structure visualization.\n\n"
                         "To install PyMOL:\n"
                        "- Using conda:\n"
                        "    conda install -c conda-forge pymol-open-source\n"
                        "- On Ubuntu/Debian:\n"
                        "    sudo apt-get install pymol\n"
                        "- On MacOS:\n"
                        "    brew install pymol",
                    justify='left',
                    wraplength=400).pack(pady=20)

    def setup_mouse_interactions(self):
        if not PYMOL_AVAILABLE:
            # Add a label explaining PyMOL is needed for interactions
            ttk.Label(self, text="PyMOL is required for interactive visualization.\n"
                            "Please install PyMOL to enable structure viewing.",
                    wraplength=400).pack(pady=20)
            return

        def on_mouse_press(event):
            self.image_label.last_x = event.x
            self.image_label.last_y = event.y

        def on_mouse_motion(event):
            if not hasattr(self.image_label, 'last_x'):
                return

            dx = event.x - self.image_label.last_x
            dy = event.y - self.image_label.last_y

            # One-finger drag for translation, shift+drag for rotation
            if event.state & 0x1:  # Shift pressed
                print("rotation")
                cmd.rotate('y', -dx * 0.5)
                cmd.rotate('x', -dy * 0.5)
                self.update_display()

            else:  # Default to translation
                print("One-finger drag translaiton")
                cmd.translate([dx * 5, -dy * 5, 0], "protein")
                cmd.move('x', dx)
                cmd.move('y', -dy)
                self.update_display()

            self.image_label.last_x = event.x
            self.image_label.last_y = event.y

        def on_touchpad_zoom(event):
            # Convert trackpad pinch gestures to zoom
            try:
                # Get the zoom factor from the event
                delta = event.delta
                if hasattr(event, 'num') and event.num == 2:
                    delta = -delta
                
                scale = 1.0 + (delta / 120.0) * 0.1
                cmd.zoom('all', scale)
                self.update_display()
            except Exception as e:
                print(f"Zoom error: {e}")

        # Bind events
        self.image_label.bind('<ButtonPress-1>', on_mouse_press)
        self.image_label.bind('<B1-Motion>', on_mouse_motion)
        
        # Bind both mousewheel and gesture events for zoom
        self.image_label.bind('<MouseWheel>', on_touchpad_zoom)
        self.image_label.bind('<Control-MouseWheel>', on_touchpad_zoom)
        self.image_label.bind('<Command-MouseWheel>', on_touchpad_zoom)

    def rotate(self, axis, angle):
        """Rotate the view around specified axis by given angle."""
        if not PYMOL_AVAILABLE:
            return
        if not self.current_pdb:
            return
        try:
            cmd.rotate(axis, angle)
            self.update_display()
        except Exception as e:
            print(f"\nError rotating view: {str(e)}")

    def zoom(self, factor):
        """Zoom the view by the given factor."""
        if not PYMOL_AVAILABLE:
            return
        if not self.current_pdb:
            return
        try:
            cmd.zoom('all', factor)
            self.update_display()
        except Exception as e:
            print(f"\nError zooming view: {str(e)}")

    def reset_view(self):
        """Reset the view to default."""
        if not PYMOL_AVAILABLE:
            return
        if not self.current_pdb:
            return
        try:
            cmd.reset()
            cmd.zoom('all')
            self.update_display()
        except Exception as e:
            print(f"\nError resetting view: {str(e)}")

    def load_file(self, filepath=None):
        if not PYMOL_AVAILABLE:
            messagebox.showinfo("PyMOL Not Available", 
                "PyMOL is not installed. Structure visualization is not available.\n\n"
                "To install PyMOL:\n"
                "- Using conda:\n    conda install -c conda-forge pymol-open-source\n"
                "- On Ubuntu/Debian:\n    sudo apt-get install pymol\n"
                "- On MacOS:\n    brew install pymol")
            return

        if not filepath:
            filepath = filedialog.askopenfilename(
            title="Select PDB File", 
            filetypes=[('PDB Files', '*.pdb')]
        )
   
        if filepath:
            try:
                self.current_pdb = filepath
                cmd.delete('all')
                cmd.load(filepath, 'protein')
#               cmd.bg_color('black')
                cmd.set('depth_cue', 0)
                cmd.set('ray_trace_fog', 0)

                cmd.set('cartoon_fancy_helices', 1)
                cmd.set('cartoon_highlight_color', 'grey60')
                cmd.set('cartoon_transparency', 0.4)

                # Basic protein representation
                cmd.hide('everything', 'all')
                cmd.show_as('cartoon', 'protein')
                cmd.color('grey60', 'protein')

                # Color definitions for heme types and ligands 
                heme_colors = {
                    # c-type hemes
                    'HCO': {'color': 'tv_red', 'name': 'His-His oxidized'},
                    'HCR': {'color': 'salmon', 'name': 'His-His reduced'},
                    'MCO': {'color': 'orange', 'name': 'His-Met oxidized'}, 
                    'MCR': {'color': 'wheat', 'name': 'His-Met reduced'},
                    'CCO': {'color': 'yellow', 'name': 'His-Cys oxidized'},
                    'CCR': {'color': 'tv_yellow', 'name': 'His-Cys reduced'},
                    'YCO': {'color': 'forest', 'name': 'His-Tyr oxidized'},
                    'YCR': {'color': 'lime', 'name': 'His-Tyr reduced'},
                    'DCO': {'color': 'marine', 'name': 'His-Asp oxidized'},
                    'DCR': {'color': 'lightblue', 'name': 'His-Asp reduced'},
                    'ECO': {'color': 'purple', 'name': 'His-Glu oxidized'},
                    'ECR': {'color': 'violet', 'name': 'His-Glu reduced'},
                    'NCO': {'color': 'magenta', 'name': 'His-Asn oxidized'}, 
                    'NCR': {'color': 'pink', 'name': 'His-Asn reduced'},
                    'QCO': {'color': 'brown', 'name': 'His-Gln oxidized'},
                    'QCR': {'color': 'sand', 'name': 'His-Gln reduced'},
                    'KCO': {'color': 'teal', 'name': 'His-Lys oxidized'},
                    'KCR': {'color': 'cyan', 'name': 'His-Lys reduced'},

                    # b-type hemes
                    'HBO': {'color': 'slate', 'name': 'His-His oxidized'},
                    'HBR': {'color': 'grey', 'name': 'His-His reduced'},
                    'MBO': {'color': 'red', 'name': 'His-Met oxidized'},
                    'MBR': {'color': 'salmon', 'name': 'His-Met reduced'}, 
                    'CBO': {'color': 'orange', 'name': 'His-Cys oxidized'},
                    'CBR': {'color': 'wheat', 'name': 'His-Cys reduced'},
                    'YBO': {'color': 'yellow', 'name': 'His-Tyr oxidized'},
                    'YBR': {'color': 'tv_yellow', 'name': 'His-Tyr reduced'},
                    'DBO': {'color': 'forest', 'name': 'His-Asp oxidized'},
                    'DBR': {'color': 'lime', 'name': 'His-Asp reduced'},
                    'EBO': {'color': 'marine', 'name': 'His-Glu oxidized'},
                    'EBR': {'color': 'lightblue', 'name': 'His-Glu reduced'},
                    'NBO': {'color': 'purple', 'name': 'His-Asn oxidized'},
                    'NBR': {'color': 'violet', 'name': 'His-Asn reduced'},
                    'QBO': {'color': 'magenta', 'name': 'His-Gln oxidized'},
                    'QBR': {'color': 'pink', 'name': 'His-Gln reduced'},
                    'KBO': {'color': 'brown', 'name': 'His-Lys oxidized'},
                    'KBR': {'color': 'cyan', 'name': 'His-Lys reduced'}
                }

                # Proximal and distal ligands mapping
                heme_ligand_pairs = {
                    # c-type hemes with proximal His
                    'HCO': {'prox': 'PHO', 'dist': 'DHO'},  # His-His
                    'HCR': {'prox': 'PHR', 'dist': 'DHR'},
                    'MCO': {'prox': 'PMO', 'dist': 'DMO'},  # His-Met  
                    'MCR': {'prox': 'PMR', 'dist': 'DMR'},
                    'CCO': {'prox': 'PCO', 'dist': 'DCO'},  # His-Cys
                    'CCR': {'prox': 'PCR', 'dist': 'DCR'},
                    'YCO': {'prox': 'PYO', 'dist': 'DYO'},  # His-Tyr
                    'YCR': {'prox': 'PYR', 'dist': 'DYR'},
                    'DCO': {'prox': 'PDO', 'dist': 'DDO'},  # His-Asp
                    'DCR': {'prox': 'PDR', 'dist': 'DDR'},
                    'ECO': {'prox': 'PEO', 'dist': 'DEO'},  # His-Glu
                    'ECR': {'prox': 'PER', 'dist': 'DER'},
                    'NCO': {'prox': 'PNO', 'dist': 'DNO'},  # His-Asn
                    'NCR': {'prox': 'PNR', 'dist': 'DNR'},
                    'QCO': {'prox': 'PQO', 'dist': 'DQO'},  # His-Gln
                    'QCR': {'prox': 'PQR', 'dist': 'DQR'},
                    'KCO': {'prox': 'PKO', 'dist': 'DKO'},  # His-Lys
                    'KCR': {'prox': 'PKR', 'dist': 'DKR'},

                    # b-type hemes with proximal His
                    'HBO': {'prox': 'FHO', 'dist': 'RHO'},  # His-His
                    'HBR': {'prox': 'FHR', 'dist': 'RHR'},
                    'MBO': {'prox': 'FMO', 'dist': 'RMO'},  # His-Met
                    'MBR': {'prox': 'FMR', 'dist': 'RMR'},
                    'CBO': {'prox': 'FCO', 'dist': 'RCO'},  # His-Cys
                    'CBR': {'prox': 'FCR', 'dist': 'RCR'},
                    'YBO': {'prox': 'FYO', 'dist': 'RYO'},  # His-Tyr
                    'YBR': {'prox': 'FYR', 'dist': 'RYR'},
                    'DBO': {'prox': 'FDO', 'dist': 'RDO'},  # His-Asp
                    'DBR': {'prox': 'FDR', 'dist': 'RDR'},
                    'EBO': {'prox': 'FEO', 'dist': 'REO'},  # His-Glu
                    'EBR': {'prox': 'FER', 'dist': 'RER'},
                    'NBO': {'prox': 'FNO', 'dist': 'RNO'},  # His-Asn
                    'NBR': {'prox': 'FNR', 'dist': 'RNR'},
                    'QBO': {'prox': 'FQO', 'dist': 'RQO'},  # His-Gln
                    'QBR': {'prox': 'FQR', 'dist': 'RQR'},
                    'KBO': {'prox': 'FKO', 'dist': 'RKO'},  # His-Lys
                    'KBR': {'prox': 'FKR', 'dist': 'RKR'}
                }

                # Color hemes and their specific ligands
                for heme_type, ligands in heme_ligand_pairs.items():
                    # Select and show heme
                    cmd.select('temp_heme', f'resn {heme_type}')
                    cmd.show_as('sticks', 'temp_heme')
                    cmd.color(heme_colors[heme_type]['color'], 'temp_heme')

                    # Select and show proximal ligand
                    cmd.select('prox_ligand', f'resn {ligands["prox"]}')
                    cmd.show_as('licorice', 'prox_ligand') 
                    cmd.color(heme_colors[heme_type]['color'], 'prox_ligand')

                    # Select and show distal ligand
                    cmd.select('dist_ligand', f'resn {ligands["dist"]}')
                    cmd.show_as('licorice', 'dist_ligand')
                    cmd.color(heme_colors[heme_type]['color'], 'dist_ligand')

                    # Create labels
#                   cmd.label('temp_heme', f'"{heme_colors[heme_type]["name"]}"')

                # Clean up selections
                cmd.delete('temp_heme')
                cmd.delete('prox_ligand')
                cmd.delete('dist_ligand')
                cmd.deselect()
            
                # Set view
                cmd.zoom('all')
                self.update_display()

            except Exception as e:
                print(f"\nError loading structure: {str(e)}")
                messagebox.showerror("Error", f"Failed to load structure: {str(e)}")

    def update_display(self):
        """Capture PyMOL view and update display."""
        if not self.current_pdb:
            return
            
        try:
            # Set view parameters
            cmd.zoom('all')
            
            # Save image
            temp_image = "temp_view.png"
            cmd.ray(800, 600)  # Ray trace for better quality
            cmd.png(temp_image, width=800, height=600)
            
            # Load and display image
            image = Image.open(temp_image)
            photo = ImageTk.PhotoImage(image)
            self.image_label.configure(image=photo)
            self.image_label.image = photo  # Keep reference
            
            # Clean up
            os.remove(temp_image)
            
            return self.image_label
                     
        except Exception as e:
            print(f"\nError updating display: {str(e)}")
            messagebox.showerror("Error", f"Failed to update display: {str(e)}")
    
    def show_ribbon(self):
        """Show ribbon representation."""
        if not PYMOL_AVAILABLE:
            return
        if not self.current_pdb:
            return
        try:
            cmd.hide('everything', 'all')
            cmd.show_as('ribbon', 'protein')
            cmd.color('marine', 'protein')
            self.update_display()
        except Exception as e:
            print(f"\nError showing ribbon: {str(e)}")
    
    def show_cartoon(self):
        """Show cartoon representation."""
        if not PYMOL_AVAILABLE:
            return
        if not self.current_pdb:
            return
        try:
            cmd.hide('everything', 'all')
            cmd.show_as('cartoon', 'protein')
            cmd.color('marine', 'protein')
            self.update_display()
        except Exception as e:
            print(f"\nError showing cartoon: {str(e)}")
    
    def show_surface(self):
        """Show surface representation."""
        if not PYMOL_AVAILABLE:
            return
        if not self.current_pdb:
            return
        try:
            cmd.hide('everything', 'all')
            cmd.show_as('surface', 'protein')
            cmd.color('white', 'protein')
            self.update_display()
        except Exception as e:
            print(f"\nError showing surface: {str(e)}")
    
    # def highlight_hemes(self):
    #     """Highlight heme groups."""
    #     if not self.current_pdb:
    #         return
    #     try:
    #         cmd.hide('everything', 'all')
    #         cmd.show_as('cartoon', 'protein')
    #         cmd.color('marine', 'protein')
    #         cmd.select('hemes', 'resn HEM')
    #         cmd.show_as('sticks', 'hemes')
    #         cmd.color('red', 'hemes')
    #         cmd.deselect()
    #         self.update_display()
    #     except Exception as e:
    #         print(f"\nError highlighting hemes: {str(e)}")

    def highlight_hemes(self):
        """Highlight heme groups."""
        if not PYMOL_AVAILABLE:
            return
        if not self.current_pdb:
            return
        try:
            cmd.hide('everything', 'all')
            cmd.show_as('cartoon', 'protein')
            cmd.color('marine', 'protein')

            # Iterate through heme types and create selections
            for heme_type, ligands in heme_ligand_pairs.items():
                cmd.select(f'hemes_{heme_type}', f'resn {heme_type}')
                cmd.show_as('sticks', f'hemes_{heme_type}')
                cmd.color(heme_colors[heme_type]['color'], f'hemes_{heme_type}')

            # Deselect all
            cmd.deselect()
            self.update_display()
        except Exception as e:
            print(f"\nError highlighting hemes: {str(e)}")    

    def __del__(self):
        """Cleanup PyMOL when frame is destroyed."""
        if PYMOL_AVAILABLE:
            try:
                cmd.quit()
            except:
                pass

def main():
    try:
        root = tk.Tk()
        root.title("BioDC: Structure Preparation & Relaxation")
        
        # Set window size
        root.geometry("1400x800")
        
        # Create and pack the main frame
        app = SPRFrame(root)
        app.pack(expand=True, fill='both')
        
        # Proper cleanup on window close
        def on_closing():
            if messagebox.askokcancel("Quit", "Do you want to quit?"):
                if PYMOL_AVAILABLE:
                    cmd.quit()
                root.destroy()
                
        root.protocol("WM_DELETE_WINDOW", on_closing)
        
        root.mainloop()
        
    except Exception as e:
        print(f"Error starting application: {str(e)}")
        sys.exit(1)

if __name__ == '__main__':
    main()
