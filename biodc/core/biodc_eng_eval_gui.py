import tkinter as tk
from tkinter import ttk, filedialog
from pathlib import Path
import tkinter as tk
from tkinter import ttk, messagebox, filedialog, scrolledtext
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from typing import Dict, Any

# Import calculators and interaction manager
from biodc.utils.interaction import InteractionManager
from biodc.core.engeval_modules.lambda_calculator import LambdaCalculator
from biodc.core.engeval_modules.dg_calculator import DeltaGCalculator
from biodc.core.engeval_modules.interaction_calculator import HemeInteractionCalculator
from biodc.core.engeval_modules.hda_calculator import CouplingCalculator
from biodc.core.engeval_modules.cooperativity_analyzer import HemeCooperativityAnalyzer

class EnergeticEvaluationGUI:
    def __init__(self, interaction_manager, launch_dir, forcefield_dir, pdb_file):
        """
        Initialize the Energetic Evaluation GUI
        
        Args:
            interaction_manager: Manager for user interactions
            launch_dir: Project launch directory
            forcefield_dir: Directory with forcefield files
            pdb_file: Path to input PDB file
        """
        self.interaction_manager = interaction_manager
        self.launch_dir = launch_dir
        self.forcefield_dir = forcefield_dir
        self.pdb_file = pdb_file
        
        # Initialize calculators
        self.lambda_calculator = LambdaCalculator(
            interaction_manager, launch_dir)
        self.delta_g_calculator = DeltaGCalculator(
            interaction_manager, launch_dir, forcefield_dir, pdb_file)
        self.interaction_calculator = HemeInteractionCalculator(
            interaction_manager, launch_dir, forcefield_dir, pdb_file)
        self.coupling_calculator = CouplingCalculator(
            interaction_manager, launch_dir)
        self.cooperativity_analyzer = HemeCooperativityAnalyzer(
            interaction_manager, launch_dir / "EE")
        
        # Workflow state tracking
        self.results = {
            'lambda': None,
            'delta_g': None,
            'interactions': None,
            'coupling': None,
            'cooperativity': None
        }
        
        # Create main window and frames
        self.create_main_window()
        
    def create_main_window(self):
        """Create the main application window with workflow steps"""
        self.root = tk.Tk()
        self.root.title("Energetic Evaluation")
        
        # Notebook for step-by-step workflow
        self.notebook = ttk.Notebook(self.root)
        
        # Create workflow frames
        self.frames = {
            'Sequence': SequenceSelectionFrame(self.notebook, self),
            'Reorganization': ReorganizationEnergyFrame(self.notebook, self),
            'ReactionFreeEnergy': ReactionFreeEnergyFrame(self.notebook, self),
            'Interactions': InteractionsFrame(self.notebook, self),
            'Coupling': CouplingFrame(self.notebook, self),
            'Cooperativity': CooperativityFrame(self.notebook, self),
            'Visualization': VisualizationFrame(self.notebook, self)
        }
        
        # Add frames to notebook
        for name, frame in self.frames.items():
            self.notebook.add(frame, text=name)
        
        self.notebook.pack(expand=True, fill='both')
        
        # Navigation buttons
        nav_frame = ttk.Frame(self.root)
        nav_frame.pack(fill='x', pady=10)
        
        ttk.Button(nav_frame, text="Previous", command=self.previous_step).pack(side='left', padx=5)
        ttk.Button(nav_frame, text="Next", command=self.next_step).pack(side='right', padx=5)
        
        # Status bar
        self.status_var = tk.StringVar(value="Ready")
        ttk.Label(nav_frame, textvariable=self.status_var).pack(side='bottom')

    def next_step(self):
        """Move to next workflow step"""
        current_tab = self.notebook.select()
        tabs = list(self.frames.keys())
        current_index = list(self.frames.values()).index(current_tab)
        
        if current_index < len(tabs) - 1:
            # Validate current step before moving
            current_frame = self.frames[tabs[current_index]]
            if hasattr(current_frame, 'validate') and not current_frame.validate():
                return
            
            # Move to next tab
            next_tab = list(self.frames.values())[current_index + 1]
            self.notebook.select(next_tab)
            self.status_var.set(f"Step: {tabs[current_index + 1]}")

    def previous_step(self):
        """Move to previous workflow step"""
        current_tab = self.notebook.select()
        tabs = list(self.frames.keys())
        current_index = list(self.frames.values()).index(current_tab)
        
        if current_index > 0:
            previous_tab = list(self.frames.values())[current_index - 1]
            self.notebook.select(previous_tab)
            self.status_var.set(f"Step: {tabs[current_index - 1]}")

class SequenceSelectionFrame(ttk.Frame):
    def __init__(self, parent, ee_gui):
        super().__init__(parent)
        self.ee_gui = ee_gui
        
        # Heme sequence selection methods
        self.sequence_var = tk.StringVar()
        self.manual_var = tk.StringVar()
        
        # Auto-detection button
        ttk.Button(self, text="Auto-detect Heme Sequence", 
                   command=self.auto_detect_sequence).pack(pady=10)
        
        # Manual entry
        manual_frame = ttk.LabelFrame(self, text="Manual Sequence Entry")
        manual_frame.pack(pady=10)
        
        ttk.Label(manual_frame, text="Enter Heme Residue IDs:").pack()
        ttk.Entry(manual_frame, textvariable=self.manual_var).pack()
        ttk.Button(manual_frame, text="Set Sequence", 
                   command=self.set_manual_sequence).pack()
        
        # Sequence display
        self.sequence_label = ttk.Label(self, text="Selected Sequence: ")
        self.sequence_label.pack(pady=10)
        
        # Structural analysis preview
        self.preview_frame = ttk.LabelFrame(self, text="Sequence Analysis")
        self.preview_frame.pack(pady=10)
        
    def auto_detect_sequence(self):
        """Automatically detect heme sequence from PDB"""
        try:
            from biodc.utils.structure_analyzer import PDBProcessor
            
            processor = PDBProcessor()
            atoms_dict = processor.read_pdb_atoms(self.ee_gui.pdb_file)
            
            # Display sequence detection options
            sequence = self.ee_gui.cooperativity_analyzer.select_heme_sequence(
                self.ee_gui.pdb_file
            )
            
            # Display sequence and update state
            self.sequence_var.set(' → '.join(map(str, sequence)))
            self.sequence_label.config(text=f"Selected Sequence: {' → '.join(map(str, sequence))}")
            
            # Update EE GUI state
            self.ee_gui.sequence = sequence
            
            # Show sequence analysis
            self._show_sequence_analysis(sequence)
            
        except Exception as e:
            messagebox.showerror("Sequence Detection Error", str(e))
    
    def set_manual_sequence(self):
        """Set manually entered sequence"""
        try:
            # Parse manual input
            sequence = [int(x) for x in self.manual_var.get().split()]
            
            # Update display and state
            self.sequence_var.set(' → '.join(map(str, sequence)))
            self.sequence_label.config(text=f"Selected Sequence: {' → '.join(map(str, sequence))}")
            
            # Update EE GUI state
            self.ee_gui.sequence = sequence
            
            # Show sequence analysis
            self._show_sequence_analysis(sequence)
            
        except ValueError:
            messagebox.showerror("Invalid Input", "Please enter space-separated integer heme IDs")
    
    def _show_sequence_analysis(self, sequence):
        """Display structural analysis of selected sequence"""
        from biodc.utils.structure_analyzer import PDBProcessor
        
        # Clear previous preview
        for widget in self.preview_frame.winfo_children():
            widget.destroy()
        
        processor = PDBProcessor()
        
        # Coordinates and analysis
        atoms_dict = processor.read_pdb_atoms(self.ee_gui.pdb_file)
        
        # Plane angle between consecutive hemes
        angles = []
        distances = []
        
        for i in range(len(sequence) - 1):
            heme1_coords = atoms_dict[sequence[i]]
            heme2_coords = atoms_dict[sequence[i+1]]
            
            angle = processor.calculate_plane_angle(heme1_coords, heme2_coords)
            distance = processor.calculate_min_distance(heme1_coords, heme2_coords)
            
            angles.append(angle)
            distances.append(distance)
        
        # Create a text display
        text_display = scrolledtext.ScrolledText(self.preview_frame, height=10, width=50, wrap=tk.WORD)
        text_display.pack(padx=5, pady=5)
        
        # Insert analysis text
        text_display.insert(tk.END, "Sequence Structural Analysis:\n")
        text_display.insert(tk.END, f"Total hemes: {len(sequence)}\n")
        text_display.insert(tk.END, f"Sequence: {' → '.join(map(str, sequence))}\n\n")
        
        text_display.insert(tk.END, "Inter-Heme Geometry:\n")
        text_display.insert(tk.END, "-------------------\n")
        
        for i, (angle, distance) in enumerate(zip(angles, distances), 1):
            text_display.insert(tk.END, 
                f"Step {i} (Heme {sequence[i-1]} → Heme {sequence[i]}):\n"
                f"  Plane Angle: {angle:.2f}°\n"
                f"  Edge-to-Edge Distance: {distance:.2f} Å\n"
            )
        
        text_display.config(state=tk.DISABLED)  # Make read-only

    def validate(self):
        """Validate sequence selection"""
        if not hasattr(self.ee_gui, 'sequence'):
            messagebox.showerror("Error", "Please select or enter a heme sequence")
            return False
        return True
    
class ReorganizationEnergyFrame(ttk.Frame):
    def __init__(self, parent, ee_gui):
        super().__init__(parent)
        self.ee_gui = ee_gui
        self.lambda_calculator = ee_gui.lambda_calculator
        
        # Create main layout
        self.create_method_selection()
        self.create_calculator_options()
        self.create_results_display()
        
    def create_method_selection(self):
        """Create method selection section"""
        method_frame = ttk.LabelFrame(self, text="Calculation Method")
        method_frame.pack(fill='x', padx=10, pady=10)
        
        # Method selection variables
        self.method_var = tk.StringVar(value='compute')
        
        # Radio buttons for method selection
        ttk.Radiobutton(method_frame, text="Compute from Structure", 
                        variable=self.method_var, value='compute').pack(side='left', padx=5)
        ttk.Radiobutton(method_frame, text="Manual Entry", 
                        variable=self.method_var, value='manual').pack(side='left', padx=5)
        
    def create_calculator_options(self):
        """Create options for lambda calculation"""
        options_frame = ttk.LabelFrame(self, text="Calculation Options")
        options_frame.pack(fill='x', padx=10, pady=10)
        
        # SASA Backend Selection
        ttk.Label(options_frame, text="SASA Backend:").pack(side='left', padx=5)
        self.backend_var = tk.StringVar(value='vmd')
        backend_dropdown = ttk.Combobox(
            options_frame, 
            textvariable=self.backend_var, 
            values=['vmd', 'freesasa'], 
            state='readonly', 
            width=10
        )
        backend_dropdown.pack(side='left', padx=5)
        
        # Reorganization Parameters
        param_frame = ttk.Frame(options_frame)
        param_frame.pack(side='left', padx=10)
        
        ttk.Label(param_frame, text="α:").pack(side='left')
        self.alpha_var = tk.StringVar(value='5.18')
        ttk.Entry(param_frame, textvariable=self.alpha_var, width=6).pack(side='left', padx=2)
        
        ttk.Label(param_frame, text="β:").pack(side='left')
        self.beta_var = tk.StringVar(value='0.016')
        ttk.Entry(param_frame, textvariable=self.beta_var, width=6).pack(side='left', padx=2)
        
        # Compute button
        ttk.Button(options_frame, text="Calculate", command=self.calculate_lambda).pack(side='right', padx=5)
        
    def create_results_display(self):
        """Create results display area"""
        results_frame = ttk.LabelFrame(self, text="Calculation Results")
        results_frame.pack(fill='both', expand=True, padx=10, pady=10)
        
        # Results text widget
        self.results_text = scrolledtext.ScrolledText(
            results_frame, 
            wrap=tk.WORD, 
            height=10
        )
        self.results_text.pack(fill='both', expand=True, padx=5, pady=5)
        
        # Export and visualization buttons
        button_frame = ttk.Frame(results_frame)
        button_frame.pack(fill='x', padx=5, pady=5)
        
        ttk.Button(button_frame, text="Export Results", command=self.export_results).pack(side='left', padx=5)
        ttk.Button(button_frame, text="Plot Results", command=self.plot_results).pack(side='left', padx=5)
        
    def calculate_lambda(self):
        """Perform lambda calculation based on selected method"""
        try:
            # Validate sequence is available
            if not hasattr(self.ee_gui, 'sequence'):
                messagebox.showerror("Error", "Please select a heme sequence first")
                return
            
            # Determine calculation method
            if self.method_var.get() == 'manual':
                # Manual entry of lambda values
                lambda_values, dielectric_constants = self.lambda_calculator._get_manual_values(
                    len(self.ee_gui.sequence) - 1
                )
            else:
                # Compute from structure
                # Prepare parameters
                params = ReorganizationParameters(
                    alpha=float(self.alpha_var.get()),
                    beta=float(self.beta_var.get())
                )
                
                # Temporarily set parameters on calculator
                original_params = self.lambda_calculator.params
                self.lambda_calculator.params = params
                
                # Compute lambda values
                lambda_values, dielectric_constants = self.lambda_calculator.compute_reorganization_energy(
                    self.ee_gui.sequence, 
                    self.ee_gui.pdb_file
                )
                
                # Restore original parameters
                self.lambda_calculator.params = original_params
            
            # Display results
            self.display_results(lambda_values, dielectric_constants)
            
            # Store results in EE GUI state
            self.ee_gui.results['lambda'] = {
                'values': lambda_values,
                'dielectrics': dielectric_constants
            }
            
        except Exception as e:
            messagebox.showerror("Calculation Error", str(e))
    
    def display_results(self, lambda_values, dielectric_constants):
        """Display calculation results in results text widget"""
        self.results_text.delete(1.0, tk.END)
        
        # Display lambda values
        self.results_text.insert(tk.END, "Reorganization Energy Analysis\n")
        self.results_text.insert(tk.END, "=" * 30 + "\n\n")
        
        for i, (lamb, dielec) in enumerate(zip(lambda_values, dielectric_constants), 1):
            self.results_text.insert(tk.END, 
                f"Step {i}:\n"
                f"  λ (Reorganization Energy): {lamb:.3f} meV\n"
                f"  Dielectric Constant: {dielec:.3f}\n\n"
            )
        
        # Calculate and display summary statistics
        self.results_text.insert(tk.END, "Summary Statistics:\n")
        self.results_text.insert(tk.END, "=" * 20 + "\n")
        self.results_text.insert(tk.END, 
            f"Mean λ: {np.mean(lambda_values):.3f} meV\n"
            f"Min λ: {min(lambda_values):.3f} meV\n"
            f"Max λ: {max(lambda_values):.3f} meV\n"
        )
        
        self.results_text.config(state=tk.DISABLED)  # Make read-only
    
    def export_results(self):
        """Export results to a file"""
        if not hasattr(self.ee_gui, 'results') or 'lambda' not in self.ee_gui.results:
            messagebox.showinfo("Export Error", "No results to export")
            return
        
        # Choose export location
        export_file = filedialog.asksaveasfilename(
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("Text files", "*.txt")]
        )
        
        if export_file:
            try:
                lamb_results = self.ee_gui.results['lambda']
                with open(export_file, 'w') as f:
                    f.write("Step,Reorganization Energy (meV),Dielectric Constant\n")
                    for i, (lamb, dielec) in enumerate(
                        zip(lamb_results['values'], lamb_results['dielectrics']), 1
                    ):
                        f.write(f"{i},{lamb:.3f},{dielec:.3f}\n")
                
                messagebox.showinfo("Export Successful", f"Results exported to {export_file}")
            except Exception as e:
                messagebox.showerror("Export Error", str(e))
    
    def plot_results(self):
        """Create visualization of lambda values"""
        if not hasattr(self.ee_gui, 'results') or 'lambda' not in self.ee_gui.results:
            messagebox.showinfo("Plot Error", "No results to plot")
            return
        
        lamb_results = self.ee_gui.results['lambda']
        lambda_values = lamb_results['values']
        dielectric_constants = lamb_results['dielectrics']
        
        plt.figure(figsize=(10, 6))
        
        # Lambda values subplot
        plt.subplot(1, 2, 1)
        plt.bar(range(1, len(lambda_values) + 1), lambda_values)
        plt.title("Reorganization Energies")
        plt.xlabel("Electron Transfer Step")
        plt.ylabel("λ (meV)")
        plt.grid(axis='y', linestyle='--', alpha=0.7)
        
        # Dielectric constants subplot
        plt.subplot(1, 2, 2)
        plt.bar(range(1, len(dielectric_constants) + 1), dielectric_constants, color='orange')
        plt.title("Dielectric Constants")
        plt.xlabel("Electron Transfer Step")
        plt.ylabel("Dielectric Constant")
        plt.grid(axis='y', linestyle='--', alpha=0.7)
        
        plt.tight_layout()
        plt.show()
    
    def validate(self):
        """Validate the reorganization energy calculation"""
        # Check if results are available
        return hasattr(self.ee_gui, 'results') and 'lambda' in self.ee_gui.results
    
class ReactionFreeEnergyFrame(ttk.Frame):
    def __init__(self, parent, ee_gui):
        super().__init__(parent)
        self.ee_gui = ee_gui
        self.delta_g_calculator = ee_gui.delta_g_calculator
        
        # Create main layout
        self.create_method_selection()
        self.create_calculation_options()
        self.create_results_display()
        
    def create_method_selection(self):
        """Create method selection section"""
        method_frame = ttk.LabelFrame(self, text="Calculation Method")
        method_frame.pack(fill='x', padx=10, pady=10)
        
        # Method selection variables
        self.method_var = tk.StringVar(value='compute')
        
        # Radio buttons for method selection
        ttk.Radiobutton(method_frame, text="Compute with PBSA", 
                        variable=self.method_var, value='compute').pack(side='left', padx=5)
        ttk.Radiobutton(method_frame, text="Manual Entry", 
                        variable=self.method_var, value='manual').pack(side='left', padx=5)
        
    def create_calculation_options(self):
        """Create options for DG calculation"""
        options_frame = ttk.LabelFrame(self, text="Calculation Options")
        options_frame.pack(fill='x', padx=10, pady=10)
        
        # Reference State Selection
        ref_frame = ttk.Frame(options_frame)
        ref_frame.pack(side='left', padx=5)
        ttk.Label(ref_frame, text="Reference State:").pack(side='left')
        self.ref_state_var = tk.StringVar(value='ox')
        ttk.Radiobutton(ref_frame, text="Oxidized", 
                        variable=self.ref_state_var, value='ox').pack(side='left')
        ttk.Radiobutton(ref_frame, text="Reduced", 
                        variable=self.ref_state_var, value='red').pack(side='left')
        
        # Calculation Type Selection
        calc_frame = ttk.Frame(options_frame)
        calc_frame.pack(side='left', padx=10)
        ttk.Label(calc_frame, text="Calculation Type:").pack(side='left')
        self.calc_type_var = tk.StringVar(value='standard')
        ttk.Combobox(calc_frame, 
                     textvariable=self.calc_type_var, 
                     values=['standard', 'delphi', 'membrane'], 
                     state='readonly', 
                     width=10).pack(side='left')
        
        # Additional PBSA Parameters
        param_frame = ttk.Frame(options_frame)
        param_frame.pack(side='left', padx=10)
        
        # Ionic Strength
        ttk.Label(param_frame, text="Ionic Strength (mM):").pack(side='left')
        self.istrng_var = tk.StringVar(value='150.0')
        ttk.Entry(param_frame, textvariable=self.istrng_var, width=8).pack(side='left', padx=2)
        
        # Compute button
        ttk.Button(options_frame, text="Calculate", command=self.calculate_delta_g).pack(side='right', padx=5)
        
    def create_results_display(self):
        """Create results display area"""
        results_frame = ttk.LabelFrame(self, text="Calculation Results")
        results_frame.pack(fill='both', expand=True, padx=10, pady=10)
        
        # Results text widget
        self.results_text = scrolledtext.ScrolledText(
            results_frame, 
            wrap=tk.WORD, 
            height=10
        )
        self.results_text.pack(fill='both', expand=True, padx=5, pady=5)
        
        # Export and visualization buttons
        button_frame = ttk.Frame(results_frame)
        button_frame.pack(fill='x', padx=5, pady=5)
        
        ttk.Button(button_frame, text="Export Results", command=self.export_results).pack(side='left', padx=5)
        ttk.Button(button_frame, text="Plot Results", command=self.plot_results).pack(side='left', padx=5)
        
    def calculate_delta_g(self):
        """Perform delta G calculation based on selected method"""
        try:
            # Validate sequence is available
            if not hasattr(self.ee_gui, 'sequence'):
                messagebox.showerror("Error", "Please select a heme sequence first")
                return
            
            # Ensure lambda values are available if computing
            if (self.method_var.get() == 'compute' and 
                (not hasattr(self.ee_gui.results, 'lambda') or 
                 'values' not in self.ee_gui.results.get('lambda', {}))):
                messagebox.showwarning("Warning", "Reorganization energy calculation recommended before DG")
            
            # Determine calculation method
            if self.method_var.get() == 'manual':
                # Manual entry of DG values
                delta_g_values = self.delta_g_calculator._get_manual_dg_values(
                    self.ee_gui.sequence
                )
            else:
                # Determine dielectric constants
                dielectric_constants = None
                if hasattr(self.ee_gui.results, 'lambda'):
                    dielectric_constants = self.ee_gui.results['lambda'].get('dielectrics')
                
                # Compute delta G
                delta_g_values = self.delta_g_calculator.compute_reaction_free_energy(
                    sequence=self.ee_gui.sequence,
                    dielectric_constants=dielectric_constants,
                    is_cyclic=False  # Default to linear pathway
                )
            
            # Display results
            self.display_results(delta_g_values)
            
            # Store results in EE GUI state
            self.ee_gui.results['delta_g'] = {
                'values': delta_g_values
            }
            
        except Exception as e:
            messagebox.showerror("Calculation Error", str(e))
    
    def display_results(self, delta_g_values):
        """Display calculation results in results text widget"""
        self.results_text.delete(1.0, tk.END)
        
        # Display DG values
        self.results_text.insert(tk.END, "Reaction Free Energy Analysis\n")
        self.results_text.insert(tk.END, "=" * 30 + "\n\n")
        
        for i, dg in enumerate(delta_g_values, 1):
            self.results_text.insert(tk.END, 
                f"Step {i}:\n"
                f"  ΔG: {dg:.3f} eV\n\n"
            )
        
        # Calculate and display summary statistics
        self.results_text.insert(tk.END, "Summary Statistics:\n")
        self.results_text.insert(tk.END, "=" * 20 + "\n")
        self.results_text.insert(tk.END, 
            f"Mean ΔG: {np.mean(delta_g_values):.3f} eV\n"
            f"Min ΔG: {min(delta_g_values):.3f} eV\n"
            f"Max ΔG: {max(delta_g_values):.3f} eV\n"
            f"Total ΔG: {sum(delta_g_values):.3f} eV\n"
        )
        
        self.results_text.config(state=tk.DISABLED)  # Make read-only
    
    def export_results(self):
        """Export results to a file"""
        if not hasattr(self.ee_gui, 'results') or 'delta_g' not in self.ee_gui.results:
            messagebox.showinfo("Export Error", "No results to export")
            return
        
        # Choose export location
        export_file = filedialog.asksaveasfilename(
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("Text files", "*.txt")]
        )
        
        if export_file:
            try:
                delta_g_values = self.ee_gui.results['delta_g']['values']
                with open(export_file, 'w') as f:
                    f.write("Step,Reaction Free Energy (eV)\n")
                    for i, dg in enumerate(delta_g_values, 1):
                        f.write(f"{i},{dg:.6f}\n")
                
                messagebox.showinfo("Export Successful", f"Results exported to {export_file}")
            except Exception as e:
                messagebox.showerror("Export Error", str(e))
    
    def plot_results(self):
        """Create visualization of delta G values"""
        if not hasattr(self.ee_gui, 'results') or 'delta_g' not in self.ee_gui.results:
            messagebox.showinfo("Plot Error", "No results to plot")
            return
        
        delta_g_values = self.ee_gui.results['delta_g']['values']
        
        plt.figure(figsize=(10, 6))
        
        # Bar plot of DG values
        plt.bar(range(1, len(delta_g_values) + 1), delta_g_values)
        plt.title("Reaction Free Energy (ΔG)")
        plt.xlabel("Electron Transfer Step")
        plt.ylabel("ΔG (eV)")
        
        # Add value labels on top of bars
        for i, v in enumerate(delta_g_values):
            plt.text(i + 1, v, f'{v:.3f}', ha='center', va='bottom')
        
        # Highlight positive and negative values
        plt.axhline(y=0, color='r', linestyle='--')
        plt.grid(axis='y', linestyle='--', alpha=0.7)
        
        plt.tight_layout()
        plt.show()
    
    def validate(self):
        """Validate the reaction free energy calculation"""
        # Check if results are available
        return hasattr(self.ee_gui, 'results') and 'delta_g' in self.ee_gui.results
    
class InteractionsFrame(ttk.Frame):
    def __init__(self, parent, ee_gui):
        super().__init__(parent)
        self.ee_gui = ee_gui
        self.interaction_calculator = ee_gui.interaction_calculator
        
        self.create_method_selection()
        self.create_calculation_options()
        self.create_results_display()
        
    def create_method_selection(self):
        """Create method selection section"""
        method_frame = ttk.LabelFrame(self, text="Interaction Calculation Method")
        method_frame.pack(fill='x', padx=10, pady=10)
        
        self.method_var = tk.StringVar(value='compute')
        
        ttk.Radiobutton(method_frame, text="Compute Interactions", 
                        variable=self.method_var, value='compute').pack(side='left', padx=5)
        ttk.Radiobutton(method_frame, text="Manual Entry", 
                        variable=self.method_var, value='manual').pack(side='left', padx=5)
        
    def create_calculation_options(self):
        """Create options for interaction calculation"""
        options_frame = ttk.LabelFrame(self, text="Calculation Options")
        options_frame.pack(fill='x', padx=10, pady=10)
        
        # Parallel processing option
        parallel_frame = ttk.Frame(options_frame)
        parallel_frame.pack(side='left', padx=5)
        
        ttk.Label(parallel_frame, text="Parallel Calculations:").pack(side='left')
        self.parallel_var = tk.StringVar(value='no')
        ttk.Radiobutton(parallel_frame, text="Yes", 
                        variable=self.parallel_var, value='yes').pack(side='left')
        ttk.Radiobutton(parallel_frame, text="No", 
                        variable=self.parallel_var, value='no').pack(side='left')
        
        # Processor count (if parallel)
        self.processors_var = tk.StringVar(value='1')
        proc_frame = ttk.Frame(options_frame)
        proc_frame.pack(side='left', padx=5)
        ttk.Label(proc_frame, text="Processors:").pack(side='left')
        ttk.Entry(proc_frame, textvariable=self.processors_var, width=5).pack(side='left')
        
        # Compute button
        ttk.Button(options_frame, text="Calculate", command=self.calculate_interactions).pack(side='right', padx=5)
        
    def create_results_display(self):
        """Create results display area"""
        results_frame = ttk.LabelFrame(self, text="Interaction Energies")
        results_frame.pack(fill='both', expand=True, padx=10, pady=10)
        
        # Results text widget
        self.results_text = scrolledtext.ScrolledText(results_frame, wrap=tk.WORD, height=10)
        self.results_text.pack(fill='both', expand=True, padx=5, pady=5)
        
        # Export and visualization buttons
        button_frame = ttk.Frame(results_frame)
        button_frame.pack(fill='x', padx=5, pady=5)
        
        ttk.Button(button_frame, text="Export Results", command=self.export_results).pack(side='left', padx=5)
        ttk.Button(button_frame, text="Visualize Matrix", command=self.plot_interaction_matrix).pack(side='left', padx=5)
        
    def calculate_interactions(self):
        """Perform interaction energy calculation"""
        try:
            # Validate sequence
            if not hasattr(self.ee_gui, 'sequence'):
                messagebox.showerror("Error", "Please select a heme sequence first")
                return
            
            # Determine dielectric constants
            dielectric_constants = None
            if hasattr(self.ee_gui.results, 'lambda'):
                dielectric_constants = self.ee_gui.results['lambda'].get('dielectrics')
            
            # Determine parallel processing
            n_parallel = None
            if self.parallel_var.get() == 'yes':
                n_parallel = int(self.processors_var.get())
            
            # Compute or manually enter interactions
            if self.method_var.get() == 'manual':
                # Placeholder for manual interaction entry
                interactions = self.interaction_calculator._get_manual_matrix(self.ee_gui.sequence)
            else:
                # Compute interactions
                interactions = self.interaction_calculator.compute_heme_interactions(
                    sequence=self.ee_gui.sequence,
                    dielectric_constants=dielectric_constants,
                    n_parallel=n_parallel
                )
            
            # Display and store results
            self.display_results(interactions)
            self.ee_gui.results['interactions'] = interactions
            
        except Exception as e:
            messagebox.showerror("Calculation Error", str(e))
    
    def display_results(self, interactions):
        """Display interaction results"""
        self.results_text.delete(1.0, tk.END)
        self.results_text.insert(tk.END, "Heme-Heme Interaction Energies\n")
        self.results_text.insert(tk.END, "=" * 30 + "\n\n")
        
        # Display interaction matrix with summary
        for (heme1, heme2), energy in interactions.items():
            self.results_text.insert(tk.END, 
                f"Heme {heme1} - Heme {heme2}: {energy:.3f} eV\n"
            )
        
        self.results_text.config(state=tk.DISABLED)
    
    def export_results(self):
        """Export interaction results"""
        if not hasattr(self.ee_gui, 'results') or 'interactions' not in self.ee_gui.results:
            messagebox.showinfo("Export Error", "No results to export")
            return
        
        export_file = filedialog.asksaveasfilename(
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv")]
        )
        
        if export_file:
            try:
                interactions = self.ee_gui.results['interactions']
                with open(export_file, 'w') as f:
                    f.write("Heme1,Heme2,Interaction Energy (eV)\n")
                    for (heme1, heme2), energy in interactions.items():
                        f.write(f"{heme1},{heme2},{energy:.6f}\n")
                
                messagebox.showinfo("Export Successful", f"Results exported to {export_file}")
            except Exception as e:
                messagebox.showerror("Export Error", str(e))
    
    def plot_interaction_matrix(self):
        """Visualize interaction matrix"""
        if not hasattr(self.ee_gui, 'results') or 'interactions' not in self.ee_gui.results:
            messagebox.showinfo("Plot Error", "No results to plot")
            return
        
        interactions = self.ee_gui.results['interactions']
        sequence = self.ee_gui.sequence
        
        # Create matrix for visualization
        matrix = np.zeros((len(sequence), len(sequence)))
        for (heme1, heme2), energy in interactions.items():
            idx1 = sequence.index(heme1)
            idx2 = sequence.index(heme2)
            matrix[idx1, idx2] = energy
            matrix[idx2, idx1] = energy  # Symmetric
        
        plt.figure(figsize=(10, 8))
        plt.imshow(matrix, cmap='coolwarm', aspect='auto')
        plt.colorbar(label='Interaction Energy (eV)')
        plt.title("Heme-Heme Interaction Energy Matrix")
        plt.xlabel("Heme Index")
        plt.ylabel("Heme Index")
        plt.xticks(range(len(sequence)), sequence)
        plt.yticks(range(len(sequence)), sequence)
        
        # Annotate each cell with energy value
        for i in range(len(sequence)):
            for j in range(len(sequence)):
                plt.text(j, i, f'{matrix[i, j]:.2f}', 
                         ha='center', va='center', 
                         color='white' if abs(matrix[i, j]) > 0.5 else 'black')
        
        plt.tight_layout()
        plt.show()
    
    def validate(self):
        """Validate interaction calculation"""
        return hasattr(self.ee_gui, 'results') and 'interactions' in self.ee_gui.results

class CouplingFrame(ttk.Frame):
    def __init__(self, parent, ee_gui):
        super().__init__(parent)
        self.ee_gui = ee_gui
        self.coupling_calculator = ee_gui.coupling_calculator
        
        self.create_method_selection()
        self.create_calculation_options()
        self.create_results_display()
        
    def create_method_selection(self):
        """Create method selection section"""
        method_frame = ttk.LabelFrame(self, text="Coupling Calculation Method")
        method_frame.pack(fill='x', padx=10, pady=10)
        
        self.method_var = tk.StringVar(value='geometry')
        
        ttk.Radiobutton(method_frame, text="Compute from Geometry", 
                        variable=self.method_var, value='geometry').pack(side='left', padx=5)
        ttk.Radiobutton(method_frame, text="Manual Entry", 
                        variable=self.method_var, value='manual').pack(side='left', padx=5)
        
    def create_calculation_options(self):
        """Create options for coupling calculation"""
        options_frame = ttk.LabelFrame(self, text="Calculation Options")
        options_frame.pack(fill='x', padx=10, pady=10)
        
        # Stacking type preset
        stacking_frame = ttk.Frame(options_frame)
        stacking_frame.pack(side='left', padx=5)
        
        ttk.Label(stacking_frame, text="Stacking Preset:").pack(side='left')
        self.stacking_var = tk.StringVar(value='default')
        stacking_dropdown = ttk.Combobox(
            stacking_frame, 
            textvariable=self.stacking_var, 
            values=['default', 'slip-stacked', 'T-stacked', 'co-planar'], 
            state='readonly', 
            width=15
        )
        stacking_dropdown.pack(side='left', padx=5)
        
        # Compute button
        ttk.Button(options_frame, text="Calculate", command=self.calculate_coupling).pack(side='right', padx=5)
        
    def create_results_display(self):
        """Create results display area"""
        results_frame = ttk.LabelFrame(self, text="Electronic Coupling Results")
        results_frame.pack(fill='both', expand=True, padx=10, pady=10)
        
        # Results text widget
        self.results_text = scrolledtext.ScrolledText(results_frame, wrap=tk.WORD, height=10)
        self.results_text.pack(fill='both', expand=True, padx=5, pady=5)
        
        # Export and visualization buttons
        button_frame = ttk.Frame(results_frame)
        button_frame.pack(fill='x', padx=5, pady=5)
        
        ttk.Button(button_frame, text="Export Results", command=self.export_results).pack(side='left', padx=5)
        ttk.Button(button_frame, text="Visualize Couplings", command=self.plot_couplings).pack(side='left', padx=5)
        
    def calculate_coupling(self):
        """Perform electronic coupling calculation"""
        try:
            # Validate sequence
            if not hasattr(self.ee_gui, 'sequence'):
                messagebox.showerror("Error", "Please select a heme sequence first")
                return
            
            # Compute or manually enter couplings
            if self.method_var.get() == 'manual':
                couplings = self.coupling_calculator._get_manual_couplings(self.ee_gui.sequence)
            else:
                # Compute from geometry
                couplings = self.coupling_calculator.compute_electronic_coupling(
                    sequence=self.ee_gui.sequence,
                    pdb_file=self.ee_gui.pdb_file
                )
            
            # Display and store results
            self.display_results(couplings)
            self.ee_gui.results['coupling'] = couplings
            
        except Exception as e:
            messagebox.showerror("Calculation Error", str(e))
    
    def display_results(self, couplings):
        """Display coupling results"""
        self.results_text.delete(1.0, tk.END)
        self.results_text.insert(tk.END, "Electronic Coupling Analysis\n")
        self.results_text.insert(tk.END, "=" * 30 + "\n\n")
        
        for i, coupling in enumerate(couplings, 1):
            self.results_text.insert(tk.END, 
                f"Step {i}: {coupling:.3f} meV\n"
            )
        
        # Summary statistics
        self.results_text.insert(tk.END, "\nSummary:\n")
        self.results_text.insert(tk.END, 
            f"Mean Coupling: {np.mean(couplings):.3f} meV\n"
            f"Min Coupling: {min(couplings):.3f} meV\n"
            f"Max Coupling: {max(couplings):.3f} meV\n"
        )
        
        self.results_text.config(state=tk.DISABLED)
    
    def export_results(self):
        """Export coupling results"""
        if not hasattr(self.ee_gui, 'results') or 'coupling' not in self.ee_gui.results:
            messagebox.showinfo("Export Error", "No results to export")
            return
        
        export_file = filedialog.asksaveasfilename(
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv")]
        )
        
        if export_file:
            try:
                couplings = self.ee_gui.results['coupling']
                with open(export_file, 'w') as f:
                    f.write("Step,Electronic Coupling (meV)\n")
                    for i, coupling in enumerate(couplings, 1):
                        f.write(f"{i},{coupling:.6f}\n")
                
                messagebox.showinfo("Export Successful", f"Results exported to {export_file}")
            except Exception as e:
                messagebox.showerror("Export Error", str(e))

    def plot_couplings(self):
        """Visualize coupling values"""
        if not hasattr(self.ee_gui, 'results') or 'coupling' not in self.ee_gui.results:
            messagebox.showinfo("Plot Error", "No results to plot")
            return
        
        couplings = self.ee_gui.results['coupling']
        
        plt.figure(figsize=(10, 6))
        plt.bar(range(1, len(couplings) + 1), couplings)
        plt.title("Electronic Coupling Values")
        plt.xlabel("Electron Transfer Step")
        plt.ylabel("Coupling (meV)")
        plt.grid(axis='y', linestyle='--', alpha=0.7)
        
        # Add value labels on top of bars
        for i, v in enumerate(couplings):
            plt.text(i + 1, v, f'{v:.3f}', ha='center', va='bottom')
        
        plt.tight_layout()
        plt.show()
    
    def validate(self):
        """Validate coupling calculation"""
        return hasattr(self.ee_gui, 'results') and 'coupling' in self.ee_gui.results

class CooperativityFrame(ttk.Frame):
    def __init__(self, parent, ee_gui):
        super().__init__(parent)
        self.ee_gui = ee_gui
        self.cooperativity_analyzer = ee_gui.cooperativity_analyzer
        
        self.create_method_selection()
        self.create_calculation_options()
        self.create_results_display()
        
    def create_method_selection(self):
        """Create method selection section"""
        method_frame = ttk.LabelFrame(self, text="Cooperativity Analysis Method")
        method_frame.pack(fill='x', padx=10, pady=10)
        
        self.method_var = tk.StringVar(value='both')
        
        ttk.Radiobutton(method_frame, text="Sequential and Geometric", 
                        variable=self.method_var, value='both').pack(side='left', padx=5)
        ttk.Radiobutton(method_frame, text="Sequential Only", 
                        variable=self.method_var, value='seq').pack(side='left', padx=5)
        ttk.Radiobutton(method_frame, text="Geometric Only", 
                        variable=self.method_var, value='geo').pack(side='left', padx=5)
        
    def create_calculation_options(self):
        """Create options for cooperativity analysis"""
        options_frame = ttk.LabelFrame(self, text="Analysis Parameters")
        options_frame.pack(fill='x', padx=10, pady=10)
        
        # Interaction level and periodic options
        interaction_frame = ttk.Frame(options_frame)
        interaction_frame.pack(side='left', padx=5)
        
        ttk.Label(interaction_frame, text="Adjacent Only:").pack(side='left')
        self.adjacent_var = tk.StringVar(value='yes')
        ttk.Radiobutton(interaction_frame, text="Yes", 
                        variable=self.adjacent_var, value='yes').pack(side='left')
        ttk.Radiobutton(interaction_frame, text="No", 
                        variable=self.adjacent_var, value='no').pack(side='left')
        
        # Periodic boundary condition
        periodic_frame = ttk.Frame(options_frame)
        periodic_frame.pack(side='left', padx=5)
        ttk.Label(periodic_frame, text="Periodic:").pack(side='left')
        self.periodic_var = tk.StringVar(value='no')
        ttk.Radiobutton(periodic_frame, text="Yes", 
                        variable=self.periodic_var, value='yes').pack(side='left')
        ttk.Radiobutton(periodic_frame, text="No", 
                        variable=self.periodic_var, value='no').pack(side='left')
        
        # Compute button
        ttk.Button(options_frame, text="Analyze", command=self.calculate_cooperativity).pack(side='right', padx=5)
        
    def create_results_display(self):
        """Create results display area"""
        results_frame = ttk.LabelFrame(self, text="Cooperativity Analysis Results")
        results_frame.pack(fill='both', expand=True, padx=10, pady=10)
        
        # Results text widget
        self.results_text = scrolledtext.ScrolledText(results_frame, wrap=tk.WORD, height=10)
        self.results_text.pack(fill='both', expand=True, padx=5, pady=5)
        
        # Export and visualization buttons
        button_frame = ttk.Frame(results_frame)
        button_frame.pack(fill='x', padx=5, pady=5)
        
        ttk.Button(button_frame, text="Export Results", command=self.export_results).pack(side='left', padx=5)
        ttk.Button(button_frame, text="Visualize Cooperativity", command=self.plot_cooperativity).pack(side='left', padx=5)
        
    def calculate_cooperativity(self):
        """Perform cooperativity analysis"""
        try:
            # Validate sequence and prerequisite calculations
            if not hasattr(self.ee_gui, 'sequence'):
                messagebox.showerror("Error", "Please select a heme sequence first")
                return
            
            # Ensure interaction energy matrix is available
            if not hasattr(self.ee_gui.results, 'interactions'):
                messagebox.showwarning("Warning", "Interaction energy calculation recommended")
                return
            
            # Prepare analysis parameters
            params = AnalysisParameters(
                model=self.method_var.get(),
                adjacent_only=self.adjacent_var.get() == 'yes',
                make_periodic=self.periodic_var.get() == 'yes',
                energy_shift=0.0,  # Default values can be refined
                interaction_scale=1.0,
                plot_option='both'
            )
            
            # Compute interactions if not already done
            interactions = self.ee_gui.results.get('interactions', 
                self.ee_gui.interaction_calculator.compute_heme_interactions(
                    sequence=self.ee_gui.sequence,
                    dielectric_constants=self.ee_gui.results.get('lambda', {}).get('dielectrics')
                )
            )
            
            # Perform cooperativity analysis
            results = self.cooperativity_analyzer.analyze_interactions(
                energies=interactions,
                sequence=self.ee_gui.sequence,
                params=params,
                output_dir=self.ee_gui.launch_dir / "EE" / "cooperativity"
            )
            
            # Display and store results
            self.display_results(results)
            self.ee_gui.results['cooperativity'] = results
            
        except Exception as e:
            messagebox.showerror("Calculation Error", str(e))
    
    def display_results(self, results):
        """Display cooperativity analysis results"""
        self.results_text.delete(1.0, tk.END)
        self.results_text.insert(tk.END, "Cooperativity Analysis\n")
        self.results_text.insert(tk.END, "=" * 30 + "\n\n")
        
        # Display results for each model type
        for model_type in ['sequential', 'geometric']:
            if model_type in results:
                self.results_text.insert(tk.END, f"{model_type.capitalize()} Model:\n")
                
                # Oxidation order
                order = results[model_type].get('oxidation_order', [])
                self.results_text.insert(tk.END, 
                    f"  Oxidation Order: {' → '.join(map(str, order))}\n"
                )
                
                # Delta G values
                delta_g = results[model_type].get('delta_G', [])
                total_dg = 0
                for step, (heme1, heme2, dg) in enumerate(delta_g, 1):
                    self.results_text.insert(tk.END, 
                        f"  Step {step}: Heme {heme1} → Heme {heme2}: {dg:.3f} eV\n"
                    )
                    total_dg += dg
                
                self.results_text.insert(tk.END, 
                    f"  Total ΔG: {total_dg:.3f} eV\n\n"
                )
        
        self.results_text.config(state=tk.DISABLED)
    
    def export_results(self):
        """Export cooperativity analysis results"""
        if not hasattr(self.ee_gui, 'results') or 'cooperativity' not in self.ee_gui.results:
            messagebox.showinfo("Export Error", "No results to export")
            return
        
        export_file = filedialog.asksaveasfilename(
            defaultextension=".json",
            filetypes=[("JSON files", "*.json"), ("Text files", "*.txt")]
        )
        
        if export_file:
            try:
                import json
                
                with open(export_file, 'w') as f:
                    json.dump(self.ee_gui.results['cooperativity'], f, indent=2)
                
                messagebox.showinfo("Export Successful", f"Results exported to {export_file}")
            except Exception as e:
                messagebox.showerror("Export Error", str(e))
    
    def plot_cooperativity(self):
        """Visualize cooperativity results"""
        if not hasattr(self.ee_gui, 'results') or 'cooperativity' not in self.ee_gui.results:
            messagebox.showinfo("Plot Error", "No results to plot")
            return
        
        results = self.ee_gui.results['cooperativity']
        
        # Create multiple subplots for different visualizations
        plt.figure(figsize=(15, 10))
        
        # Subplot 1: Delta G Landscape
        plt.subplot(2, 2, 1)
        for model_type in ['sequential', 'geometric']:
            if model_type in results:
                delta_g = [x[2] for x in results[model_type].get('delta_G', [])]
                plt.bar(range(1, len(delta_g) + 1), delta_g, 
                        label=model_type.capitalize(), alpha=0.6)
        
        plt.title("ΔG Landscape")
        plt.xlabel("Electron Transfer Step")
        plt.ylabel("ΔG (eV)")
        plt.legend()
        
        # Subplot 2: Redox Curves
        plt.subplot(2, 2, 2)
        if 'redox_data' in results:
            redox_data = results['redox_data']
            plt.plot(redox_data['E_range'], redox_data['f_ox']['independent'], 
                     label='Independent', linestyle='--')
            plt.plot(redox_data['E_range'], redox_data['f_ox']['sequential'], 
                     label='Sequential')
            plt.title("Oxidized Fraction")
            plt.xlabel("Potential (V)")
            plt.ylabel("Oxidized Fraction")
            plt.legend()
        
        plt.tight_layout()
        plt.show()
    
    def validate(self):
        """Validate cooperativity analysis"""
        return hasattr(self.ee_gui, 'results') and 'cooperativity' in self.ee_gui.results

class VisualizationFrame(ttk.Frame):
    def __init__(self, parent, ee_gui):
        super().__init__(parent)
        self.ee_gui = ee_gui
        
        # Create main layout
        self.create_summary_section()
        self.create_comparative_visualization()
        self.create_export_options()
    
    def create_summary_section(self):
        """Create summary of all calculation results"""
        summary_frame = ttk.LabelFrame(self, text="Analysis Summary")
        summary_frame.pack(fill='x', padx=10, pady=10)
        
        # Text widget to display summary
        self.summary_text = scrolledtext.ScrolledText(
            summary_frame, 
            wrap=tk.WORD, 
            height=10
        )
        self.summary_text.pack(fill='both', expand=True, padx=5, pady=5)
    
    def create_comparative_visualization(self):
        """Create comparative visualization options"""
        viz_frame = ttk.LabelFrame(self, text="Comparative Visualizations")
        viz_frame.pack(fill='x', padx=10, pady=10)
        
        # Visualization type selection
        self.viz_var = tk.StringVar(value='multi-param')
        
        ttk.Radiobutton(viz_frame, text="Multi-Parameter", 
                        variable=self.viz_var, value='multi-param').pack(side='left', padx=5)
        ttk.Radiobutton(viz_frame, text="Electron Transfer Landscape", 
                        variable=self.viz_var, value='landscape').pack(side='left', padx=5)
        ttk.Radiobutton(viz_frame, text="Redox Behavior", 
                        variable=self.viz_var, value='redox').pack(side='left', padx=5)
        
        # Visualization button
        ttk.Button(viz_frame, text="Generate Visualization", 
                   command=self.generate_visualization).pack(side='right', padx=5)
    
    def create_export_options(self):
        """Create export and reporting options"""
        export_frame = ttk.LabelFrame(self, text="Export and Report")
        export_frame.pack(fill='x', padx=10, pady=10)
        
        # Export format selection
        export_frame_inner = ttk.Frame(export_frame)
        export_frame_inner.pack(side='left', padx=5)
        
        ttk.Label(export_frame_inner, text="Export Format:").pack(side='left')
        self.export_var = tk.StringVar(value='pdf')
        export_dropdown = ttk.Combobox(
            export_frame_inner, 
            textvariable=self.export_var, 
            values=['pdf', 'json', 'csv'], 
            state='readonly', 
            width=10
        )
        export_dropdown.pack(side='left', padx=5)
        
        # Export and report buttons
        ttk.Button(export_frame, text="Generate Report", 
                   command=self.generate_report).pack(side='right', padx=5)
        ttk.Button(export_frame, text="Export Results", 
                   command=self.export_results).pack(side='right', padx=5)
    
    def generate_visualization(self):
        """Generate comparative visualization based on selected type"""
        # Validate all required results are present
        if not self.validate_results():
            messagebox.showerror("Visualization Error", "Please complete all previous calculation steps")
            return
        
        # Select visualization method
        if self.viz_var.get() == 'multi-param':
            self.visualize_multi_parameter()
        elif self.viz_var.get() == 'landscape':
            self.visualize_electron_transfer_landscape()
        elif self.viz_var.get() == 'redox':
            self.visualize_redox_behavior()
    
    def visualize_multi_parameter(self):
        """Create multi-parameter visualization"""
        plt.figure(figsize=(15, 10))
        
        # Lambda values
        plt.subplot(2, 3, 1)
        lambda_values = self.ee_gui.results.get('lambda', {}).get('values', [])
        plt.bar(range(1, len(lambda_values) + 1), lambda_values)
        plt.title("Reorganization Energy")
        plt.xlabel("Step")
        plt.ylabel("λ (meV)")
        
        # Delta G values
        plt.subplot(2, 3, 2)
        delta_g_values = self.ee_gui.results.get('delta_g', {}).get('values', [])
        plt.bar(range(1, len(delta_g_values) + 1), delta_g_values)
        plt.title("Reaction Free Energy")
        plt.xlabel("Step")
        plt.ylabel("ΔG (eV)")
        
        # Coupling values
        plt.subplot(2, 3, 3)
        coupling_values = self.ee_gui.results.get('coupling', [])
        plt.bar(range(1, len(coupling_values) + 1), coupling_values)
        plt.title("Electronic Coupling")
        plt.xlabel("Step")
        plt.ylabel("Coupling (meV)")
        
        # Interaction matrix heatmap
        plt.subplot(2, 3, 4)
        interactions = self.ee_gui.results.get('interactions', {})
        sequence = self.ee_gui.sequence
        matrix = np.zeros((len(sequence), len(sequence)))
        
        for (h1, h2), energy in interactions.items():
            i, j = sequence.index(h1), sequence.index(h2)
            matrix[i, j] = energy
            matrix[j, i] = energy
        
        plt.imshow(matrix, cmap='coolwarm')
        plt.colorbar(label='Interaction Energy (eV)')
        plt.title("Interaction Energy Matrix")
        plt.xlabel("Heme Index")
        plt.ylabel("Heme Index")
        
        plt.tight_layout()
        plt.show()
    
    def visualize_electron_transfer_landscape(self):
        """Visualize electron transfer landscape"""
        # Implement detailed landscape visualization
        pass
    
    def visualize_redox_behavior(self):
        """Visualize redox behavior"""
        cooperativity_results = self.ee_gui.results.get('cooperativity', {})
        
        if 'redox_data' not in cooperativity_results:
            messagebox.showwarning("Visualization Error", "Redox data not available")
            return
        
        redox_data = cooperativity_results['redox_data']
        
        plt.figure(figsize=(12, 6))
        
        # Independent model oxidized fraction
        plt.subplot(1, 2, 1)
        plt.plot(redox_data['E_range'], redox_data['f_ox']['independent'], label='Independent')
        plt.title("Oxidized Fraction")
        plt.xlabel("Potential (V)")
        plt.ylabel("Oxidized Fraction")
        plt.legend()
        
        # Sequential model oxidized fraction
        plt.subplot(1, 2, 2)
        plt.plot(redox_data['E_range'], redox_data['f_ox']['sequential'], label='Sequential')
        plt.title("Sequential Oxidized Fraction")
        plt.xlabel("Potential (V)")
        plt.ylabel("Oxidized Fraction")
        plt.legend()
        
        plt.tight_layout()
        plt.show()
    
    def generate_report(self):
        """Generate comprehensive analysis report"""
        # Implement report generation logic
        pass
    
    def export_results(self):
        """Export results in selected format"""
        export_format = self.export_var.get()
        
        if export_format == 'pdf':
            self.export_pdf()
        elif export_format == 'json':
            self.export_json()
        elif export_format == 'csv':
            self.export_csv()
    
    def validate_results(self):
        """Validate that all necessary results are present"""
        required_results = [
            'lambda', 'delta_g', 'interactions', 
            'coupling', 'cooperativity'
        ]
        
        for result in required_results:
            if result not in self.ee_gui.results:
                return False
        
        return True   

def main():
    # Choose launch directory
    launch_dir = Path.cwd() / "EE_Analysis"
    launch_dir.mkdir(parents=True, exist_ok=True)
    
    # Choose forcefield directory (adjust as needed)
    forcefield_dir = Path("/path/to/your/forcefields")
    
    # Select PDB file
    pdb_file = filedialog.askopenfilename(
        title="Select PDB File", 
        filetypes=[('PDB Files', '*.pdb')]
    )
    
    if not pdb_file:
        print("No PDB file selected. Exiting.")
        return
    
    # Initialize interaction manager
    interaction_manager = InteractionManager(launch_dir=launch_dir)
    
    # Initialize calculators
    lambda_calculator = LambdaCalculator(
        interaction_manager, launch_dir)
    
    delta_g_calculator = DeltaGCalculator(
        interaction_manager, launch_dir, forcefield_dir, pdb_file)
    
    interaction_calculator = HemeInteractionCalculator(
        interaction_manager, launch_dir, forcefield_dir, pdb_file)
    
    coupling_calculator = CouplingCalculator(
        interaction_manager, launch_dir)
    
    cooperativity_analyzer = HemeCooperativityAnalyzer(
        interaction_manager, launch_dir / "EE")
    
    # Create and run the GUI
    ee_gui = EnergeticEvaluationGUI(
        interaction_manager=interaction_manager,
        launch_dir=launch_dir,
        forcefield_dir=forcefield_dir,
        pdb_file=pdb_file
    )
    
    # Run the GUI main loop
    ee_gui.root.mainloop()

if __name__ == "__main__":
    main()






