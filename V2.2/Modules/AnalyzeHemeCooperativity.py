import os
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patheffects import withStroke
from matplotlib.colors import LinearSegmentedColormap, to_rgba
from tabulate import tabulate
from argparse import ArgumentParser

# Constants
R = 8.314  # J/(mol*K)
T = 300    # K (25°C)
F = 96485  # C/mol

def ensure_output_dir(directory):
    """Create output directory if it doesn't exist"""
    output_dir = os.path.join(directory, 'output')
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    return output_dir

def read_energy_matrix(filename, energy_shift=0, interaction_scale=1.0):
    """Read and adjust the energy matrix from a file.
    
    Parameters
    ----------
    filename : str
        Path to the energy matrix file
    energy_shift : float
        Shift to apply to diagonal elements (site energies)
    interaction_scale : float
        Scaling factor for off-diagonal elements (interactions)
        
    Returns
    -------
    tuple
        (original_matrix, adjusted_matrix)
    """
    import numpy as np
    
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    # Find where matrix starts (first line starting with '[')
    matrix_start = next(i for i, line in enumerate(lines) if line.strip().startswith('['))
    
    # Read matrix
    matrix = []
    for line in lines[matrix_start:]:
        if line.strip().startswith('['):
            # Clean up the line and convert to list of floats
            cleaned = line.strip('[] \n').replace(' ', '')
            values = [float(x)/1000 for x in cleaned.split(',')]  # Convert from meV to eV
            matrix.append(values)
    
    original_matrix = np.array(matrix)
    adjusted_matrix = original_matrix.copy()
    
    # Apply energy shift to diagonal elements (site energies)
    np.fill_diagonal(adjusted_matrix, np.diag(adjusted_matrix) + energy_shift)
    
    # Apply interaction scaling to off-diagonal elements
    off_diagonal_mask = ~np.eye(adjusted_matrix.shape[0], dtype=bool)
    adjusted_matrix[off_diagonal_mask] *= interaction_scale
    
    return original_matrix, adjusted_matrix

def calculate_independent_oxidation(energy_matrix):
    """Calculate oxidation energies without interactions"""
    print("\nCalculating independent model (no interactions):")

    # Get site energies from diagonal
    site_energies = np.diagonal(energy_matrix).copy()
    print(f"Site energies: {site_energies}")

    # Create potentials dictionary
    ind_potentials = {i+1: energy for i, energy in enumerate(site_energies)}

    # Create history - each heme keeps its original energy until oxidized
    n = len(site_energies)
    energy_history = []

    # Initial state - all original energies
    current = site_energies.copy()
    energy_history.append(current.copy())

    # For each oxidation step
    order = np.argsort(site_energies)  # Sort by energy
    for heme in order:
        current = site_energies.copy()
        # Set previously oxidized hemes to inf
        current[order[:np.where(order == heme)[0][0]+1]] = float('inf')
        energy_history.append(current)

    print("\nFinal results:")
    print(f"Independent potentials: {ind_potentials}")

    return order, energy_history, ind_potentials

def calculate_heme_oxidation(energy_matrix, interaction_level=None):
    print("\nCalculating heme oxidation sequence:")
    print(f"Input energy matrix:\n{energy_matrix}")
    print(f"Interaction level: {interaction_level}")

    n = len(energy_matrix)
    oxidation_order = []
    current_energies = np.diagonal(energy_matrix).copy()
    energy_history = [current_energies.copy()]

    print(f"\nInitial site energies: {current_energies}")

    independent_potentials = {i+1: energy for i, energy in enumerate(current_energies)}
    sequential_potentials = {}

    for step in range(n):
        print(f"\nStep {step + 1}:")
        heme_to_oxidize = np.argmin(current_energies)
        print(f"  Heme to oxidize: {heme_to_oxidize + 1}")
        print(f"  Energy: {current_energies[heme_to_oxidize]:.3f} eV")

        oxidation_order.append(heme_to_oxidize)
        sequential_potentials[heme_to_oxidize + 1] = current_energies[heme_to_oxidize]

        print("  Updating remaining heme energies:")
        for i in range(n):
            if i not in oxidation_order:
                old_energy = current_energies[i]
                if interaction_level is None or abs(i - heme_to_oxidize) < interaction_level:
                    current_energies[i] += energy_matrix[heme_to_oxidize][i]
                    print(f"    Heme {i + 1}: {old_energy:.3f} eV -> {current_energies[i]:.3f} eV")

        current_energies[heme_to_oxidize] = float('inf')
        energy_history.append(current_energies.copy())

    print("\nFinal results:")
    print(f"Oxidation order: {[x + 1 for x in oxidation_order]}")
    print(f"Independent potentials: {independent_potentials}")
    print(f"Sequential potentials: {sequential_potentials}")

    return oxidation_order, energy_history, independent_potentials, sequential_potentials

def calculate_geometric_oxidation(energy_matrix, interaction_level=None):
    """Calculate oxidation sequence following geometric order of hemes"""
    print("\nCalculating geometric oxidation sequence:")
    print(f"Input energy matrix:\n{energy_matrix}")
    print(f"Interaction level: {interaction_level}")

    n = len(energy_matrix)
    geometric_order = list(range(n))  # 0,1,2,3,4,5,6
    current_energies = np.diagonal(energy_matrix).copy()
    energy_history = [current_energies.copy()]

    print(f"\nInitial site energies: {current_energies}")

    geometric_potentials = {}
    oxidized_hemes = []

    for step, heme in enumerate(geometric_order):
        print(f"\nStep {step + 1}:")
        print(f"  Oxidizing heme: {heme + 1}")
        print(f"  Energy: {current_energies[heme]:.3f} eV")

        # Record the potential at which this heme was oxidized
        geometric_potentials[heme + 1] = current_energies[heme]
        oxidized_hemes.append(heme)

        print("  Updating remaining heme energies:")
        for i in range(n):
            if i not in oxidized_hemes:
                old_energy = current_energies[i]
                if interaction_level is None or abs(i - heme) < interaction_level:
                    current_energies[i] += energy_matrix[heme][i]
                    print(f"    Heme {i + 1}: {old_energy:.3f} eV -> {current_energies[i]:.3f} eV")

        current_energies[heme] = float('inf')
        energy_history.append(current_energies.copy())

    print("\nFinal results:")
    print(f"Geometric oxidation order: {[x + 1 for x in geometric_order]}")
    print(f"Geometric potentials: {geometric_potentials}")

    return geometric_order, energy_history, geometric_potentials

def calculate_thermodynamic_oxidation(energy_matrix, interaction_level=None):
    """Calculate oxidation sequence following thermodynamic favorability"""
    print("\nCalculating thermodynamic oxidation sequence:")
    print(f"Input energy matrix:\n{energy_matrix}")
    print(f"Interaction level: {interaction_level}")

    n = len(energy_matrix)
    oxidation_order = []
    current_energies = np.diagonal(energy_matrix).copy()
    energy_history = [current_energies.copy()]

    print(f"\nInitial site energies: {current_energies}")

    thermodynamic_potentials = {}

    for step in range(n):
        print(f"\nStep {step + 1}:")
        heme_to_oxidize = np.argmin(current_energies)
        print(f"  Heme to oxidize: {heme_to_oxidize + 1}")
        print(f"  Energy: {current_energies[heme_to_oxidize]:.3f} eV")

        oxidation_order.append(heme_to_oxidize)
        thermodynamic_potentials[heme_to_oxidize + 1] = current_energies[heme_to_oxidize]

        print("  Updating remaining heme energies:")
        for i in range(n):
            if i not in oxidation_order:
                old_energy = current_energies[i]
                if interaction_level is None or abs(i - heme_to_oxidize) < interaction_level:
                    current_energies[i] += energy_matrix[heme_to_oxidize][i]
                    print(f"    Heme {i + 1}: {old_energy:.3f} eV -> {current_energies[i]:.3f} eV")

        current_energies[heme_to_oxidize] = float('inf')
        energy_history.append(current_energies.copy())

    print("\nFinal results:")
    print(f"Thermodynamic oxidation order: {[x + 1 for x in oxidation_order]}")
    print(f"Thermodynamic potentials: {thermodynamic_potentials}")

    return oxidation_order, energy_history, thermodynamic_potentials

def analyze_pathway_differences(geometric_order, geometric_potentials,
                             thermo_order, thermo_potentials):
    """Analyze differences between geometric and thermodynamic pathways"""
    print("\nPathway Analysis:")
    print("="*50)

    # Compare oxidation orders
    print("\nOxidation Order Comparison:")
    print(f"Geometric pathway: {[x+1 for x in geometric_order]}")
    print(f"Thermodynamic pathway: {[x+1 for x in thermo_order]}")

    # Find where paths diverge
    first_difference = None
    for i, (g, t) in enumerate(zip(geometric_order, thermo_order)):
        if g != t:
            first_difference = i
            break

    if first_difference is not None:
        print(f"\nPaths first diverge at step {first_difference + 1}")
        print(f"Geometric oxidizes heme {geometric_order[first_difference] + 1}")
        print(f"Thermodynamic oxidizes heme {thermo_order[first_difference] + 1}")
    else:
        print("\nPaths are identical")

    # Compare total energy required
    geo_total = sum(geometric_potentials.values())
    thermo_total = sum(thermo_potentials.values())
    energy_difference = geo_total - thermo_total

    print("\nEnergy Analysis:")
    print(f"Total energy for geometric pathway: {geo_total:.3f} eV")
    print(f"Total energy for thermodynamic pathway: {thermo_total:.3f} eV")
    print(f"Energy penalty for geometric pathway: {energy_difference:.3f} eV")

    # Detailed step-by-step comparison
    print("\nStep-by-step comparison:")
    print("Step  Geometric  Thermodynamic  Difference")
    print("-" * 45)
    for i in range(len(geometric_order)):
        geo_heme = geometric_order[i] + 1
        thermo_heme = thermo_order[i] + 1
        geo_energy = geometric_potentials[geo_heme]
        thermo_energy = thermo_potentials[thermo_heme]
        diff = geo_energy - thermo_energy
        print(f"{i+1:4d}  {geo_heme:4d}({geo_energy:7.3f})  "
              f"{thermo_heme:4d}({thermo_energy:7.3f})  {diff:10.3f}")

    return {
        'divergence_step': first_difference,
        'total_energy_difference': energy_difference,
        'geometric_total': geo_total,
        'thermodynamic_total': thermo_total
    }

def calculate_delta_G(potentials, adjacent_only=True, make_periodic=True):
    print("\nCalculating ΔG values:")
    print(f"Input potentials: {potentials}")
    print(f"Adjacent only mode: {adjacent_only}")
    print(f"Periodic mode: {make_periodic}")

    n = len(potentials)
    delta_G = []

    if adjacent_only:
        print("\nCalculating ΔG for adjacent pairs only:")
        for i in range(n):
            # For last pair, only calculate if make_periodic is True
            if i == n-1 and not make_periodic:
                break
                
            j = (i + 1) % n  # Circular sequence if make_periodic, otherwise stops at n-1
            E_i = potentials[i+1]
            E_j = potentials[j+1]
            dG = -1 * ((-1 * E_i) + E_j)
            print(f"\nPair {i+1} -> {j+1}:")
            print(f"  E{i+1} = {E_i:.3f} eV")
            print(f"  E{j+1} = {E_j:.3f} eV")
            print(f"  ΔG = -1 * ((-1 * {E_i:.3f}) + {E_j:.3f}) = {dG:.3f} eV")
            delta_G.append((i+1, j+1, E_i, E_j, dG))
    else:
        print("\nCalculating ΔG for all possible pairs:")
        for i in range(n):
            for j in range(i+1, n):
                E_i = potentials[i+1]
                E_j = potentials[j+1]
                dG = -1 * ((-1 * E_i) + E_j)
                print(f"\nPair {i+1} -> {j+1}:")
                print(f"  E{i+1} = {E_i:.3f} eV")
                print(f"  E{j+1} = {E_j:.3f} eV")
                print(f"  ΔG = -1 * ((-1 * {E_i:.3f}) + {E_j:.3f}) = {dG:.3f} eV")
                delta_G.append((i+1, j+1, E_i, E_j, dG))

    print("\nFinal ΔG results:")
    for result in delta_G:
        print(f"Heme {result[0]} -> Heme {result[1]}: ΔG = {result[4]:.3f} eV")

    return delta_G

def write_delta_G(filename, delta_G):
    with open(filename, 'w') as f:
        for item in delta_G:
            f.write(f"(HEM-{item[0]} = {item[2]:.3f} eV) -> (HEM-{item[1]} = {item[3]:.3f} eV); DG = {item[4]:.3f} eV\n")

def calculate_K_sequential(E_values):
    K = [1]
    cumulative_sum = 0
    for E in E_values:
        cumulative_sum += E
        K.append(math.exp(-cumulative_sum * F / (R * T)))
    return K

def calculate_f_red_sequential(E, E_values):
    n = len(E_values)
    X = np.exp(E * F / (R * T))
    K = calculate_K_sequential(E_values)
    
    numerator = sum((n-i) * K[i] * X**i for i in range(n+1))
    denominator = sum(K[i] * X**i for i in range(n+1))
    
    return numerator / (n * denominator)

def calculate_f_ox_sequential(E, E_values):
    n = len(E_values)
    X = np.exp(E * F / (R * T))
    K = calculate_K_sequential(E_values)
    
    numerator = sum(i * K[i] * X**i for i in range(n+1))
    denominator = sum(K[i] * X**i for i in range(n+1))
    
    return numerator / (n * denominator)

def calculate_b_independent(E_values):
    return [np.exp(-E_i * F / (R * T)) for E_i in E_values]

def calculate_f_ox_independent(E, E_values):
    n = len(E_values)
    X = np.exp(E * F / (R * T))
    b = calculate_b_independent(E_values)
    
    return (1/n) * sum((b_i * X) / (b_i * X + 1) for b_i in b)

def calculate_f_red_independent(E, E_values):
    n = len(E_values)
    X = np.exp(E * F / (R * T))
    b = calculate_b_independent(E_values)
    
    return (1/n) * sum(1 / (b_i * X + 1) for b_i in b)

def calculate_fractions(E_range, potentials, source_name):
    """Calculate oxidized and reduced fractions for a set of potentials"""
    print(f"\nCalculating fractions for {source_name}")
#   print("Potentials keys:", potentials.keys())
    
    data = {}
    
    if 'independent' in potentials:
        print("Processing independent model")
        ind_values = list(potentials['independent'].values())
        data[f'F_ox ({source_name} Independent)'] = [calculate_f_ox_independent(E, ind_values) for E in E_range]
        data[f'F_red ({source_name} Independent)'] = [calculate_f_red_independent(E, ind_values) for E in E_range]
    
    if 'sequential' in potentials:
        print("Processing sequential model")
        for i, seq_potentials in enumerate(potentials['sequential']):
            seq_values = list(seq_potentials.values())
            data[f'F_ox ({source_name} Sequential {i+1})'] = [calculate_f_ox_sequential(E, seq_values) for E in E_range]
            data[f'F_red ({source_name} Sequential {i+1})'] = [calculate_f_red_sequential(E, seq_values) for E in E_range]
    
#   print("Generated data keys:", data.keys())
    return data

def plot_curves(E_range, all_data, filename, plot_option):
    plt.figure(figsize=(3.3, 3.3), dpi=600)
    
    colors = {'BioDC': 'purple', 'QM': 'green', 'Exp': 'black'}
    
    # Process each source's data
    for source, data in all_data.items():
        base_color = colors[source]
        
        # Only count sequential models for one type (ox) to get correct number
        num_sequential = sum(1 for key in data.keys() 
                           if 'Sequential' in key and 'ox' in key.lower())
        
        # Create color gradients for each source
        if source == 'BioDC':
            sequential_colors = plt.cm.Purples(np.linspace(0.8, 0.2, num_sequential))[::-1]
        elif source == 'QM':
            sequential_colors = plt.cm.Greens(np.linspace(0.8, 0.2, num_sequential))[::-1]
        elif source == 'Exp':
            sequential_colors = plt.cm.Greys(np.linspace(0.8, 0.2, num_sequential))[::-1]
        
        sequential_counter = 0
        for label, values in data.items():
            if (plot_option == 'ox' and 'ox' in label.lower()) or \
               (plot_option == 'red' and 'red' in label.lower()) or \
               plot_option == 'both':
                
                if 'Independent' in label:
                    color = 'white'
                elif 'Sequential' in label:
                    color = sequential_colors[sequential_counter % num_sequential]
                    sequential_counter += 1
                
                linestyle = '--' if 'Independent' in label else '-'
                line, = plt.plot(E_range, values, color=color, linestyle=linestyle)
                line.set_path_effects([withStroke(linewidth=2, foreground=base_color)])
                
                # Only add annotations for BioDC curves
                if source == 'BioDC':
                    if 'Independent' in label:
                        plt.annotate("Ind.", xy=(0.1, 0.15), xycoords='axes fraction',
                                    color='white', fontweight='bold', ha='left', va='center',
                                    bbox=dict(boxstyle='round,pad=0.2', fc='none', ec='none'),
                                    path_effects=[withStroke(linewidth=1, foreground=base_color)],
                                    fontsize=14)
                    
                    if 'Sequential' in label and (sequential_counter % num_sequential) == 0:
                        plt.annotate("Seq.", xy=(0.1, 0.08), xycoords='axes fraction',
                                    color='black', fontweight='bold', ha='left', va='center',
                                    bbox=dict(boxstyle='round,pad=0.2', fc='none', ec='none'),
                                    path_effects=[withStroke(linewidth=1, foreground=base_color)],
                                    fontsize=14)
    
    plt.xlabel('E (V vs. SHE)')
    plt.ylabel('Fraction')
    plt.tick_params(direction='in', which='both', top=True, right=True)
    
    plt.tight_layout()
    plt.savefig(filename, dpi=600, bbox_inches='tight')
    plt.close()

def plot_DG_landscape(directory, delta_G_ind, delta_G_seq_all, filename):
    plt.figure(figsize=(3.3, 3.3), dpi=600)

    num_steps = len(delta_G_ind)
    x = range(1, num_steps + 1)

    # Create color gradient
    num_sequential = len(delta_G_seq_all)
    colors = list(plt.cm.Greys(np.linspace(0.8, 0.2, num_sequential)))  # Light to dark grays

    # Plot sequential models first
    for i, delta_G_seq in enumerate(delta_G_seq_all):
        y_seq = [dg[4] for dg in delta_G_seq]
        plt.plot(x, y_seq, 'k:', linewidth=1, zorder=1)
        plt.scatter(x, y_seq, marker='s', s=50, facecolors=colors[i],
                   edgecolors='black', linewidth=1, zorder=2, label=f'Seq{i+1}')

    # Plot independent model last
    y_ind = [dg[4] for dg in delta_G_ind]
    plt.plot(x, y_ind, 'k:', linewidth=1, zorder=3)
    plt.scatter(x, y_ind, marker='s', s=50, facecolors='white',
               edgecolors='black', linewidth=1, zorder=4, label='Ind.')

    plt.xlabel('Electron Transfer Step')
    plt.ylabel('ΔG (eV)')
    plt.tick_params(direction='in', which='both', top=True, right=True)

    # Set major ticks with interval of 1
    plt.xticks(range(1, num_steps + 1))

    # Place legend inside plot area, in upper right corner
    plt.legend(fontsize='small', ncol=2, columnspacing=0.5, handletextpad=0.1,
              loc='upper right')

    plt.tight_layout()
    plt.savefig(os.path.join(directory, filename), dpi=600, bbox_inches='tight')
    plt.close()

def create_rainbow_colormap(num_hemes):
    colors = ['purple', 'blue', 'green', 'yellow', 'orange', 'red']
    cmap = LinearSegmentedColormap.from_list("rainbow", colors, N=num_hemes)
    return cmap

def get_pale_color(color, factor=0.3):
    """Create a paler version of the given color"""
    rgba = to_rgba(color)
    return (rgba[0] + (1 - rgba[0]) * (1 - factor),
            rgba[1] + (1 - rgba[1]) * (1 - factor),
            rgba[2] + (1 - rgba[2]) * (1 - factor),
            rgba[3])

def check_overlap(pos1, pos2, threshold=0.02):
    """Check if two y-positions are too close"""
    return abs(pos1 - pos2) < threshold

def adjust_label_positions(final_energies):
    """
    Adjust label positions to avoid overlaps by moving overlapping labels in opposite directions
    Returns a list of y-offsets for each label
    """
    # Create list of (energy, index) pairs for sorting
    positions = [(energy, i) for i, energy in enumerate(final_energies)]
    positions.sort()  # Sort by energy

    y_offsets = [0.0] * len(final_energies)  # Start with no offsets

    # Check each pair of adjacent positions
    for i in range(1, len(positions)):
        curr_energy, curr_idx = positions[i]
        prev_energy, prev_idx = positions[i-1]

        if check_overlap(curr_energy, prev_energy):
            # If overlap detected, move both labels apart
            shift = 0.01  # Half of total desired separation
            y_offsets[curr_idx] = shift  # Move higher label up
            y_offsets[prev_idx] = -shift  # Move lower label down

    return y_offsets

def plot_potential_progression(name, energy_history, filename):
    plt.figure(figsize=(3.3, 3.3), dpi=600)

    num_hemes = len(energy_history[0])
    num_stages = len(energy_history)
    x = np.arange(num_stages)

    # Create custom rainbow colormap
    cmap = create_rainbow_colormap(num_hemes)
    colors = [cmap(i / (num_hemes - 1)) for i in range(num_hemes)]

    # Get final non-inf energies for each heme
    final_energies = []
    for heme_num in range(num_hemes):
        last_non_inf = None
        for stage in reversed(range(num_stages)):
            energy = energy_history[stage][heme_num]
            if not np.isinf(energy):
                last_non_inf = energy
                break
        final_energies.append(last_non_inf)

    # Calculate label position adjustments
    y_offsets = adjust_label_positions(final_energies)

    for heme_num in range(num_hemes):
        y = []
        last_non_inf = None

        for stage in range(num_stages):
            energy = energy_history[stage][heme_num]
            if np.isinf(energy):
                if last_non_inf is not None:
                    y.append(last_non_inf)
                else:
                    y.append(np.nan)
            else:
                y.append(energy)
                last_non_inf = energy

        # Plot lines
        plt.plot(x, y, '--', color=colors[heme_num], alpha=0.5)

        # Plot markers
        for i, energy in enumerate(y):
            if np.isnan(energy):
                continue
            if np.isinf(energy_history[i][heme_num]):
                # Unfilled marker: pale heme color fill, black edge
                plt.plot(i, energy, 's', markerfacecolor=get_pale_color(colors[heme_num]),
                         markeredgecolor='black', markersize=6, markeredgewidth=1)
            else:
                # Filled marker: heme color fill, black edge
                plt.plot(i, energy, 's', markerfacecolor=colors[heme_num],
                         markeredgecolor='black', markersize=6, markeredgewidth=1)

        # Add site number next to the last marker with adjusted position
        if last_non_inf is not None:
            label_y = last_non_inf + y_offsets[heme_num]
            plt.text(x[-1] + 0.15, label_y, f'#{heme_num+1}', color=colors[heme_num], fontsize=8,
                     ha='left', va='center', fontweight='normal',
                     path_effects=[withStroke(linewidth=1, foreground='black')])

    plt.xlabel('Oxidation Stage')
    plt.ylabel(r'$\Delta E\ \left(eV\right)$')
    plt.xticks(x, [f'S{i}' for i in range(num_stages)])

    plt.tick_params(direction='in', which='both', top=True, right=True)

    # Extend x-axis to accommodate site numbers
    x_min, x_max = plt.xlim()
    plt.xlim(x_min, x_max + 0.5)

    plt.tight_layout()
    plt.savefig(filename, dpi=600, bbox_inches='tight')
    plt.close()

def process_matrix(name, filepath, model, energy_shift, interaction_scale, adjacent_only, make_periodic, output_dir):
    """Process a single matrix using one of two analysis modes.

    Parameters
    ----------
    name : str
        Name identifier for output files
    filepath : str
        Path to input matrix file
    model : str
        Analysis mode to use ('geo' or 'seq')
    energy_shift : float
        Shift to apply to site energies (diagonal elements)
    interaction_scale : float
        Scaling factor for interaction energies (off-diagonal elements)
    adjacent_only : bool
        If True, only calculate ΔG between adjacent hemes
    make_periodic : bool
        If True, include step from last heme back to first
    output_dir : str
        Directory for output files

    Returns
    -------
    dict
        Results dictionary containing analysis data
    """
    try:
        if model not in ['geo', 'seq']:
            raise ValueError("Model must be either 'geo' or 'seq'")

        if not os.path.exists(filepath):
            print(f"File not found: {filepath}")
            return None

        print(f"\n{'='*50}")
        print(f"Processing matrix for {name}")
        print(f"{'='*50}")
        print(f"Parameters:")
        print(f"  Model: {model}")
        print(f"  Energy shift: {energy_shift}")
        print(f"  Interaction scale: {interaction_scale}")
        print(f"  Adjacent only: {adjacent_only}")
        print(f"  Make periodic: {make_periodic}")

        _, adjusted_matrix = read_energy_matrix(filepath, energy_shift, interaction_scale)
        results = {}

        # Always calculate independent case first
        print("\nCalculating independent model (reference case)...")
        _, ind_history, ind_potentials = calculate_independent_oxidation(adjusted_matrix)
        results['independent'] = ind_potentials
        results['energy_history_ind'] = ind_history

        # Calculate ΔG for independent case
        delta_G_ind = calculate_delta_G(ind_potentials, adjacent_only, make_periodic)
        if not adjacent_only:
            delta_G_ind.sort(key=lambda x: abs(x[4]), reverse=True)
        results['delta_G_ind'] = delta_G_ind

        # Write independent results
        write_delta_G(os.path.join(output_dir, f'DG_ind_{name}.txt'), delta_G_ind)
        plot_potential_progression(name, ind_history,
                                 os.path.join(output_dir, f'potential_progression_ind_{name}.png'))

        if model == 'seq':
            print("\nCalculating sequential model...")
            sequential_potentials_all = []
            delta_G_seq_all = []
            energy_histories = []

            for interaction_level in range(1, len(adjusted_matrix) + 1):
                print(f"\nCalculating with interaction level {interaction_level}...")
                oxidation_order, history, _, seq_potentials = calculate_heme_oxidation(
                    adjusted_matrix, interaction_level)
                sequential_potentials_all.append(seq_potentials)
                energy_histories.append(history)

                delta_G_seq = calculate_delta_G(seq_potentials, adjacent_only, make_periodic)
                if not adjacent_only:
                    delta_G_seq.sort(key=lambda x: abs(x[4]), reverse=True)
                delta_G_seq_all.append(delta_G_seq)

                write_delta_G(os.path.join(output_dir, f'DG_seq_{name}_level_{interaction_level}.txt'),
                             delta_G_seq)
                plot_potential_progression(
                    name, history,
                    os.path.join(output_dir, f'potential_progression_seq_{name}_level_{interaction_level}.png'))

            results['sequential'] = sequential_potentials_all
            results['energy_histories_seq'] = energy_histories
            results['delta_G_seq'] = delta_G_seq_all
            results['energy_history_seq'] = energy_histories[-1]  # Store the final level history
            
            # Generate ΔG landscape plot if adjacent_only is True
            if adjacent_only:
                try:
                    plot_DG_landscape(output_dir, delta_G_ind, delta_G_seq_all,
                                    f'DG_landscape_{name}_seq.png')
                    print("\nΔG landscape plot generated")
                except Exception as e:
                    print(f"Error generating landscape plot: {str(e)}")

#           analyze_potential_changes(energy_histories[-1], "Sequential (highest interaction level)")
            analyze_potential_changes(energy_histories[-1], 
                                    "Sequential (highest interaction level)",
                                    output_dir)
        elif model == 'geo':
            print("\nCalculating geometric model...")
            geometric_order, geo_history, geo_potentials = calculate_geometric_oxidation(adjusted_matrix)
            results['sequential'] = [geo_potentials]  # Keep format compatible with sequential
            results['geometric_order'] = geometric_order
            results['energy_history_geo'] = geo_history

            delta_G_geo = calculate_delta_G(geo_potentials, adjacent_only, make_periodic)
            if not adjacent_only:
                delta_G_geo.sort(key=lambda x: abs(x[4]), reverse=True)
            results['delta_G_geo'] = delta_G_geo

            write_delta_G(os.path.join(output_dir, f'DG_geo_{name}.txt'), delta_G_geo)
            plot_potential_progression(name, geo_history,
                                    os.path.join(output_dir, f'potential_progression_geo_{name}.png'))

            # Generate ΔG landscape plot if adjacent_only is True
            if adjacent_only:
                try:
                    plot_DG_landscape(output_dir, delta_G_ind, [delta_G_geo],
                                    f'DG_landscape_{name}_geo.png')
                    print("\nΔG landscape plot generated")
                except Exception as e:
                    print(f"Error generating landscape plot: {str(e)}")

#           analyze_potential_changes(geo_history, "Geometric")
            analyze_potential_changes(geo_history, "Geometric", output_dir)

        # Save final results for potential progression visualization
        if model == 'seq':
            save_potential_progression(
                f'{name}_sequential_level_{len(adjusted_matrix)}',
                results['energy_history_seq'],
                output_dir)
        elif model == 'geo':
            save_potential_progression(
                f'{name}_geometric',
                results['energy_history_geo'],
                output_dir)

        return results

    except Exception as e:
        print(f"Error processing matrix {name}: {str(e)}")
        raise

def find_global_potential_range(all_potentials):
    """
    Find the global min and max potentials across all models.
    
    Parameters
    ----------
    all_potentials : list
        List of dictionaries containing potential values for different models
        Each dict has 'independent' and/or 'sequential' keys
    
    Returns
    -------
    tuple
        (min_potential - buffer, max_potential + buffer)
    """
    min_potential = float('inf')
    max_potential = float('-inf')
    
    for potentials in all_potentials:
        if 'independent' in potentials:
            values = list(potentials['independent'].values())
            min_potential = min(min_potential, min(values))
            max_potential = max(max_potential, max(values))
        if 'sequential' in potentials:
            for seq_pots in potentials['sequential']:  # seq_pots is already a dict
                values = list(seq_pots.values())
                min_potential = min(min_potential, min(values))
                max_potential = max(max_potential, max(values))
    
    buffer = 0.2
    return min_potential - buffer, max_potential + buffer

def save_data_to_text(filename, E_range, data_dict, E_values_seq_all, E_values_ind):
    with open(filename, 'w') as textfile:
        textfile.write("# Independent potentials: " + " ".join(f"{pot:.3f}" for pot in sorted(E_values_ind.values())) + "\n")
        for i, E_values_seq in enumerate(E_values_seq_all):
            textfile.write(f"# Sequential potentials (level {i+1}): " + " ".join(f"{pot:.3f}" for pot in E_values_seq.values()) + "\n")
        
        headers = ["E"] + list(data_dict.keys())
        col_widths = [max(len(header), 12) for header in headers]
        
        header_line = "  ".join(f"{header:>{width}}" for header, width in zip(headers, col_widths))
        textfile.write(header_line + "\n")
        
        for i, E in enumerate(E_range):
            row_data = [f"{E:.6f}"] + [f"{data_dict[key][i]:.6f}" for key in data_dict.keys()]
            row = "  ".join(f"{data:>{width}}" for data, width in zip(row_data, col_widths))
            textfile.write(row + "\n")

def analyze_matrix(filepath, model, energy_shift=0, interaction_scale=1.0, adjacent_only=True, make_periodic=True):
    """Analyze a single energy matrix with specified parameters"""
    _, adjusted_matrix = read_energy_matrix(filepath, energy_shift, interaction_scale)
    results = {}

    if model in ['ind', 'both']:
        _, energy_history, ind_potentials, _ = calculate_heme_oxidation(adjusted_matrix)
        results['independent'] = ind_potentials
        # Calculate ΔG for independent potentials
        delta_G_ind = calculate_delta_G(ind_potentials, adjacent_only, make_periodic)
        if not adjacent_only:
            delta_G_ind.sort(key=lambda x: abs(x[4]), reverse=True)
        results['delta_G_ind'] = delta_G_ind

    if model in ['seq', 'both']:
        sequential_potentials_all = []
        delta_G_seq_all = []
        energy_histories = []

        for interaction_level in range(1, len(adjusted_matrix) + 1):
            oxidation_order, history, _, seq_potentials = calculate_heme_oxidation(adjusted_matrix, interaction_level)
            sequential_potentials_all.append(seq_potentials)
            energy_histories.append(history)

            # Calculate ΔG for sequential potentials
            delta_G_seq = calculate_delta_G(seq_potentials, adjacent_only, make_periodic)
            if not adjacent_only:
                delta_G_seq.sort(key=lambda x: abs(x[4]), reverse=True)
            delta_G_seq_all.append(delta_G_seq)

        results['sequential'] = sequential_potentials_all
        results['delta_G_seq'] = delta_G_seq_all
        results['energy_histories'] = energy_histories[-1]  # Keep only the final history

    return results, adjusted_matrix

def print_potential_progression(name, energy_history):
    """Print the potential progression for each heme at each oxidation stage.

    Parameters:
    -----------
    name : str
        Name identifier for the analysis
    energy_history : list
        List of numpy arrays containing the potentials for each heme at each stage
    """
    num_hemes = len(energy_history[0])
    num_stages = len(energy_history)

    print(f"\nPotential Progression Data for {name}")
    print("=" * 50)

    # Create headers for the table
    headers = ["Heme"] + [f"Stage {i}" for i in range(num_stages)]

    # Create rows for the table
    table_data = []
    for heme in range(num_hemes):
        row = [f"Heme {heme + 1}"]
        for stage in range(num_stages):
            value = energy_history[stage][heme]
            if np.isinf(value):
                row.append("---")  # Replace inf with dashes for readability
            else:
                row.append(f"{value:.3f}")
        table_data.append(row)

    # Print the table using tabulate
    print(tabulate(table_data, headers=headers, tablefmt="grid"))

    # Print additional statistics
    print("\nPotential Changes Summary:")
    print("-" * 50)
    for heme in range(num_hemes):
        print(f"\nHeme {heme + 1}:")
        last_finite = None
        total_change = 0
        changes = []

        for stage in range(num_stages):
            current = energy_history[stage][heme]
            if not np.isinf(current):
                if last_finite is not None:
                    change = current - last_finite
                    changes.append(change)
                    print(f"  Stage {stage-1} → {stage}: {change:+.3f} eV")
                    total_change += change
                last_finite = current

        if changes:
            print(f"  Total change: {total_change:+.3f} eV")
            print(f"  Average change: {total_change/len(changes):+.3f} eV")

def save_potential_progression_csv(name, energy_history, output_dir):
    """Save potential progression data to CSV file.

    Parameters:
    -----------
    name : str
        Name identifier for the analysis
    energy_history : list
        List of numpy arrays containing the potentials for each heme at each stage
    output_dir : str
        Directory to save the output files
    """
    import csv
    import os

    filename = os.path.join(output_dir, f'potential_progression_{name}.csv')

    num_hemes = len(energy_history[0])
    num_stages = len(energy_history)

    # Create CSV file
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)

        # Write header
        writer.writerow(['Heme'] + [f'Stage_{i}' for i in range(num_stages)])

        # Write data for each heme
        for heme in range(num_hemes):
            row = [f'Heme_{heme + 1}']
            for stage in range(num_stages):
                value = energy_history[stage][heme]
                if np.isinf(value):
                    row.append('NA')  # Use 'NA' for infinite values
                else:
                    row.append(f'{value:.6f}')
            writer.writerow(row)

def save_potential_progression_json(name, energy_history, output_dir):
    """Save potential progression data to JSON file.

    Parameters:
    -----------
    name : str
        Name identifier for the analysis
    energy_history : list
        List of numpy arrays containing the potentials for each heme at each stage
    output_dir : str
        Directory to save the output files
    """
    import json
    import os

    filename = os.path.join(output_dir, f'potential_progression_{name}.json')

    num_hemes = len(energy_history[0])
    num_stages = len(energy_history)

    # Create data structure
    data = {
        'name': name,
        'num_hemes': num_hemes,
        'num_stages': num_stages,
        'hemes': {}
    }

    # Add data for each heme
    for heme in range(num_hemes):
        heme_data = []
        for stage in range(num_stages):
            value = energy_history[stage][heme]
            if np.isinf(value):
                heme_data.append(None)  # Use None for infinite values
            else:
                heme_data.append(float(f'{value:.6f}'))
        data['hemes'][f'Heme_{heme + 1}'] = heme_data

    # Calculate and store changes between stages
    changes = {}
    for heme in range(num_hemes):
        heme_changes = []
        last_finite = None
        for stage in range(num_stages):
            current = energy_history[stage][heme]
            if not np.isinf(current):
                if last_finite is not None:
                    change = float(f'{(current - last_finite):.6f}')
                    heme_changes.append(change)
                last_finite = current
        changes[f'Heme_{heme + 1}'] = heme_changes
    data['changes'] = changes

    # Save to JSON file
    with open(filename, 'w') as jsonfile:
        json.dump(data, jsonfile, indent=2)

def save_potential_progression(name, energy_history, output_dir):
    """Save potential progression data in both CSV and JSON formats.

    Parameters:
    -----------
    name : str
        Name identifier for the analysis
    energy_history : list
        List of numpy arrays containing the potentials for each heme at each stage
    output_dir : str
        Directory to save the output files
    """
    save_potential_progression_csv(name, energy_history, output_dir)
    save_potential_progression_json(name, energy_history, output_dir)

def create_comparison_plot(sequential_results, geometric_results, directory):
    """Create combined analysis plot comparing sequential and geometric pathways."""
    try:
        # Create figure
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(7, 7))

        # Get number of hemes from the matrix size
        num_hemes = len(sequential_results['energy_history_ind'][0])
        heme_positions = range(1, num_hemes + 1)

        # Get initial (independent) potentials
        initial_potentials = sequential_results['energy_history_ind'][0]

        # Safely calculate y-axis limits for top panels
        def get_max_finite_value(history):
            max_values = []
            for stage in history:
                finite_values = [e for e in stage if not np.isinf(e)]
                if finite_values:
                    max_values.append(max(finite_values))
            return max(max_values) if max_values else max(initial_potentials)

        y_min = min(initial_potentials) - abs(min(initial_potentials) * 0.1)
        max_seq = get_max_finite_value(sequential_results['energy_history_seq'])
        max_geo = get_max_finite_value(geometric_results['energy_history_geo'])
        y_max = max(max_seq, max_geo) + abs(max(max_seq, max_geo) * 0.1)

        # Add panel labels and configure ticks
        for ax, label in zip([ax1, ax2, ax3, ax4], ['A', 'B', 'C', 'D']):
            ax.text(-0.2, 1.05, f'({label})', transform=ax.transAxes,
                    fontsize=10, fontweight='bold')
            ax.tick_params(direction='in', which='both', top=True, right=True)

        # Set up top panels (potential progression)
        for ax in (ax1, ax2):
            ax.grid(True, linestyle='--', alpha=0.3)
            ax.set_xlabel('Heme Position')
            ax.set_ylim(y_min, y_max)

        ax1.set_ylabel('Potential (eV)')
        ax2.tick_params(axis='y', which='both', labelleft=False)

        # Add titles
        ax1.text(0.98, 0.94, 'Thermodyn. Pathway',
                 transform=ax1.transAxes, fontsize=10, va='top', ha='right',
                 fontweight='bold')
        ax2.text(0.98, 0.94, 'Geometric Pathway',
                 transform=ax2.transAxes, fontsize=10, va='top', ha='right',
                 fontweight='bold')

        # Create blue gradient for oxidation order
        blue_colors = plt.cm.Blues(np.linspace(0.3, 0.9, num_hemes))

        # Plot initial potentials on top panels
        for ax in (ax1, ax2):
            ax.scatter(heme_positions, initial_potentials, color='gray',
                      edgecolor='black', linewidth=1, s=50, zorder=2)

        # Process Sequential Data (Top Left Panel)
        seq_history = sequential_results['energy_history_seq']
        seq_order = []
        for stage in range(len(seq_history)):
            energies = seq_history[stage]
            for heme_num, energy in enumerate(energies):
                if heme_num not in seq_order and np.isinf(energy):
                    seq_order.append(heme_num)
                    break

        # Plot sequence for sequential pathway
        for i, heme_num in enumerate(seq_order):
            stages = [stage[heme_num] for stage in seq_history]
            # Find the last finite value and its position
            last_finite_idx = -1
            last_finite_val = None
            for j, val in enumerate(stages):
                if not np.isinf(val):
                    last_finite_idx = j
                    last_finite_val = val
                else:
                    break

            if last_finite_val is not None:
                # Plot endpoint with number label
                ax1.scatter([heme_num + 1], [last_finite_val], color=blue_colors[i],
                           edgecolor='black', linewidth=1, s=50, zorder=4)
                text = ax1.text(heme_num + 1.2, last_finite_val, f'{i+1}',
                               color=blue_colors[i], fontsize=10,
                               fontweight='bold', va='center')
                text.set_path_effects([withStroke(linewidth=1.1, foreground='black')])

                # Draw arrow if there's a change from initial
                initial_val = stages[0]
                if abs(last_finite_val - initial_val) > 0.001:
                    ax1.annotate('',
                        xy=(heme_num + 1, last_finite_val),
                        xytext=(heme_num + 1, initial_val),
                        arrowprops=dict(arrowstyle='->', color=blue_colors[i],
                                      linewidth=2, mutation_scale=15),
                        zorder=3
                    )

        # Process Geometric Data (Top Right Panel)
        geo_history = geometric_results['energy_history_geo']
        geo_order = geometric_results['geometric_order']

        for i, heme_num in enumerate(geo_order):
            stages = [stage[heme_num] for stage in geo_history]
            # Find the last finite value and its position
            last_finite_idx = -1
            last_finite_val = None
            for j, val in enumerate(stages):
                if not np.isinf(val):
                    last_finite_idx = j
                    last_finite_val = val
                else:
                    break

            if last_finite_val is not None:
                # Plot endpoint with number label
                ax2.scatter([heme_num + 1], [last_finite_val], color=blue_colors[i],
                           edgecolor='black', linewidth=1, s=50, zorder=4)
                text = ax2.text(heme_num + 1.2, last_finite_val, f'{i+1}',
                               color=blue_colors[i], fontsize=10,
                               fontweight='bold', va='center')
                text.set_path_effects([withStroke(linewidth=1, foreground='black')])

                # Draw arrow if there's a change from initial
                initial_val = stages[0]
                if abs(last_finite_val - initial_val) > 0.001:
                    ax2.annotate('',
                        xy=(heme_num + 1, last_finite_val),
                        xytext=(heme_num + 1, initial_val),
                        arrowprops=dict(arrowstyle='->', color=blue_colors[i],
                                      linewidth=2, mutation_scale=15),
                        zorder=3
                    )

        # Adjust top panel x limits
        for ax in (ax1, ax2):
            ax.set_xlim(0.5, num_hemes + 0.5)
            ax.set_xticks(heme_positions)

        # Set up bottom panels (ΔG landscapes)
        dg_ind = [dg[4] for dg in sequential_results['delta_G_ind']]
        dg_seq = [dg[4] for dg in sequential_results['delta_G_seq'][-1]]
        dg_geo = [dg[4] for dg in geometric_results['delta_G_geo']]
        et_steps = range(1, len(dg_ind) + 1)

        for ax in (ax3, ax4):
            ax.grid(True, linestyle='--', alpha=0.3)
            ax.set_xlabel('Electron Transfer Step')

        ax3.set_ylabel('ΔG (eV)')
        ax4.tick_params(axis='y', which='both', labelleft=False)

        # Calculate y-axis limits for bottom panels with buffer
        dg_min = min(min(dg_ind), min(dg_seq), min(dg_geo))
        dg_max = max(max(dg_ind), max(dg_seq), max(dg_geo))
        dg_range = dg_max - dg_min
        y_buffer = dg_range * 0.1

        for ax in (ax3, ax4):
            ax.set_ylim(dg_min - y_buffer, dg_max + y_buffer)

        # Plot ΔG landscapes
        # Independent (both panels)
        for ax in (ax3, ax4):
            ax.plot(et_steps, dg_ind, 'k:', zorder=1)
            ax.scatter(et_steps, dg_ind, marker='s', s=100, facecolor='white',
                      edgecolor='black', zorder=2)

        # Sequential (left panel)
        ax3.plot(et_steps, dg_seq, 'k:', zorder=1)
        ax3.scatter(et_steps, dg_seq, marker='s', s=100, color='black', zorder=2)

        # Geometric (right panel)
        ax4.plot(et_steps, dg_geo, 'k:', zorder=1)
        ax4.scatter(et_steps, dg_geo, marker='s', s=100, color='black', zorder=2)

        # Add styled annotations for bottom panels
        def style_two_line_text(ax, line1, x, y):
            ax.text(x, y, line1, transform=ax.transAxes,
                   fontsize=10, va='top', ha='right', fontweight='bold')
            ind_text = ax.text(x, y-0.05, 'Independent',
                             transform=ax.transAxes, fontsize=10,
                             va='top', ha='right', fontweight='bold',
                             color='white')
            ind_text.set_path_effects([withStroke(linewidth=1, foreground='black')])

        style_two_line_text(ax3, 'Thermodyn. vs.', 0.95, 0.95)
        style_two_line_text(ax4, 'Geometric vs.', 0.95, 0.95)

        plt.tight_layout()
        output_path = os.path.join(directory, 'oxidation_pathway_comparison.png')
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()

        print(f"\nComparison plot saved to: {output_path}")

#       print("\nAnalyzing potential changes in each pathway...")
#       analyze_potential_changes(sequential_results['energy_history_seq'], "Sequential")
#       analyze_potential_changes(geometric_results['energy_history_geo'], "Geometric")

    except Exception as e:
        print(f"Error in create_comparison_plot: {str(e)}")
        print("Sequential results keys:", sequential_results.keys())
        print("Geometric results keys:", geometric_results.keys())
        raise

def analyze_potential_changes(energy_history, name="", output_dir=None):
    """Analyze how each heme's potential changes through the oxidation sequence."""
    
    # Prepare output
    output_lines = []
    def add_line(line=""):
        print(line)
        output_lines.append(line)

    add_line(f"\n{name} Pathway Analysis:")
    add_line("="*50)

    num_hemes = len(energy_history[0])

    # Determine oxidation order
    oxidation_order = []
    for stage in range(1, len(energy_history)):
        for heme in range(num_hemes):
            if not np.isinf(energy_history[stage-1][heme]) and np.isinf(energy_history[stage][heme]):
                oxidation_order.append(heme)
                break

    for heme in range(num_hemes):
        initial_potential = energy_history[0][heme]
        add_line(f"\nHeme {heme + 1}:")
        add_line(f"Initial potential: {initial_potential:.3f} eV")

        last_finite = initial_potential
        total_change = 0
        changes = []

        for stage in range(1, len(energy_history)):
            current = energy_history[stage][heme]
            if not np.isinf(current):
                change = current - last_finite
                changes.append(change)
                oxidized_heme = oxidation_order[stage-1] + 1
                add_line(f"Step {stage}: {change:+.3f} eV (When Heme {oxidized_heme} is oxidized)")
                total_change += change
                last_finite = current
            elif np.isinf(current):
                if changes:
                    add_line(f"Cumulative shift before oxidation: {total_change:+.3f} eV")
                add_line(f"Oxidized at step {stage}")
                break

    add_line("\nAnalysis complete.")

    # Save to file if output_dir is provided
    if output_dir:
        clean_name = name.lower().replace(" ", "_") \
                        .replace("(", "") \
                        .replace(")", "") \
                        .replace(",", "")
        filename = os.path.join(output_dir, f'potential_changes_{clean_name}.txt')

        with open(filename, 'w') as f:
            for line in output_lines:
                f.write(line + '\n')
        print(f"\nAnalysis saved to: {filename}")

def find_midpoint_potential(E_range, f_red):
    """Find the potential where reduced fraction equals 0.5.
    
    Parameters:
    -----------
    E_range : array-like
        Array of potential values
    f_red : array-like
        Array of reduced fractions corresponding to E_range
        
    Returns:
    --------
    float
        Potential value where f_red = 0.5
    """
    # Find the index where f_red is closest to 0.5
    idx = np.argmin(np.abs(np.array(f_red) - 0.5))
    
    # If we're exactly at 0.5, return that potential
    if abs(f_red[idx] - 0.5) < 1e-6:
        return E_range[idx]
    
    # Otherwise, interpolate between closest points
    if idx > 0 and idx < len(E_range) - 1:
        # Find points on either side of 0.5
        if f_red[idx] > 0.5:
            idx_before, idx_after = idx, idx + 1
        else:
            idx_before, idx_after = idx - 1, idx
            
        # Linear interpolation
        E1, E2 = E_range[idx_before], E_range[idx_after]
        f1, f2 = f_red[idx_before], f_red[idx_after]
        
        return E1 + (E2 - E1) * (0.5 - f1) / (f2 - f1)
    
    return E_range[idx]

def analyze_midpoint_potentials(E_range, potentials, source_name):
    """Analyze midpoint potentials for all models.
    
    Parameters:
    -----------
    E_range : array-like
        Array of potential values
    potentials : dict
        Dictionary containing potential values for different models
    source_name : str
        Name identifier for the data source
        
    Returns:
    --------
    dict
        Dictionary containing midpoint potentials for each model
    """
    results = {}
    
    if 'independent' in potentials:
        print(f"\nAnalyzing midpoint potentials for {source_name} (Independent model)")
        ind_values = list(potentials['independent'].values())
        f_red = [calculate_f_red_independent(E, ind_values) for E in E_range]
        E_half = find_midpoint_potential(E_range, f_red)
        results['independent'] = E_half
        print(f"Independent model E1/2: {E_half:.3f} V")
    
    if 'sequential' in potentials:
        results['sequential'] = []
        for i, seq_potentials in enumerate(potentials['sequential']):
            print(f"\nAnalyzing {source_name} Sequential model (Level {i+1})")
            seq_values = list(seq_potentials.values())
            f_red = [calculate_f_red_sequential(E, seq_values) for E in E_range]
            E_half = find_midpoint_potential(E_range, f_red)
            results['sequential'].append(E_half)
            print(f"Sequential model (Level {i+1}) E1/2: {E_half:.3f} V")
    
    return results

def main():
    """Main function for heme cooperativity analysis."""
    parser = ArgumentParser(description='Analyze heme cooperativity')
    parser.add_argument('plot_option', help='Plot option (ox, red, or both)')
    parser.add_argument('adjacent_only', type=str, help='Calculate only adjacent hemes (true/false)')
    parser.add_argument('make_periodic', type=str, help='Include step from last heme back to first (true/false)')
    parser.add_argument('biodc_mat', help='BioDC matrix file')
    parser.add_argument('biodc_model', help='BioDC model type (geo, seq, or both)')
    parser.add_argument('biodc_eng_shift', type=float, help='BioDC energy shift')
    parser.add_argument('biodc_int_scale', type=float, help='BioDC interaction scale')
    parser.add_argument('qm_mat', nargs='?', help='QM matrix file')
    parser.add_argument('qm_model', nargs='?', help='QM model type')
    parser.add_argument('qm_eng_shift', type=float, nargs='?', help='QM energy shift')
    parser.add_argument('qm_int_scale', type=float, nargs='?', help='QM interaction scale')
    parser.add_argument('exp_mat', nargs='?', help='Experimental matrix file')
    parser.add_argument('exp_model', nargs='?', help='Experimental model type')
    parser.add_argument('exp_eng_shift', type=float, nargs='?', help='Experimental energy shift')
    parser.add_argument('exp_int_scale', type=float, nargs='?', help='Experimental interaction scale')

    args = parser.parse_args()

    # Convert arguments
    adjacent_only = args.adjacent_only.lower() == 'true'
    make_periodic = args.make_periodic.lower() == 'true'

    # Create output directory
    output_dir = ensure_output_dir('.')

    print("\nDebug - Command line arguments:")
    print(f"biodc_mat: {args.biodc_mat}")
    print(f"biodc_model: {args.biodc_model}")
    print(f"biodc_eng_shift: {args.biodc_eng_shift}")
    print(f"biodc_int_scale: {args.biodc_int_scale}")
    print(f"qm_mat: {args.qm_mat}")
    print(f"qm_model: {args.qm_model}")
    print(f"qm_eng_shift: {args.qm_eng_shift}")
    print(f"qm_int_scale: {args.qm_int_scale}")
    print(f"exp_mat: {args.exp_mat}")
    print(f"exp_model: {args.exp_model}")
    print(f"exp_eng_shift: {args.exp_eng_shift}")
    print(f"exp_int_scale: {args.exp_int_scale}")

    # Process matrices
    all_results = {}
    source_names = ['BioDC', 'QM', 'Exp']
    matrices_info = [
        (args.biodc_mat, args.biodc_model, args.biodc_eng_shift, args.biodc_int_scale),
        (args.qm_mat, args.qm_model, args.qm_eng_shift, args.qm_int_scale),
        (args.exp_mat, args.exp_model, args.exp_eng_shift, args.exp_int_scale)
    ]

    # Filter out incomplete matrix entries
    matrices = []
    valid_names = []
    for i, (mat, model, shift, scale) in enumerate(matrices_info):
        if all(x is not None for x in (mat, model, shift, scale)):
            matrices.append((mat, model, shift, scale))
            valid_names.append(source_names[i])

    print(f"\nProcessing {len(matrices)} matrices with names: {valid_names}")

    for (filepath, model, energy_shift, interaction_scale), name in zip(matrices, valid_names):
        print(f"\nProcessing {name} matrix:")
        print(f"  File: {filepath}")
        print(f"  Model: {model}")
        print(f"  Energy shift: {energy_shift}")
        print(f"  Interaction scale: {interaction_scale}")

        try:
            if model.lower() == 'both':
                # Process both models
                seq_result = process_matrix(name, filepath, 'seq', energy_shift,
                                         interaction_scale, adjacent_only, make_periodic, output_dir)
                geo_result = process_matrix(name, filepath, 'geo', energy_shift,
                                         interaction_scale, adjacent_only, make_periodic, output_dir)

                if seq_result and geo_result:
                    # Create comparison plot
                    create_comparison_plot(seq_result, geo_result, output_dir)

                    # Store results
                    all_results[f"{name}_seq"] = seq_result
                    all_results[f"{name}_geo"] = geo_result
            else:
                # Process single model
                result = process_matrix(name, filepath, model, energy_shift, interaction_scale,
                                      adjacent_only, make_periodic, output_dir)
                if result:
                    all_results[name] = result

        except Exception as e:
            print(f"Error processing {name} matrix: {str(e)}")
            import traceback
            print(traceback.format_exc())
            continue

    # Find global potential range for redox plots
    all_potentials = []
    for result in all_results.values():
        if 'independent' in result:
            all_potentials.append({'independent': result['independent']})
        if 'sequential' in result:
            all_potentials.append({'sequential': result['sequential']})

    lower_bound, upper_bound = find_global_potential_range(all_potentials)
    E_range = np.linspace(lower_bound, upper_bound, 1600)

    print("\nAnalyzing Midpoint Potentials")
    print("=" * 50)
    
    midpoint_results = {}
    for name, result in all_results.items():
        if result:
            midpoint_results[name] = analyze_midpoint_potentials(E_range, result, name)
    
    # Save midpoint potential results to file
    output_file = os.path.join(output_dir, 'midpoint_potentials.txt')
    with open(output_file, 'w') as f:
        f.write("Midpoint Potential Analysis Results\n")
        f.write("=" * 50 + "\n\n")
        
        for source, results in midpoint_results.items():
            f.write(f"\n{source}:\n")
            f.write("-" * len(source) + "\n")
            
            if 'independent' in results:
                f.write(f"Independent model E1/2: {results['independent']:.3f} V vs. SHE\n")
            
            if 'sequential' in results:
                for i, E_half in enumerate(results['sequential']):
                    f.write(f"Sequential model (Level {i+1}) E1/2: {E_half:.3f} V vs. SHE\n")
            
            f.write("\n")
    
    print(f"\nMidpoint potential analysis results saved to: {output_file}")

    # Calculate fractions and create plots for each model type
    if args.biodc_model.lower() == 'both':
        # Handle both models case
        for model_type in ['seq', 'geo']:
            try:
                print(f"\nProcessing {model_type} model...")
                print(f"All results keys: {all_results.keys()}")
                
                # Get data for all sources for this model type
                plot_data = {}
                for name, result in all_results.items():
                    print(f"Processing result for {name}")
                    # Handle BioDC's special case (ends with _seq or _geo)
                    if name.endswith(f"_{model_type}"):
                        source_name = name.replace(f"_{model_type}", "")
                        print(f"Adding {source_name} to plot_data")
                        plot_data[source_name] = calculate_fractions(E_range, result, source_name)
                    # Handle QM and Exp (no suffix)
                    elif name in ['QM', 'Exp']:
                        print(f"Adding {name} to plot_data")
                        plot_data[name] = calculate_fractions(E_range, result, name)
                
                print(f"Final plot_data keys: {plot_data.keys()}")
                
                if plot_data:  # Only proceed if we have data
                    # Save data for each source
                    for source_name, source_data in plot_data.items():
                        output_data_file = os.path.join(output_dir, f'redox_data_{source_name}_{model_type}.txt')
                        save_data_to_text(output_data_file, E_range, source_data, [], {})
                        print(f"Saved data for {source_name} to {output_data_file}")
                    
                    # Create plot with all sources
                    output_plot_file = os.path.join(output_dir, f'redox_plot_{model_type}.png')
                    plot_curves(E_range, plot_data, output_plot_file, args.plot_option)
                    print(f"Saved plot to {output_plot_file}")
                else:
                    print(f"No data found for {model_type} model")

            except Exception as e:
                print(f"Error processing {model_type} model: {str(e)}")
                import traceback
                print(traceback.format_exc())
                continue

    else:
        # Handle single model case - use all results
        try:
            print(f"\nProcessing single model: {args.biodc_model}")
            plot_data = {}
            for name, result in all_results.items():
                # Create nested structure required by plot_curves
                plot_data[name] = calculate_fractions(E_range, result, name)
            
            # Save data for each source
            for source_name, source_data in plot_data.items():
                output_data_file = os.path.join(output_dir, f'redox_data_{source_name}_{args.biodc_model}.txt')
                save_data_to_text(output_data_file, E_range, source_data, [], {})
                print(f"Saved data for {source_name} to {output_data_file}")
            
            output_plot_file = os.path.join(output_dir, f'redox_plot_{args.biodc_model}.png')
            plot_curves(E_range, plot_data, output_plot_file, args.plot_option)
            print(f"Saved plot to {output_plot_file}")

        except Exception as e:
            print(f"Error processing single model: {str(e)}")
            import traceback
            print(traceback.format_exc())

if __name__ == "__main__":
    sys.exit(main())
