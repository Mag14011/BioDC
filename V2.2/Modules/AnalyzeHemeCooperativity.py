import os
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patheffects import withStroke
from matplotlib.colors import LinearSegmentedColormap, to_rgba
from tabulate import tabulate

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
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    matrix_start = next(i for i, line in enumerate(lines) if line.strip().startswith('['))
    matrix = []
    for line in lines[matrix_start:]:
        if line.strip().startswith('['):
            row = [float(x)/1000 for x in line.strip()[1:-1].split(',')]
            matrix.append(row)
    
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

def calculate_delta_G(potentials, adjacent_only=True):
    print("\nCalculating ΔG values:")
    print(f"Input potentials: {potentials}")
    print(f"Adjacent only mode: {adjacent_only}")

    n = len(potentials)
    delta_G = []

    if adjacent_only:
        print("\nCalculating ΔG for adjacent pairs only:")
        for i in range(n):
            j = (i + 1) % n  # Circular sequence
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
    data = {}
    
    if 'independent' in potentials:
        ind_values = list(potentials['independent'].values())
        data[f'F_ox ({source_name} Independent)'] = [calculate_f_ox_independent(E, ind_values) for E in E_range]
        data[f'F_red ({source_name} Independent)'] = [calculate_f_red_independent(E, ind_values) for E in E_range]
    
    if 'sequential' in potentials:
        for i, seq_potentials in enumerate(potentials['sequential']):
            seq_values = list(seq_potentials.values())
            data[f'F_ox ({source_name} Sequential {i+1})'] = [calculate_f_ox_sequential(E, seq_values) for E in E_range]
            data[f'F_red ({source_name} Sequential {i+1})'] = [calculate_f_red_sequential(E, seq_values) for E in E_range]
    
    return data

def plot_curves(E_range, all_data, filename, plot_option):
    plt.figure(figsize=(3.3, 3.3), dpi=600)
    
    colors = {'BioDC': 'black', 'QM': 'red', 'Exp': 'blue'}
    
    # Process each source's data
    for source, data in all_data.items():
        base_color = colors[source]
        
        # Only count sequential models for one type (ox) to get correct number
        num_sequential = sum(1 for key in data.keys() 
                           if 'Sequential' in key and 'ox' in key.lower())
        grays = plt.cm.Greys(np.linspace(0.8, 0.2, num_sequential))[::-1]  # Reversed color gradient
        
        sequential_counter = 0
        for label, values in data.items():
            if (plot_option == 'ox' and 'ox' in label.lower()) or \
               (plot_option == 'red' and 'red' in label.lower()) or \
               plot_option == 'both':
                if 'Independent' in label:
                    color = 'white'
                elif 'Sequential' in label:
                    color = grays[sequential_counter % num_sequential]  # Use modulo to handle both ox and red
                    sequential_counter += 1
                
                linestyle = '--' if 'Independent' in label else '-'
                line, = plt.plot(E_range, values, color=color, linestyle=linestyle)
                line.set_path_effects([withStroke(linewidth=2, foreground=base_color)])
                
                # Only add annotations for BioDC curves
                if source == 'BioDC':
                    if 'Independent' in label:
                        mid_index = np.argmin(np.abs(np.array(values) - 0.5))
                        x_mid = E_range[mid_index]
                        y_mid = values[mid_index]
                        
                        plt.annotate("Ind.", xy=(x_mid - 0.08, y_mid), xycoords='data',
                                    color='white', fontweight='bold', ha='center', va='center',
                                    bbox=dict(boxstyle='round,pad=0.2', fc='none', ec='none'),
                                    path_effects=[withStroke(linewidth=1, foreground=base_color)],
                                    fontsize=14)
                    
                    if 'Sequential' in label and (sequential_counter % num_sequential) == 0:
                        mid_index = np.argmin(np.abs(np.array(values) - 0.5))
                        x_mid = E_range[mid_index]
                        y_mid = values[mid_index]
                        
                        plt.annotate("Seq.", xy=(x_mid + 0.1, y_mid), xycoords='data',
                                    color='black', fontweight='bold', ha='center', va='center',
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

def process_matrix(name, filepath, model, energy_shift, interaction_scale, adjacent_only, output_dir):
    """Process a single matrix using one of two analysis modes.

    Parameters
    ----------
    name : str
        Name identifier for output files
    filepath : str
        Path to input matrix file
    model : str
        Analysis mode to use:
        - 'geo': Compare independent model to geometric sequence (structural order)
        - 'seq': Compare independent model to sequential model (varying interaction levels)
    energy_shift : float
        Shift to apply to site energies (diagonal elements)
    interaction_scale : float
        Scaling factor for interaction energies (off-diagonal elements)
    adjacent_only : bool
        If True, only calculate ΔG between adjacent hemes
    output_dir : str
        Directory for output files

    Returns
    -------
    dict
        Results dictionary containing analysis data for both independent and chosen model
    """
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

    _, adjusted_matrix = read_energy_matrix(filepath, energy_shift, interaction_scale)
    results = {}

    # Always calculate independent case first
    print("\nCalculating independent model (reference case)...")
    _, ind_history, ind_potentials = calculate_independent_oxidation(adjusted_matrix)
    results['independent'] = ind_potentials

    # Calculate ΔG for independent case
    delta_G_ind = calculate_delta_G(ind_potentials, adjacent_only)
    if not adjacent_only:
        delta_G_ind.sort(key=lambda x: abs(x[4]), reverse=True)
    results['delta_G_ind'] = delta_G_ind
    results['energy_history_ind'] = ind_history

    # Write independent results
    write_delta_G(os.path.join(output_dir, f'DG_ind_{name}.txt'), delta_G_ind)
    plot_potential_progression(name, ind_history,
                             os.path.join(output_dir, f'potential_progression_ind_{name}.png'))

    # Geometric sequence calculations
    if model == 'geo':
        print("\nCalculating geometric model (structural sequence)...")
        geometric_order, geo_history, geo_potentials = calculate_geometric_oxidation(adjusted_matrix)

        # Store results in format compatible with find_global_potential_range
        results['sequential'] = [geo_potentials]
        results['geometric_order'] = geometric_order
        results['energy_history_geo'] = geo_history

        # Calculate ΔG for geometric case
        delta_G_geo = calculate_delta_G(geo_potentials, adjacent_only)
        if not adjacent_only:
            delta_G_geo.sort(key=lambda x: abs(x[4]), reverse=True)
        results['delta_G_geo'] = delta_G_geo

        # Write geometric results
        write_delta_G(os.path.join(output_dir, f'DG_geo_{name}.txt'), delta_G_geo)

        # Generate plots comparing independent vs geometric
        if adjacent_only:
            try:
                plot_DG_landscape(output_dir, delta_G_ind, [delta_G_geo],
                                f'DG_landscape_{name}.png')
                print("\nΔG landscape plot generated")
            except Exception as e:
                print(f"Error generating landscape plot: {str(e)}")

        plot_potential_progression(name, geo_history,
                                 os.path.join(output_dir, f'potential_progression_geo_{name}.png'))

    # Sequential model calculations
    elif model == 'seq':
        print("\nCalculating sequential model (varying interaction levels)...")
        sequential_potentials_all = []
        delta_G_seq_all = []
        energy_histories = []

        for interaction_level in range(1, len(adjusted_matrix) + 1):
            print(f"\nCalculating with interaction level {interaction_level}...")
            oxidation_order, history, _, seq_potentials = calculate_heme_oxidation(
                adjusted_matrix, interaction_level)
            sequential_potentials_all.append(seq_potentials)
            energy_histories.append(history)

            delta_G_seq = calculate_delta_G(seq_potentials, adjacent_only)
            if not adjacent_only:
                delta_G_seq.sort(key=lambda x: abs(x[4]), reverse=True)
            delta_G_seq_all.append(delta_G_seq)

            # Write sequential results for each level
            write_delta_G(os.path.join(output_dir, f'DG_seq_{name}_level_{interaction_level}.txt'),
                         delta_G_seq)

        # Store results in format compatible with find_global_potential_range
        results['sequential'] = sequential_potentials_all
        results['energy_histories_seq'] = energy_histories
        results['delta_G_seq'] = delta_G_seq_all

        # Generate plots comparing independent vs sequential
        if adjacent_only:
            try:
                plot_DG_landscape(output_dir, delta_G_ind, delta_G_seq_all,
                                f'DG_landscape_{name}.png')
                print("\nΔG landscape plot generated")
            except Exception as e:
                print(f"Error generating landscape plot: {str(e)}")

        plot_potential_progression(name, energy_histories[-1],
                                 os.path.join(output_dir, f'potential_progression_seq_{name}.png'))

    return results

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

def analyze_matrix(filepath, model, energy_shift=0, interaction_scale=1.0, adjacent_only=True):
    """Analyze a single energy matrix with specified parameters"""
    _, adjusted_matrix = read_energy_matrix(filepath, energy_shift, interaction_scale)
    results = {}

    if model in ['ind', 'both']:
        _, energy_history, ind_potentials, _ = calculate_heme_oxidation(adjusted_matrix)
        results['independent'] = ind_potentials
        # Calculate ΔG for independent potentials
        delta_G_ind = calculate_delta_G(ind_potentials, adjacent_only)
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
            delta_G_seq = calculate_delta_G(seq_potentials, adjacent_only)
            if not adjacent_only:
                delta_G_seq.sort(key=lambda x: abs(x[4]), reverse=True)
            delta_G_seq_all.append(delta_G_seq)

        results['sequential'] = sequential_potentials_all
        results['delta_G_seq'] = delta_G_seq_all
        results['energy_histories'] = energy_histories[-1]  # Keep only the final history

    return results, adjusted_matrix

if __name__ == "__main__":
    import sys

    # Validate command line arguments
    if len(sys.argv) < 7:
        print("Usage: script.py <plot_option> <adjacent_only> "
              "<biodc_mat> <biodc_model> <biodc_eng_shift> <biodc_int_scale> "
              "[<qm_mat> <qm_model> <qm_eng_shift> <qm_int_scale> "
              "[<exp_mat> <exp_model> <exp_eng_shift> <exp_int_scale>]]")
        sys.exit(1)

    # Parse base arguments
    plot_option = sys.argv[1]
    adjacent_only = sys.argv[2].lower() == 'true'

    # Create output directory
    output_dir = ensure_output_dir('.')

    # Process matrices
    matrices = []
    i = 3
    while i < len(sys.argv):
        if i + 3 >= len(sys.argv):
            break

        filepath = sys.argv[i]
        model = sys.argv[i+1]
        energy_shift = float(sys.argv[i+2])
        interaction_scale = float(sys.argv[i+3])

        matrices.append((filepath, model, energy_shift, interaction_scale))
        i += 4

    # Ensure we have complete sets of arguments
    if len(matrices) not in [1, 2, 3]:
        print("Error: Incomplete matrix parameters")
        sys.exit(1)

    # Process each matrix
    results = {}
    source_names = ['BioDC', 'QM', 'Exp']
    for (filepath, model, energy_shift, interaction_scale), name in zip(matrices, source_names):
        result = process_matrix(name, filepath, model, energy_shift, interaction_scale,
                              adjacent_only, output_dir)
        if result:
            results[name] = result

    # Find global potential range and calculate fractions
    all_potentials = []
    for result in results.values():
        if 'independent' in result:
            all_potentials.append({'independent': result['independent']})
        if 'sequential' in result:
            all_potentials.append({'sequential': result['sequential']})

    lower_bound, upper_bound = find_global_potential_range(all_potentials)
    E_range = np.linspace(lower_bound, upper_bound, 1600)

    # Calculate fractions and prepare for plotting
    all_data = {}
    for name, result in results.items():
        all_data[name] = calculate_fractions(E_range, result, name)

    # Save all data to single file in output directory
    data_dict = {}
    for source_data in all_data.values():
        data_dict.update(source_data)
    save_data_to_text(os.path.join(output_dir, 'redox_data.txt'), E_range, data_dict, [], {})

    # Plot combined curves
    plot_curves(E_range, all_data, os.path.join(output_dir, 'redox_plot.png'), plot_option)
