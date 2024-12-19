import sys
import datetime
import math
import numpy as np
import re
from tabulate import tabulate
import matplotlib.pyplot as plt
import hopping_generalized as hp

def parse_matrix_line(line):
    """Parse a line containing matrix elements in the format [x, y, z, ...]"""
    try:
        clean_line = line.strip('[] \n')
        return [float(x.strip()) for x in clean_line.split(',')]
    except ValueError:
        return None

def read_E_matrix(filename):
    """Read E matrix from file."""
    matrix = []
    with open(filename, 'r') as f:
        for line in f:
            row = parse_matrix_line(line)
            if row is not None:
                matrix.append(row)
    
    # Verify matrix is square
    n = len(matrix)
    if not all(len(row) == n for row in matrix):
        raise ValueError("Input matrix must be square")
    
    return np.array(matrix)

def replicate_array(arr, n_replications):
    """
    Replicate a 1D array n times by repeating all values.

    Args:
        arr (numpy.ndarray): Original 1D array
        n_replications (int): Number of times to replicate

    Returns:
        numpy.ndarray: Replicated array
    """
    result = arr.copy()

    for _ in range(n_replications):
        current_size = len(result)
        # For H and lambda arrays, we repeat the entire array
        result = np.concatenate([result, arr])

    return result

def replicate_matrix(matrix, n_replications):
    """
    Replicate a square matrix n times by sliding the original matrix.

    Args:
        matrix (numpy.ndarray): Original square matrix
        n_replications (int): Number of times to replicate

    Returns:
        numpy.ndarray: Replicated matrix
    """
    orig_size = matrix.shape[0]

    for n in range(n_replications):
        current_size = matrix.shape[0]
        new_size = current_size + (orig_size - 1)
        new_matrix = np.zeros((new_size, new_size))

        # Copy the existing matrix to top-left corner
        new_matrix[:current_size, :current_size] = matrix

        # Fill in the slid matrix (starting from the last position of original)
        for i in range(orig_size):
            for j in range(orig_size):
                if i == 0 and j == 0:
                    # Skip the first element when sliding
                    continue

                new_i = i + (current_size - 1)
                new_j = j + (current_size - 1)

                if new_i < new_size and new_j < new_size:
                    new_matrix[new_i, new_j] = matrix[i, j]

        matrix = new_matrix

    return matrix

def modify_E_matrix(E, diagonal_shift=0.0, offdiagonal_scale=1.0, n_replications=0):
    """
    Modify E matrix according to specified transformations and replications.

    Args:
        E (numpy.ndarray): Input matrix
        diagonal_shift (float): Value to add to diagonal elements (in meV)
        offdiagonal_scale (float): Factor to scale off-diagonal elements
        n_replications (int): Number of times to replicate

    Returns:
        numpy.ndarray: Modified matrix
    """
    E_modified = E.copy()
    n = len(E)

    # Shift diagonal elements
    for i in range(n):
        E_modified[i,i] += diagonal_shift

    # Scale off-diagonal elements
    for i in range(n):
        for j in range(n):
            if i != j:
                E_modified[i,j] *= offdiagonal_scale

    # Perform replication if requested
    if n_replications > 0:
        E_modified = replicate_matrix(E_modified, n_replications)

    return E_modified

def read_H_values(filename):
    """Read H values from file."""
    H_values = []
    pattern = r'Hda\s*=\s*(\d+\.?\d*)\s*meV'
    
    with open(filename, 'r') as f:
        for line in f:
            match = re.search(pattern, line)
            if match:
                H_values.append(float(match.group(1)))
    
    return np.array(H_values)

def load_results(filename):
    """Load results from a file containing JSON data."""
    import json

    with open(filename, 'r') as f:
        content = f.read()
        json_start = content.find('### JSON DATA ###\n') + len('### JSON DATA ###\n')
        json_data = content[json_start:]
        data = json.loads(json_data)

    # Extract the mixed case populations for plotting
    populations = data['flux_results']['mixed']['forward']['populations']

    return (data['flux_results'],
            np.array(data['redox_potentials']['reduced']),
            np.array(data['redox_potentials']['oxidized']),
            np.array(data['redox_potentials']['mixed']),
            populations)

def read_lambda_values(filename):
    """Read reorganization energies from file."""
    lambda_values = []
    pattern = r'Reorg\.\s*Eng\.\s*=\s*(\d+\.?\d*)'
    
    with open(filename, 'r') as f:
        content = f.read()
        matches = re.finditer(pattern, content)
        for match in matches:
            lambda_values.append(float(match.group(1)) * 1000)  # Convert to meV
    
    return np.array(lambda_values)

def format_float(x, precision=3):
    """Format float to string with consistent precision and right alignment."""
    return f"{x:>{precision+5}.{precision}f}"

def print_matrix(matrix, title):
    """Print a matrix in a readable format with consistent number formatting."""
    print(f"\n{title}:")
    headers = [f"Heme {i+1}" for i in range(len(matrix))]
    table = [[f"Heme {i+1}"] + [format_float(x) for x in row] for i, row in enumerate(matrix)]
    print(tabulate(table, headers=[""] + headers, tablefmt="grid", stralign="right"))

def print_rates_table(kfor, kback):
    """Print forward and backward rates in a table with consistent formatting."""
    print("\nElectron Transfer Rates (10⁶ s⁻¹):")
    headers = ["Heme Pair", "Forward Rate", "Backward Rate"]
    table = []
    for i in range(len(kfor)):
        pair = f"{i+1} → {i+2}"
        table.append([
            pair,
            format_float(kfor[i]/1e6),
            format_float(kback[i]/1e6)
        ])
    print(tabulate(table, headers=headers, tablefmt="grid", stralign="right"))

def print_redox_potentials(R_red, R_ox, R_sc):
    """Print redox potentials in a comparative table with consistent formatting."""
    print("\nRedox Potentials (eV):")
    headers = ["Heme", "Fully Reduced", "Fully Oxidized", "Mixed"]
    table = []
    for i in range(len(R_red)):
        table.append([
            f"Heme {i+1}",
            format_float(R_red[i]/1000),
            format_float(R_ox[i]/1000),
            format_float(R_sc[i]/1000)
        ])
    print(tabulate(table, headers=headers, tablefmt="grid", stralign="right"))

def create_plots(R_red, R_ox, R_mixed, populations, output_base, label_step=1,
                 ox_pos=(7.2, -0.2), red_pos=(7.2, -0.45), mixed_pos=(7.2, -0.3)):
    """
    Create and save publication-quality plots with absolute annotation positions.
    X-axis limits are set based on data range only.
    """
    # Set figure size and DPI
    fig = plt.figure(figsize=(3.3, 3.3))

    # Create subplot grid with specific height ratios and spacing
    gs = fig.add_gridspec(2, 1, height_ratios=[1, 1], hspace=0.15)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex=ax1)

    # Get x values (heme indices)
    x = np.arange(1, len(R_red) + 1)

    # Set x limits based only on data
    x_min = x[0] - 0.2
    x_max = x[-1] + 0.2

    # Plot redox potentials (top subplot)
    l1 = ax1.plot(x, R_red/1000, ':', marker='s', mfc='none', color='red')
    l2 = ax1.plot(x, R_ox/1000, ':', marker='s', mfc='none', color='blue')
    l3 = ax1.plot(x, R_mixed/1000, ':', marker='s', mfc='none', color='green')

    # Set x limits for data visibility
    ax1.set_xlim(x_min, x_max)

    # Calculate final points for connection lines
    y_ox = R_ox[-1]/1000
    y_red = R_red[-1]/1000
    y_mixed = R_mixed[-1]/1000

    # Position annotations at absolute coordinates
    ax1.annotate('Ox.', xy=(x[-1], y_ox),
                xytext=ox_pos,
                color='blue', weight='bold', ha='left', va='center')

    ax1.annotate('Red.', xy=(x[-1], y_red),
                xytext=red_pos,
                color='red', weight='bold', ha='left', va='center')

    ax1.annotate('Mixed', xy=(x[-1], y_mixed),  
                xytext=mixed_pos,
                color='green', weight='bold', ha='left', va='center')

    # Configure top subplot
    ax1.set_ylabel(r'$E^{\circ}$ (eV)')
    ax1.tick_params(direction='in', which='both', top=True)
    ax1.spines['bottom'].set_visible(True)  # Restore bottom spine
    plt.setp(ax1.get_xticklabels(), visible=False)  # Hide x-tick labels for top subplot

    # Plot populations (bottom subplot)
    ax2.plot(x, populations, ':k')

    # Plot markers with continuous grayscale
    for i, pop in enumerate(populations):
        gray_value = 1.0 - pop
        ax2.plot(x[i], pop, 'o',
                markerfacecolor=f'{gray_value:.3f}',
                markeredgecolor='black',
                markersize=6)

    # Configure bottom subplot
    ax2.set_xlabel('Heme Index')
    ax2.set_ylabel('Population')
    ax2.tick_params(direction='in', which='both', top=True)
    ax2.set_ylim(-0.05, 1.05)
    ax2.set_xlim(x_min, x_max)  # Same limits as top subplot

    # Configure x-axis labels
    ax2.set_xticks(x[::label_step])
    if label_step > 1:
        current_ticks = list(ax2.get_xticks())
        if 1 not in current_ticks:
            current_ticks = [1] + current_ticks
            ax2.set_xticks(current_ticks)

    # Save the plot
    plt.savefig(f"{output_base}_plots.png", dpi=300, bbox_inches='tight')
    plt.close()

def calculate_directional_flux(redox_potentials, H, lamb, direction='forward'):
    """
    Calculate rates and flux for a given set of redox potentials.

    Parameters:
    redox_potentials: array of redox potentials (in meV)
    H: array of coupling values (in meV)
    lamb: array of reorganization energies (in meV)
    direction: 'forward' (1→N) or 'reverse' (N→1)

    Returns:
    flux, forward_rates, backward_rates, populations
    """
    n_hemes = len(redox_potentials)
    kfor = np.zeros(n_hemes-1)
    kback = np.zeros(n_hemes-1)

    if direction == 'forward':
        # Calculate rates in forward direction
        for i in range(n_hemes-1):
            deltaA = redox_potentials[i] - redox_potentials[i+1]
            kfor[i] = calculate_kji(H[i], lamb[i], deltaA)
            kback[i] = calculate_kji(H[i], lamb[i], -deltaA)
    else:
        # Calculate rates in reverse direction (flip redox potentials)
        redox_potentials = redox_potentials[::-1]
        H = H[::-1]
        lamb = lamb[::-1]
        for i in range(n_hemes-1):
            deltaA = redox_potentials[i] - redox_potentials[i+1]
            kfor[i] = calculate_kji(H[i], lamb[i], deltaA)
            kback[i] = calculate_kji(H[i], lamb[i], -deltaA)

    # Calculate flux
    sol = hp.solve_flux(kfor.tolist(), kback.tolist(), verbose=False)
    return sol[-1], kfor, kback, sol[:-1]  # Return flux, rates, and populations

def write_log_header(log_file, case, direction):
    """Write header for a new calculation to the log file."""
    log_file.write(f"\n{'='*80}\n")
    log_file.write(f"Case: {case}, Direction: {direction}\n")
    log_file.write(f"{'='*80}\n\n")
    log_file.write("Iteration  Max_Change(meV)  Flux(10⁶ s⁻¹)  Redox_Potentials(eV)  Populations\n")
    log_file.write("-" * 80 + "\n")

def write_log_entry(log_file, iteration, max_diff, flux, redox_potentials, populations):
    """Write iteration data to log file."""
    # Format redox potentials and populations as space-separated strings
    potentials_str = " ".join([f"{p/1000:8.3f}" for p in redox_potentials])
    populations_str = " ".join([f"{p:8.3f}" for p in populations])

    log_file.write(f"{iteration:9d}  {max_diff:13.6f}  {flux/1e6:13.6f}  {potentials_str}  {populations_str}\n")

def write_output_file(outfl, flux_results, R_red, R_ox, R_mixed):
    """Write detailed results to output file with consistent formatting."""
    outfl.write("Electron Transfer Analysis Results\n")
    outfl.write("=================================\n\n")

    # Write redox potentials
    outfl.write("Redox Potentials (eV):\n")
    outfl.write("---------------------\n")
    headers = ["Heme", "Fully Reduced", "Fully Oxidized", "Mixed State"]
    table = []
    for i in range(len(R_red)):
        table.append([
            f"Heme {i+1}",
            format_float(R_red[i]/1000),
            format_float(R_ox[i]/1000),
            format_float(R_mixed[i]/1000)
        ])
    outfl.write(tabulate(table, headers=headers, tablefmt="grid", stralign="right"))
    outfl.write("\n\n")

    # Write flux comparison
    outfl.write("Flux Comparison:\n")
    outfl.write("----------------\n")
    headers = ["Case", "Direction", "Flux (10⁶ s⁻¹)", "Rate Constants (10⁶ s⁻¹)"]
    table = []
    for case in ['reduced', 'oxidized', 'mixed']:
        case_name = case.title()
        for direction in ['forward', 'reverse']:
            rates = flux_results[case][direction]['rates']
            rates_str = ", ".join([format_float(r/1e6) for r in rates['forward']] +
                                [format_float(r/1e6) for r in rates['backward']])
            table.append([
                case_name,
                "1→N" if direction == 'forward' else "N→1",
                format_float(flux_results[case][direction]['flux']/1e6),
                rates_str
            ])
    outfl.write(tabulate(table, headers=headers, tablefmt="grid", stralign="right"))
    outfl.write("\n\n")

    # Write mixed state populations
    outfl.write("Mixed State Populations:\n")
    outfl.write("----------------------\n")
    headers = ["Heme", "Population"]
    table = []
    for i, pop in enumerate(flux_results['mixed']['forward']['populations']):
        table.append([f"Heme {i+1}", format_float(pop)])
    outfl.write(tabulate(table, headers=headers, tablefmt="grid", stralign="right"))
    outfl.write("\n")

def save_results(filename, flux_results, R_red, R_ox, R_mixed):
    """Save all results in JSON format after the human-readable output."""
    import json
    
    # Helper function remains the same
    def convert_arrays(obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, dict):
            return {k: convert_arrays(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [convert_arrays(item) for item in obj]
        return obj
    
    flux_results = convert_arrays(flux_results)
    data = {
        'redox_potentials': {
            'reduced': R_red.tolist(),
            'oxidized': R_ox.tolist(),
            'mixed': R_mixed.tolist()
        },
        'flux_results': flux_results
    }
    
    with open(filename, 'w') as f:
        write_output_file(f, flux_results, R_red, R_ox, R_mixed)
        f.write('\n### JSON DATA ###\n')
        json.dump(data, f, indent=2)

def calculate_kji(Hji, lambji, delAji):
    """Calculate electron transfer rate."""
    hbar = 6.582119514e-16  # eV*s
    kB = 8.6173304e-5      # eV/K
    T = 300.0
    kbT = kB * T
    
    # Convert to eV
    Hji = Hji * 0.001
    lambji = lambji * 0.001
    delAji = delAji * 0.001
    
    fourkt = 4 * kbT
    hemeeff = 2 * math.pi/(hbar * math.sqrt(4 * math.pi * kbT))
    exponent = (delAji + lambji)**2/lambji
    exponent = -1.0 * exponent/fourkt
    
    return hemeeff * Hji**2 * math.exp(exponent) / math.sqrt(lambji)

def calculate_max_difference(current, previous):
    """Calculate maximum difference between current and previous values."""
    return np.max(np.abs(current - previous))

def compute_redox_potentials(E, case, populations=None):
    """
    Compute redox potentials based on case and populations.

    Parameters:
    E: interaction matrix
    case: 'reduced', 'oxidized', or 'mixed'
    populations: site populations (only used for 'mixed' case)
    """
    n_hemes = len(E)
    R = np.zeros(n_hemes)

    for i in range(n_hemes):
        if case == 'reduced':
            R[i] = E[i,i]
        elif case == 'oxidized':
            R[i] = E[i,i] + sum(E[i,j] for j in range(n_hemes) if j != i)
        elif case == 'mixed':
            R[i] = E[i,i] + sum(E[i,j] * (1-populations[j]) for j in range(n_hemes) if j != i)

    return R

def compute_redox_potentials_with_mixing(E, case, populations, prev_R, mixing_param=0.3):
    """
    Compute redox potentials with linear mixing to help convergence.

    Parameters:
    E: interaction matrix
    case: 'reduced', 'oxidized', or 'mixed'
    populations: site populations
    prev_R: previous redox potentials
    mixing_param: mixing parameter (0 < α < 1), smaller values give more damping
    """
    n_hemes = len(E)
    R_new = np.zeros(n_hemes)

    # Calculate new potentials
    for i in range(n_hemes):
        if case == 'reduced':
            R_new[i] = E[i,i]
        elif case == 'oxidized':
            R_new[i] = E[i,i] + sum(E[i,j] for j in range(n_hemes) if j != i)
        elif case == 'mixed':
            R_new[i] = E[i,i] + sum(E[i,j] * (1-populations[j]) for j in range(n_hemes) if j != i)

    # Apply linear mixing
    if prev_R is not None:
        R_mixed = mixing_param * R_new + (1 - mixing_param) * prev_R
        return R_mixed
    return R_new

def adaptive_mixing_parameter(iteration, max_diff, prev_max_diff, current_param):
    """
    Adaptively adjust the mixing parameter based on convergence behavior.

    Parameters:
    iteration: current iteration number
    max_diff: current maximum difference
    prev_max_diff: previous maximum difference
    current_param: current mixing parameter
    """
    if iteration < 2:
        return current_param

    # If error is increasing, reduce mixing parameter
    if max_diff > prev_max_diff:
        return max(0.1, current_param * 0.8)
    # If error is decreasing, gradually increase mixing parameter
    else:
        return min(0.8, current_param * 1.1)

def main():
    import argparse
    import os

    parser = argparse.ArgumentParser(description='Calculate and plot electron transfer properties')
    parser.add_argument('E_matrix_file', help='File containing E matrix')
    parser.add_argument('H_file', help='File containing H values')
    parser.add_argument('lambda_file', help='File containing lambda values')
    parser.add_argument('output_base', help='Base name for output files')
    parser.add_argument('--ox-pos', nargs=2, type=float, default=[7.2, -0.2],
                        help='X and Y coordinates for Ox. annotation (e.g., 7.2 -0.2)')
    parser.add_argument('--red-pos', nargs=2, type=float, default=[7.2, -0.45],
                        help='X and Y coordinates for Red. annotation (e.g., 7.2 -0.45)')
    parser.add_argument('--mixed-pos', nargs=2, type=float, default=[7.2, -0.3],
                        help='X and Y coordinates for Mixed annotation (e.g., 7.2 -0.3)')
    parser.add_argument('--label-step', type=int, default=1,
                        help='Step size for x-axis labels (default: 1)')
    parser.add_argument('--shift-diagonal', type=float, default=0.0,
                        help='Uniform shift to apply to diagonal elements of E matrix (meV)')
    parser.add_argument('--scale-offdiagonal', type=float, default=1.0,
                        help='Scale factor to apply to off-diagonal elements of E matrix')
    parser.add_argument('--replicate', type=int, default=0,
                        help='Number of times to replicate the E matrix (default: 0)')
    parser.add_argument('--mixing-param', type=float, default=0.3,
                       help='Initial mixing parameter for convergence (default: 0.3)')
    parser.add_argument('--adaptive-mixing', action='store_true',
                       help='Enable adaptive mixing parameter')
    
    args = parser.parse_args()
    
    # Determine output filenames
    output_txt = f"{args.output_base}.txt"
    log_txt = f"{args.output_base}_iterations.log"

    # Check if output file exists
    if os.path.exists(output_txt):
        print(f"Loading existing results from {output_txt}")
        flux_results, R_red, R_ox, R_mixed, populations = load_results(output_txt)

        # Create plot from existing data
        create_plots(R_red, R_ox, R_mixed, populations, args.output_base, args.label_step,
                    ox_pos=args.ox_pos, red_pos=args.red_pos, mixed_pos=args.mixed_pos)
        print(f"Plot has been saved as {args.output_base}_plots.png")
        return

    # If output doesn't exist, perform calculations
    try:
        # Load the matrices and arrays
        E = read_E_matrix(args.E_matrix_file)
        H = read_H_values(args.H_file)
        lamb = read_lambda_values(args.lambda_file)

        # If replication is requested, replicate all arrays
        if args.replicate > 0:
            # First replicate H and lambda arrays
            H = replicate_array(H, args.replicate)
            lamb = replicate_array(lamb, args.replicate)
            # Then replicate and modify E matrix
            E = modify_E_matrix(E, args.shift_diagonal, args.scale_offdiagonal, args.replicate)
        else:
            # Just modify E matrix without replication
            E = modify_E_matrix(E, args.shift_diagonal, args.scale_offdiagonal)

        n_hemes = len(E)
        Eion = np.zeros(n_hemes)

        print(f"\nSuccessfully loaded parameters for {n_hemes} hemes:")
        print_matrix(E, "Modified E Matrix (meV)")
        print(f"\nH values (meV): {H.tolist()}")
        print(f"Lambda values (meV): {lamb.tolist()}")

    except Exception as e:
        print(f"Error reading input files: {e}")
        sys.exit(1)

    # Dictionary to store results for all cases
    flux_results = {
        'reduced': {'forward': {}, 'reverse': {}},
        'oxidized': {'forward': {}, 'reverse': {}},
        'mixed': {'forward': {}, 'reverse': {}}
    }

    # Store final redox potentials for each case
    final_potentials = {}

    # Open log file
    with open(log_txt, 'w') as log_file:
        log_file.write("Electron Transfer Iterations Log\n")
        log_file.write("==============================\n")
        log_file.write(f"Generated on: {datetime.datetime.now()}\n\n")

        # Process each case
        for case in ['reduced', 'oxidized', 'mixed']:
            print(f"\nCalculating {case} case...")

            # Process both directions
            for direction in ['forward', 'reverse']:
                print(f"  {direction.capitalize()} direction...")
                write_log_header(log_file, case, direction)

                # Set up parameters for this direction
                if direction == 'forward':
                    current_H = H.copy()
                    current_lamb = lamb.copy()
                else:
                    current_H = H[::-1]
                    current_lamb = lamb[::-1]

                # Initialize with zero populations
                populations = np.zeros(n_hemes)

                # Compute initial redox potentials
                R = compute_redox_potentials(E, case, populations)
                if direction == 'reverse':
                    R = R[::-1]

                # Enter self-consistent loop
                iteration = 0
                max_iterations = 1000
                threshold = 1

                while True:
                    prev_R = R.copy()

                    # Calculate rates
                    kfor = np.zeros(n_hemes-1)
                    kback = np.zeros(n_hemes-1)
                    for i in range(n_hemes-1):
                        deltaA = R[i] - R[i+1]
                        kfor[i] = calculate_kji(current_H[i], current_lamb[i], deltaA)
                        kback[i] = calculate_kji(current_H[i], current_lamb[i], -deltaA)

                    # Get flux and populations
                    sol = hp.solve_flux(kfor.tolist(), kback.tolist())
                    populations = sol[:-1]
                    flux = sol[-1]

                    # Update redox potentials according to case
                    R = compute_redox_potentials(E, case, populations)
                    if direction == 'reverse':
                        R = R[::-1]

                    # Check convergence
                    max_diff = calculate_max_difference(R, prev_R)

                    # Write log entry
                    write_log_entry(log_file, iteration, max_diff, flux, R, populations)

                    if max_diff < threshold:
                        print(f"    Converged after {iteration+1} iterations!")
                        print(f"    Final maximum change: {max_diff:.6f} meV")
                        break

                    iteration += 1
                    if iteration >= max_iterations:
                        print("\n    Warning: Maximum iterations reached without convergence")
                        print(f"    Current maximum change: {max_diff:.6f} meV")
                        break

                # Store results
                flux_results[case][direction] = {
                    'flux': flux,
                    'rates': {'forward': kfor, 'backward': kback},
                    'populations': populations
                }

                # Store final potentials (only need one direction)
                if direction == 'forward':
                    final_potentials[case] = R.copy()

                # Add a blank line in log after each calculation
                log_file.write("\n")

    # Save results and create plots
    save_results(output_txt, flux_results,
                final_potentials['reduced'],
                final_potentials['oxidized'],
                final_potentials['mixed'])

    # Create plots using mixed case forward populations
    populations = flux_results['mixed']['forward']['populations']
    create_plots(final_potentials['reduced'],
                final_potentials['oxidized'],
                final_potentials['mixed'],
                populations,
                args.output_base,
                args.label_step,
                ox_pos=args.ox_pos,
                red_pos=args.red_pos,
                mixed_pos=args.mixed_pos)

    print(f"\nResults have been written to {output_txt}")
    print(f"Iteration log has been written to {log_txt}")
    print(f"Plot has been saved as {args.output_base}_plots.png")

if __name__ == "__main__":
    main()
