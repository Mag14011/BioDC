import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import sys

def read_dg_file(filename):
    dg_values = []
    steps = []
    with open(filename, 'r') as f:
        for line in f:
            parts = line.split(';')
            dg = float(parts[1].split('=')[1].strip().split()[0])
            step = parts[0].strip()
            dg_values.append(dg)
            steps.append(step)
    return np.array(dg_values), steps

def read_lambda_file(filename):
    lambda_values = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith(" Reorg. Eng."):
                lambda_values.append(float(line.split('=')[1].strip()))
    return np.array(lambda_values)

def read_hda_file(filename):
    hda_values = []
    with open(filename, 'r') as f:
        for line in f:
            hda = float(line.split('=')[2].strip().split()[0])
            hda_values.append(hda)
    return np.array(hda_values)

def read_rates_file(filename):
    rates = []
    with open(filename, 'r') as f:
        next(f)  # Skip header
        for line in f:
            ketf, ketb = map(float, line.strip().split(','))
            rates.append(ketf)
    return np.array(rates)

def calculate_activation_energy(dg, lamda):
    return (lamda + dg)**2 / (4 * lamda)

def extract_step_label(step_string):
    start, end = step_string.split('->')
    start_num = start.split('-')[1].strip()[:-1]  # Remove the last character (closing parenthesis)
    end_num = end.split('-')[1].strip()[:-1]  # Remove the last character (closing parenthesis)
    return f"{start_num}->{end_num}"

def create_plot(dg_file, lambda_file, hda_file, rates_file, use_resid=True, output_file=None):
    # Read data from files
    dg_values, steps = read_dg_file(dg_file)
    lambda_values = read_lambda_file(lambda_file)
    hda_values = read_hda_file(hda_file)
    rate_values = read_rates_file(rates_file)

    # Calculate activation energy
    activation_energy = calculate_activation_energy(dg_values, lambda_values)

    # Table:
    # Define column headers and widths
    headers = ["Step", "ΔG (eV)", "λ (eV)", "E_a (eV)", "H_DA (meV)", "Rate (s^-1)"]
    widths = [10, 9, 8, 10, 11, 10]

    # Print header
    header_line = " | ".join(f"{h:{w}s}" for h, w in zip(headers, widths))
    print(header_line)
    print("-" * len(header_line))

    # Print data rows
    for i, step in enumerate(steps):
        if use_resid:
            step_label = extract_step_label(step)
        else:
            step_label = str(i + 1)
        
        row = [
            f"{step_label:{widths[0]}s}",
            f"{dg_values[i]:>{widths[1]}.3f}",
            f"{lambda_values[i]:>{widths[2]}.3f}",
            f"{activation_energy[i]:>{widths[3]}.3f}",
            f"{hda_values[i]:>{widths[4]}.3f}",
            f"{rate_values[i]:>{widths[5]}.2e}"
        ]
        print(" | ".join(row))

    # Create plot
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(3.3, 3.3), sharex=True)
    
    # Remove space between subplots
    fig.subplots_adjust(hspace=0)
    
    x = range(len(steps))

    # Colors
    rate_color = 'darkred'
    activation_color = 'darkgreen'
    coupling_color = 'darkblue'

    # Marker size
    marker_size = 4

    # Top plot: Electron transfer rate
    ax1.semilogy(x, rate_values, 'o:', color=rate_color, mfc='lightcoral', mec=rate_color, markersize=marker_size)
    ax1.set_ylabel('$K_{ET}$ (s$^{-1}$)', color=rate_color)
    ax1.tick_params(axis='y', colors=rate_color, direction='in', which='both')
#   ax1.yaxis.set_minor_locator(AutoMinorLocator())
    ax1.yaxis.set_minor_locator(plt.LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1))
    for spine in ax1.spines.values():
        spine.set_edgecolor(rate_color)
    
    # Remove bottom ticks of top subplot
    ax1.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)

    # Bottom plot: Activation energy and Electronic coupling
    ax2.plot(x, activation_energy, 'o:', color=activation_color, mfc='lightgreen', mec=activation_color, markersize=marker_size)
    ax2.set_ylabel('$E_a$ (eV)', color=activation_color)
    ax2.tick_params(axis='y', colors=activation_color, direction='in', which='both')
    ax2.yaxis.set_minor_locator(AutoMinorLocator())
    for spine in ax2.spines.values():
        spine.set_edgecolor(activation_color)

    ax3 = ax2.twinx()
    ax3.plot(x, hda_values, 'o:', color=coupling_color, mfc='lightblue', mec=coupling_color, markersize=marker_size)
    ax3.set_ylabel('$H_{DA}$ (meV)', rotation=270, labelpad=15, color=coupling_color)
    ax3.tick_params(axis='y', colors=coupling_color, direction='in', which='both')
    ax3.yaxis.set_minor_locator(AutoMinorLocator())
    ax3.spines['right'].set_edgecolor(coupling_color)

    # Set x-axis labels
    ax2.set_xlabel('Electron Transfer Step')
    ax2.set_xticks(x)
    
    if use_resid:
        ax2.set_xticklabels([extract_step_label(step) for step in steps], rotation=90, ha='center')
    else:
        ax2.set_xticklabels([str(i+1) for i in x], rotation=0, ha='center')
    
    ax2.tick_params(axis='x', direction='in')

    # Adjust layout to prevent clipping of labels
    plt.tight_layout()
    
    # Fine-tune the plot to remove any remaining space
    plt.subplots_adjust(hspace=0)
    
    if output_file is None:
        output_file = f'electron_transfer_plot_{"resid" if use_resid else "stepnum"}.png'
    
    plt.savefig(output_file, dpi=600, bbox_inches='tight')
    plt.close(fig)  # Close the figure to free up memory

def main():
    if len(sys.argv) < 6:
        print("Usage: python electron_transfer_plot.py /path/DG.txt /path/Lambda.txt /path/Hda.txt /path/rates.txt use_resid=True/False [/path/output.png]")
        sys.exit(1)
    
    dg_file = sys.argv[1]
    lambda_file = sys.argv[2]
    hda_file = sys.argv[3]
    rates_file = sys.argv[4]
    
    use_resid = sys.argv[5].split("=")[1].lower() == "true" if sys.argv[5].startswith("use_resid=") else True
    
    output_file = sys.argv[6] if len(sys.argv) > 6 else None
    
    create_plot(dg_file, lambda_file, hda_file, rates_file, use_resid, output_file)

if __name__ == "__main__":
    main()
