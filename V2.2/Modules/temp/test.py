import sys
import numpy as np
from datetime import datetime
import derrida  # Import the provided module
from tabulate import tabulate  # For nice table formatting

# Physical constants
T = 300.0
PI = 3.141592654
KB = 8.6173304E-5
HBAR = 6.582119514E-16

def parse_arguments(args):
    """Parse command line arguments of the form key=value."""
    params = {}
    for arg in args[1:]:  # Skip script name
        key, value = arg.split('=')
        params[key] = value
    return params

def validate_inputs(num_tot, num_s, num_t, sequence):
    """Validate input parameters."""
    if num_s + num_t != num_tot:
        raise ValueError(f"num_s ({num_s}) + num_t ({num_t}) must equal num_tot ({num_tot})")
    
    if len(sequence) != num_tot:
        raise ValueError(f"Sequence length ({len(sequence)}) must equal num_tot ({num_tot})")
    
    s_count = sequence.count('S')
    t_count = sequence.count('T')
    if s_count != num_s or t_count != num_t:
        raise ValueError(f"Sequence S/T counts ({s_count}/{t_count}) don't match num_s/num_t ({num_s}/{num_t})")

def generate_coupling(type_):
    """Generate a random coupling value based on heme type."""
    if type_ == 'S':
        return np.random.uniform(3, 13)  # S-type: 3-13 meV
    else:  # T-type
        return np.random.uniform(1, 4)   # T-type: 1-4 meV

def generate_parameters(sequence):
    """Generate all random parameters needed for the calculation."""
    params = {
        'couplings': [generate_coupling(type_) for type_ in sequence],
        'lambdas': [np.random.uniform(0.1, 1.0) for _ in sequence],  # eV
        'deltaG': [np.random.uniform(-0.3, 0.3) for _ in sequence]   # eV
    }
    return params

def compute_marcus_rate(coupling, lambda_val, deltaG=None):
    """Compute Marcus rate with or without free energy optimization."""
    # coupling is in meV, convert to eV by dividing by 1000
    prefactor = (2 * PI * (coupling/1000)**2) / (HBAR * np.sqrt(4 * PI * lambda_val * KB * T))
    
    if deltaG is None:  # free energy optimized case
        return prefactor
    else:  # full Marcus equation
        return prefactor * np.exp(-(lambda_val + deltaG)**2 / (4 * lambda_val * KB * T))

def create_headers(sequence, free_eng_opt):
    """Create headers based on sequence and calculation type."""
    headers = ["Index"]
    
    # Coupling headers
    headers.extend(f"V{i+1}({t})" for i, t in enumerate(sequence))
    
    # Lambda headers
    headers.extend(f"L{i+1}({t})" for i, t in enumerate(sequence))
    
    if not free_eng_opt:
        # Delta G headers
        headers.extend(f"dG{i+1}({t})" for i, t in enumerate(sequence))
        # Forward and backward rate headers
        headers.extend(f"Rf{i+1}" for i in range(len(sequence)))
        headers.extend(f"Rb{i+1}" for i in range(len(sequence)))
    else:
        # Single rate headers for free energy optimized case
        headers.extend(f"R{i+1}" for i in range(len(sequence)))
    
    headers.append("D")
    return headers

def format_row_data(row):
    """Format each value in the row with consistent spacing."""
    formatted = []
    for i, val in enumerate(row):
        if i == 0:  # Index
            formatted.append(f"{val:6d}")
        elif isinstance(val, str):
            if 'e' in val.lower():  # Scientific notation
                num = float(val)
                formatted.append(f"{num:12.4e}")
            else:  # Regular float
                num = float(val)
                formatted.append(f"{num:8.4f}")
        else:
            formatted.append(str(val))
    return formatted

def main():
    # Parse command line arguments
    params = parse_arguments(sys.argv)
    num_tot = int(params['num_tot'])
    num_s = int(params['num_s'])
    num_t = int(params['num_t'])
    sequence = params['seq']
    num_sets = int(float(params['num_sets']))
    free_eng_opt = params['free_eng_opt'].lower() == 'true'

    # Validate inputs
    validate_inputs(num_tot, num_s, num_t, sequence)

    # Create output filename with current date
    date_str = datetime.now().strftime("%d-%m-%Y")
    filename = f"{sequence}_{date_str}.txt"

    # Create headers
    headers = create_headers(sequence, free_eng_opt)

    # Write parameter header to both console and file
    param_line = f"# num_tot={num_tot} num_s={num_s} num_t={num_t} seq={sequence} free_eng_opt={free_eng_opt}"
    print(param_line)
    
    with open(filename, 'w') as f:
        f.write(param_line + '\n')

    # Store all data for table
    table_data = []

    # Generate all sets and compute results
    for i in range(num_sets):
        # Generate parameters
        params = generate_parameters(sequence)
        row_data = [i + 1]  # Index
        
        # Add couplings
        row_data.extend(f"{v:.4f}" for v in params['couplings'])
        
        # Add lambdas
        row_data.extend(f"{l:.4f}" for l in params['lambdas'])
        
        if free_eng_opt:
            # Compute single set of rates for free energy optimized case
            rates = [compute_marcus_rate(c, l) for c, l in 
                    zip(params['couplings'], params['lambdas'])]
            row_data.extend(f"{r:.4e}" for r in rates)
            
            # Compute diffusion with same forward and backward rates
            V, D = derrida.VD(rates, rates)
            
        else:
            # Add delta G values
            row_data.extend(f"{dg:.4f}" for dg in params['deltaG'])
            
            # Compute forward rates
            forward_rates = [compute_marcus_rate(c, l, dg) for c, l, dg in 
                           zip(params['couplings'], params['lambdas'], params['deltaG'])]
            
            # Compute backward rates (flip sign of deltaG)
            backward_rates = [compute_marcus_rate(c, l, -dg) for c, l, dg in 
                            zip(params['couplings'], params['lambdas'], params['deltaG'])]
            
            row_data.extend(f"{r:.4e}" for r in forward_rates)
            row_data.extend(f"{r:.4e}" for r in backward_rates)
            
            # Compute diffusion with different forward and backward rates
            V, D = derrida.VD(forward_rates, backward_rates)
        
        # If computation failed, skip this set
        if D is None:
            continue
            
        # Add diffusion constant
        row_data.append(f"{D:.4e}")
        
        # Store raw D value for sorting
        table_data.append((D, row_data))
        
        # Print progress
        if (i + 1) % 10000 == 0:
            print(f"Generated {i+1}/{num_sets} sets...")

    # Sort all data by diffusion coefficient in descending order
    print("Sorting data by diffusion coefficient (largest to smallest)...")
    table_data.sort(key=lambda x: x[0], reverse=True)  # Sort by D value, descending
    sorted_data = [row for _, row in table_data]  # Extract just the row data

    # Print and save sorted results in batches
    print("\nWriting sorted results...")
    batch_size = 100
    num_batches = (len(sorted_data) + batch_size - 1) // batch_size

    for batch_idx in range(num_batches):
        start_idx = batch_idx * batch_size
        end_idx = min((batch_idx + 1) * batch_size, len(sorted_data))
        batch = sorted_data[start_idx:end_idx]
        
        # Update indices and format data
        formatted_batch = []
        for i, row in enumerate(batch, start=start_idx+1):
            row[0] = i
            formatted_batch.append(format_row_data(row))
        
        # Create table for this batch with custom formatting
        batch_table = tabulate(
            formatted_batch,
            headers,
            tablefmt="simple",
            numalign="right",
            stralign="right"
        )
        
        # Print and save this batch
        print(batch_table)
        with open(filename, 'a') as f:
            f.write(batch_table + '\n')

if __name__ == "__main__":
    main()
