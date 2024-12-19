import sys
import numpy as np
from datetime import datetime
import derrida  # Import the provided module
from tabulate import tabulate  # For nice table formatting
import time  # For tracking runtime
from typing import List, Tuple, Dict

class PDBProcessor:
    def __init__(self):
        self.heme_atoms = [
            'FE', 'NA', 'C1A', 'C2A', 'C3A', 'C4A', 'CHB', 'C1B', 'NB',
            'C2B', 'C3B', 'C4B', 'CHC', 'C1C', 'NC', 'C2C', 'C3C', 'C4C',
            'CHD', 'C1D', 'ND', 'C2D', 'C3D', 'C4D', 'CHA'
        ]

    def create_pairs_from_sequence(self, seqids: str) -> List[Tuple[int, int]]:
        """Create pairs of consecutive residue IDs from a comma-separated string"""
        residues = [int(resid.strip()) for resid in seqids.split(',')]
        return list(zip(residues[:-1], residues[1:]))

    def read_pdb_atoms(self, pdb_file: str) -> Dict[int, Dict[str, np.ndarray]]:
        """Read atom coordinates from PDB file"""
        atoms_dict = {}
        
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM  ') or line.startswith('HETATM'):
                    try:
                        atom_name = line[12:16].strip()
                        resid = int(line[22:26])
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        
                        if atom_name in self.heme_atoms:
                            if resid not in atoms_dict:
                                atoms_dict[resid] = {}
                            atoms_dict[resid][atom_name] = np.array([x, y, z])
                    except (ValueError, IndexError):
                        continue
        
        return atoms_dict

    def calculate_min_distance(self, coords1: Dict[str, np.ndarray],
                             coords2: Dict[str, np.ndarray]) -> float:
        """Calculate minimum distance between two sets of coordinates"""
        min_dist = float('inf')
        
        for atom1 in self.heme_atoms:
            if atom1 not in coords1:
                continue
            for atom2 in self.heme_atoms:
                if atom2 not in coords2:
                    continue
                    
                dist = np.linalg.norm(coords1[atom1] - coords2[atom2])
                min_dist = min(min_dist, dist)
        
        return min_dist

    def measure_avg_dist(self, pdb_file: str, seqids: str) -> List[float]:
        """Calculate average distances between consecutive heme pairs"""
        residue_pairs = self.create_pairs_from_sequence(seqids)
        atoms_dict = self.read_pdb_atoms(pdb_file)
        
        distances = []
        for res1, res2 in residue_pairs:
            if res1 in atoms_dict and res2 in atoms_dict:
                min_dist = self.calculate_min_distance(atoms_dict[res1], atoms_dict[res2])
                distances.append(min_dist)
            else:
                print(f"Warning: Missing residue {res1} or {res2} in PDB file")
        
        return distances

    def get_average_spacing(self, pdb_file: str, seqids: str) -> float:
        """Calculate the average spacing between hemes"""
        distances = self.measure_avg_dist(pdb_file, seqids)
        if not distances:
            raise ValueError("No valid distances calculated")
        return np.mean(distances)

def parse_arguments(args):
    """Parse command line arguments of the form key=value."""
    params = {}
    for arg in args[1:]:  # Skip script name
        key, value = arg.split('=')
        params[key] = value
    return params

def format_row_data(row):
    """Format each value in the row with consistent spacing."""
    formatted = []
    for i, val in enumerate(row):
        if i == 0:  # Index
            formatted.append(f"{val:6d}")
        elif isinstance(val, str):
            if 'e' in val.lower():  # Scientific notation for rates and D
                num = float(val)
                formatted.append(f"{num:12.2e}")
            else:  # Regular float for couplings, lambdas, deltaG
                num = float(val)
                formatted.append(f"{num:8.3f}")
        else:
            formatted.append(str(val))
    return formatted

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
    PI = 3.141592654
    KB = 8.6173304E-5
    HBAR = 6.582119514E-16
    T = 300.0
    
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

def print_help():
    """Print help information about the script usage."""
    help_text = """
Parameter Exploration Script for Electron Transfer Calculations
------------------------------------------------------------

Usage:
    python script.py [parameters]

Required Parameters:
    num_tot=<int>        : Total number of hemes in sequence
    num_s=<int>          : Number of S-type hemes
    num_t=<int>          : Number of T-type hemes
    seq=<string>         : Sequence of S and T hemes (e.g., 'TSTST')
    num_sets=<int>       : Number of parameter sets to generate
    free_eng_opt=<bool>  : Whether to use free energy optimization ('true' or 'false')
    pdbname=<string>     : Name of PDB file to read
    seqids=<string>      : Comma-separated list of residue IDs for hemes (e.g., '1280,1274,1268,1262,1256')

Optional Parameters:
    d_target=<float>     : Target diffusion coefficient (if not specified, all results are kept)
    d_tolerance=<float>  : Tolerance for matching target D (default: 0.1, meaning ±10%)
    max_attempts=<float> : Maximum number of attempts to find matching sets (default: 1e6)

Parameter Ranges:
    S-type coupling : 3-13 meV
    T-type coupling : 1-4 meV
    Lambda         : 0.1-1.0 eV
    Delta G        : -0.3 to 0.3 eV (only used when free_eng_opt=false)

Examples:
    1. Basic usage without target D:
       python script.py num_tot=5 num_s=2 num_t=3 seq=TSTST num_sets=100000 free_eng_opt=false pdbname=min.pdb seqids=1280,1274,1268,1262,1256

    2. With target D:
       python script.py num_tot=5 num_s=2 num_t=3 seq=TSTST num_sets=1000 free_eng_opt=false pdbname=min.pdb seqids=1280,1274,1268,1262,1256 d_target=1.0e-4 d_tolerance=0.1

Output:
    - Results are written to a file named: <sequence>_<date>.txt
    - Diffusion coefficients are automatically scaled by the square of the average heme spacing
    - Results are sorted by diffusion coefficient (largest to smallest)
    - Progress updates and statistics are printed during execution

Notes:
    - The script will validate that num_s + num_t = num_tot
    - The sequence must match the specified numbers of S and T hemes
    - When using d_target, the acceptance rate may be low if the target is rare in the parameter space
    """
    print(help_text)

def parse_arguments(args):
    """Parse command line arguments of the form key=value."""
    if len(args) == 1 or args[1] in ['--help', '-h', 'help']:
        print_help()
        sys.exit(0)

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

def main():
    # Parse command line arguments
    params = parse_arguments(sys.argv)
    num_tot = int(params['num_tot'])
    num_s = int(params['num_s'])
    num_t = int(params['num_t'])
    sequence = params['seq']
    num_sets = int(float(params['num_sets']))
    free_eng_opt = params['free_eng_opt'].lower() == 'true'
    pdbname = params['pdbname']
    seqids = params['seqids']
    
    # Get average spacing first
    print("Calculating heme spacing from PDB...")
    pdb_processor = PDBProcessor()
    avg_spacing_angstrom = pdb_processor.get_average_spacing(pdbname, seqids)
    # Convert from Å to cm
    avg_spacing_cm = avg_spacing_angstrom * 1e-8
    spacing_factor = avg_spacing_cm * avg_spacing_cm
    print(f"Average heme spacing: {avg_spacing_cm:.2e} cm")
    print(f"Spacing factor (squared): {spacing_factor:.2e} cm²")

    # Parse d_target if provided
    d_target = float(params.get('d_target', 0))  # 0 means no target
    d_tolerance = float(params.get('d_tolerance', 0.1))  # Default 10% tolerance
    max_attempts = int(float(params.get('max_attempts', 1e6)))  # Maximum sampling attempts

    # Validate inputs
    validate_inputs(num_tot, num_s, num_t, sequence)

    # Create output filename with current date
    date_str = datetime.now().strftime("%d-%m-%Y")
    filename = f"{sequence}_{date_str}.txt"

    # Create headers
    headers = create_headers(sequence, free_eng_opt)

    # Write parameter header to both console and file
    param_line = (f"# num_tot={num_tot} num_s={num_s} num_t={num_t} "
                 f"seq={sequence} free_eng_opt={free_eng_opt} "
                 f"pdbname={pdbname} avg_spacing={avg_spacing_cm:.2e}")
    if d_target > 0:
        param_line += f" d_target={d_target:.4e} d_tolerance={d_tolerance:.4e}"
    
    print(param_line)
    with open(filename, 'w') as f:
        f.write(param_line + '\n')

    # Store all data for table
    table_data = []
    
    # Track acceptance statistics
    attempts = 0
    accepted = 0
    start_time = time.time()
    
    # Generate sets until we have enough accepted ones
    while accepted < num_sets and attempts < max_attempts:
        attempts += 1
        
        # Generate parameters
        params = generate_parameters(sequence)
        row_data = [accepted + 1]  # Index
        
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

        # Scale D by spacing factor
        D_scaled = D * spacing_factor
        
        # Check if scaled D is within tolerance of target (if targeting)
        if d_target > 0:
            relative_error = abs(D_scaled - d_target) / d_target
            if relative_error > d_tolerance:
                # Print progress periodically
                if attempts % 10000 == 0:
                    elapsed = time.time() - start_time
                    accept_rate = (accepted / attempts) * 100 if attempts > 0 else 0
                    print(f"Attempts: {attempts}, Accepted: {accepted}, "
                          f"Rate: {accept_rate:.2f}%, Time: {elapsed:.1f}s")
                continue
        
        # If we get here, the set is accepted
        accepted += 1
        
        # Add scaled diffusion constant
        row_data.append(f"{D_scaled:.4e}")
        
        # Store scaled D value for sorting
        table_data.append((D_scaled, row_data))
        
        # Print progress periodically
        if accepted % 100 == 0:
            elapsed = time.time() - start_time
            accept_rate = (accepted / attempts) * 100 if attempts > 0 else 0
            print(f"Found {accepted}/{num_sets} sets, "
                  f"Acceptance rate: {accept_rate:.2f}%, "
                  f"Time: {elapsed:.1f}s")

    # Print final statistics
    elapsed = time.time() - start_time
    accept_rate = (accepted / attempts) * 100 if attempts > 0 else 0
    print(f"\nFinal Statistics:")
    print(f"Total attempts: {attempts}")
    print(f"Accepted sets: {accepted}")
    print(f"Acceptance rate: {accept_rate:.2f}%")
    print(f"Total time: {elapsed:.1f}s")

    if accepted < num_sets:
        print(f"\nWarning: Only found {accepted} sets out of {num_sets} requested")
        print(f"Consider adjusting tolerance or increasing max_attempts")

    # Sort all data by diffusion coefficient in descending order
    print("\nSorting data by diffusion coefficient (largest to smallest)...")
    table_data.sort(key=lambda x: x[0], reverse=True)  # Sort by D value
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
