import sys
import numpy as np
from datetime import datetime, timedelta
import derrida  # Import the provided module
from tabulate import tabulate  # For nice table formatting
import time  # For tracking runtime
from typing import List, Tuple, Dict
from multiprocessing import Pool, Manager, Value, Lock
import ctypes
import os

# Physical constants
T = 300.0
PI = 3.141592654
KB = 8.6173304E-5
HBAR = 6.582119514E-16

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

class SharedCounter:
    """A counter that can be shared between processes"""
    def __init__(self):
        self.val = Value('i', 0)
        self.lock = Lock()

    def increment(self):
        with self.lock:
            self.val.value += 1

    def value(self):
        with self.lock:
            return self.val.value

def parse_arguments(args):
    """Parse command line arguments of the form key=value."""
    if len(args) == 1 or args[1] in ['--help', '-h', 'help']:
        print_help()
        sys.exit(0)

    params = {}
    for arg in args[1:]:
        key, value = arg.split('=')
        # Handle boolean values
        if value.lower() in ['true', 'false']:
            params[key] = value
        # Handle scientific notation and float values
        elif 'E' in value.upper() or 'e' in value:
            params[key] = float(value)
        elif value.replace('.','').isdigit():
            params[key] = float(value)
        else:
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

def worker_process(shared_data):
    """Worker function for parallel processing"""
    print(f"Worker started with process ID: {os.getpid()}")

    try:
        (sequence, free_eng_opt, spacing_factor, d_target, d_tolerance,
         shared_results, attempts_counter, accepted_counter, num_sets, max_attempts) = shared_data

        print(f"Worker {os.getpid()} parameters:")
        print(f"  free_eng_opt: {free_eng_opt}")
        print(f"  spacing_factor: {spacing_factor}")
        print(f"  d_target: {d_target}")
        print(f"  max_attempts: {max_attempts}")

        local_attempts = 0

        while local_attempts < max_attempts:
            local_attempts += 1
            attempts_counter.increment()

            if accepted_counter.value() >= num_sets:
                print(f"Worker {os.getpid()} stopping: target number of sets reached")
                break

            # Generate parameters
            params = generate_parameters(sequence)

            # Calculate D value
            try:
                if free_eng_opt:
                    rates = [compute_marcus_rate(c, l) for c, l in
                            zip(params['couplings'], params['lambdas'])]
                    V, D = derrida.VD(rates, rates)
                else:
                    forward_rates = [compute_marcus_rate(c, l, dg) for c, l, dg in
                                   zip(params['couplings'], params['lambdas'], params['deltaG'])]
                    backward_rates = [compute_marcus_rate(c, l, -dg) for c, l, dg in
                                    zip(params['couplings'], params['lambdas'], params['deltaG'])]
                    V, D = derrida.VD(forward_rates, backward_rates)

                if D is not None:
                    D_scaled = D * spacing_factor

                    # Accept all results if no target, or check against target
                    accept = True
                    if d_target > 0:
                        relative_error = abs(D_scaled - d_target) / d_target
                        accept = relative_error <= d_tolerance

                    if accept:
                        # Format row data
                        row_data = [0]  # Index will be updated later
                        row_data.extend(f"{v:.3f}" for v in params['couplings'])
                        row_data.extend(f"{l:.3f}" for l in params['lambdas'])

                        if not free_eng_opt:
                            row_data.extend(f"{dg:.3f}" for dg in params['deltaG'])
                            row_data.extend(f"{r:.2e}" for r in forward_rates)
                            row_data.extend(f"{r:.2e}" for r in backward_rates)
                        else:
                            row_data.extend(f"{r:.2e}" for r in rates)

                        row_data.append(f"{D_scaled:.2e}")

                        accepted_counter.increment()
                        shared_results.append((D_scaled, row_data))

                        if accepted_counter.value() % 1000 == 0:
                            print(f"Worker {os.getpid()}: Found {accepted_counter.value()} sets")

            except Exception as e:
                print(f"Error in worker {os.getpid()}: {e}")

        print(f"Worker {os.getpid()} finished: attempts={local_attempts}, "
              f"accepted={accepted_counter.value()}")

    except Exception as e:
        print(f"Fatal error in worker {os.getpid()}: {e}")
        raise

def run_parallel_sampling(num_processes, num_sets, sequence, free_eng_opt,
                        spacing_factor, d_target, d_tolerance, max_attempts):
    """Run sampling in parallel with progress tracking"""
    print_flush("Starting parallel sampling...")
    
    with Manager() as manager:
        shared_results = manager.list()
        attempts_counter = SharedCounter()
        accepted_counter = SharedCounter()

        max_attempts_per_worker = int(max_attempts // num_processes)
        print_flush(f"Max attempts per worker: {max_attempts_per_worker}")
        
        shared_data = (sequence, free_eng_opt, spacing_factor, d_target, d_tolerance,
                      shared_results, attempts_counter, accepted_counter,
                      num_sets, max_attempts_per_worker)

        print_flush("Creating worker pool...")
        with Pool(num_processes) as pool:
            try:
                result = pool.apply_async(worker_process, (shared_data,))  # Changed from map_async
                
                while not result.ready():
                    attempts = attempts_counter.value()
                    accepted = accepted_counter.value()
                    accept_rate = (accepted / attempts * 100) if attempts > 0 else 0.0
                    print_flush(f"Progress: Attempts={attempts}, Accepted={accepted}, "
                              f"Rate={accept_rate:.2f}%")
                    time.sleep(1.0)

                result.get()  # Will raise any exceptions that occurred
                print_flush("All workers finished")
                
            except Exception as e:
                print_flush(f"Error in parallel processing: {str(e)}")
                raise

        results = list(shared_results)
        results.sort(key=lambda x: x[0], reverse=True)
        return results, attempts_counter.value()

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
    num_sets=<float>     : Number of parameter sets to generate (can use scientific notation, e.g., 1E6)
    free_eng_opt=<bool>  : Whether to use free energy optimization ('true' or 'false')
    pdbname=<string>     : Name of PDB file to read
    seqids=<string>      : Comma-separated list of residue IDs for hemes (e.g., '1280,1274,1268,1262,1256')

Optional Parameters:
    d_target=<float>     : Target diffusion coefficient (can use scientific notation)
    d_tolerance=<float>  : Tolerance for matching target D (default: 0.1, meaning ±10%)
    max_attempts=<float> : Maximum number of attempts to find matching sets (default: 1e6)
    num_processes=<int>  : Number of CPU processes to use (default: all available)

Parameter Ranges:
    S-type coupling : 3-13 meV
    T-type coupling : 1-4 meV
    Lambda         : 0.1-1.0 eV
    Delta G        : -0.3 to 0.3 eV (only used when free_eng_opt=false)

Examples:
    1. Basic usage without target D:
       python script.py num_tot=5 num_s=2 num_t=3 seq=TSTST num_sets=1E6 free_eng_opt=false \\
           pdbname=min.pdb seqids=1280,1274,1268,1262,1256

    2. With target D and parallel processing:
       python script.py num_tot=5 num_s=2 num_t=3 seq=TSTST num_sets=1E6 free_eng_opt=false \\
           pdbname=min.pdb seqids=1280,1274,1268,1262,1256 d_target=2.05E-4 \\
           d_tolerance=0.1 max_attempts=1E6 num_processes=4

Notes:
    - The script automatically scales diffusion coefficients by the square of the average heme spacing
    - Results are sorted by diffusion coefficient (largest to smallest)
    - Progress updates show acceptance rate and estimated time remaining
    - When using d_target, the acceptance rate may be low if the target is rare in the parameter space
    """
    print(help_text)

def main():
    # Parse command line arguments
    params = parse_arguments(sys.argv)
    num_tot = int(params['num_tot'])
    num_s = int(params['num_s'])
    num_t = int(params['num_t'])
    sequence = params['seq']
    num_sets = int(params['num_sets'])
    free_eng_opt = str(params['free_eng_opt']).lower() == 'true'
    pdbname = params['pdbname']
    seqids = params['seqids']
    max_attempts = int(params.get('max_attempts', 1e6))
    num_processes = int(params.get('num_processes', os.cpu_count()))
    
    print(f"Running with {num_processes} processes...")

    # Validate inputs
    validate_inputs(num_tot, num_s, num_t, sequence)

    # Get average spacing first
    print("Calculating heme spacing from PDB...")
    pdb_processor = PDBProcessor()
    avg_spacing_angstrom = pdb_processor.get_average_spacing(pdbname, seqids)
    # Convert from Å to cm
    avg_spacing_cm = avg_spacing_angstrom * 1e-8
    spacing_factor = avg_spacing_cm * avg_spacing_cm
    print(f"Average heme spacing: {avg_spacing_angstrom:.2f} Å")
    print(f"Spacing factor (squared): {spacing_factor:.2e} cm²")

    # Create output filename with current date
    date_str = datetime.now().strftime("%d-%m-%Y")
    filename = f"{sequence}_{date_str}.txt"

    # Create headers
    headers = create_headers(sequence, free_eng_opt)

    # Parse d_target if provided
    d_target = float(params.get('d_target', 0))
    d_tolerance = float(params.get('d_tolerance', 0.1))

    # Write parameter header to both console and file
    param_line = (f"# num_tot={num_tot} num_s={num_s} num_t={num_t} "
                 f"seq={sequence} free_eng_opt={free_eng_opt} "
                 f"pdbname={pdbname} avg_spacing={avg_spacing_angstrom:.2f}")
    if d_target > 0:
        param_line += f" d_target={d_target:.2e} d_tolerance={d_tolerance:.2e}"
    
    print(param_line)
    with open(filename, 'w') as f:
        f.write(param_line + '\n')

    start_time = time.time()

    try:
        # Run parallel sampling
        results, total_attempts = run_parallel_sampling(
            num_processes=num_processes,
            num_sets=num_sets,
            sequence=sequence,
            free_eng_opt=free_eng_opt,
            spacing_factor=spacing_factor,
            d_target=d_target,
            d_tolerance=d_tolerance,
            max_attempts=max_attempts
        )

        # Print final statistics
        elapsed = time.time() - start_time
        if total_attempts > 0:
            accept_rate = (len(results) / total_attempts * 100)
        else:
            accept_rate = 0.0

        print(f"\nFinal Statistics:")
        print(f"Total attempts: {total_attempts:,d}")
        print(f"Accepted sets: {len(results):,d}")
        print(f"Acceptance rate: {accept_rate:.2f}%")
        print(f"Total time: {elapsed:.1f}s")
        
        if len(results) > 0:
            print(f"Average time per accepted set: {elapsed/len(results):.3f}s")
            print(f"Processing rate: {total_attempts/elapsed:.1f} attempts/second")

            # Print headers once at the start
            header_table = tabulate([], headers, tablefmt="simple",
                                  numalign="right", stralign="right")
            print("\n" + header_table)
            with open(filename, 'a') as f:
                f.write('\n' + header_table + '\n')

            # Sort results by diffusion coefficient
            sorted_data = [row for _, row in results]

            # Print and save sorted results in batches
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
                
                # Create table for this batch WITHOUT headers
                batch_table = tabulate(
                    formatted_batch,
                    tablefmt="simple",
                    numalign="right",
                    stralign="right"
                )
                
                # Print and save this batch
                print(batch_table)
                with open(filename, 'a') as f:
                    f.write(batch_table)
                    if batch_idx < num_batches - 1:
                        f.write('\n')

            print(f"\nResults have been saved to: {filename}")
        else:
            print("\nNo valid results found.")
            if d_target > 0:
                print(f"Consider adjusting d_target ({d_target:.2e}) or tolerance ({d_tolerance:.2f})")
            else:
                print("Check your input parameters and constraints.")

    except Exception as e:
        print(f"\nError during execution: {str(e)}")
        raise

if __name__ == "__main__":
    main()

