import os
import argparse
from pathlib import Path

def read_linear_sequence():
    with open('../EE/LinearizedHemeSequence.txt', 'r') as f:
        return [line.split()[1] for line in f]

def read_res_indexing():
    indexing_data = {}
    with open('../EE/ResIndexing.txt', 'r') as f:
        for line in f:
            cols = line.split()
            if len(cols) >= 5:
                indexing_data[cols[4]] = (cols[2], cols[3])
    return indexing_data

def generate_config():
    fe_numbers = read_linear_sequence()
    res_indexing = read_res_indexing()
    
    template = """# Select the topology and trajectory
# Any format recognized by MDAnalysis
# is acceptable. The trajectory should
# be re-imaged.
topology       = BioDC_new.prmtop
trajectory     = min.nc

# Select the type of analysis
# - total         : Compute total electric field vector.
#                   Calculate is parallelized
# - contributions : Compute per-residue and protein/solvent
#                   contributions to the electric field.
#                   Calculation is NOT parallelized.
# - all          : The parallelized total electric field
#                   calculation is performed first, and then
#                   the single-processor contribution analysis
#                   is performed.
analysis       = total

# Select which frame(s) or all to analyze
frames         = all

# Select one or more probe atoms. Add a
# separate probe entry for each probe atom
# NOTE: For the XXX replacement in environment
# to work, probe selections MUST include a
# 'resid' specification (e.g., resid 653)

# Select a distance dependent non-solvent environment (1st command)
# or all non-solvent residues excluding the probe atom (2nd command).
# Either way, XXX is a placeholder that will be replaced by the probe residue
# ID by the program automatically. This feature allows multiple probe atoms
# for separate calculations to be specified in the same input file, while
# ensuring that the probe each time is always excluded from the environment
# definition.
"""
    
    for fe_num in fe_numbers:
        template += f'probe          = (resname HCR and resid {fe_num} and name FE)\n'
        if fe_num in res_indexing:
            phr_num, dhr_num = res_indexing[fe_num]
            template += f'environment    = (not resname WAT and not resname Na+ and not resname Cl-) and not (resname HCR and resid {fe_num}) and not (resname PHR and resid {phr_num}) and not (resname DHR and resid {dhr_num})\n\n'

    template += """
# Select water and solvent ions
solvent        = resname WAT or resname Na+

# Define a cutoff for including solvent molecules
# or a range of cutoffs (START:STOP:STEP Å) that
# will each be assessed.
#solvent_cutoff = 10
cutoff_range  = 1:30:2

# Select the number of processors for the calculation
processors     = 100"""

    with open('input.config', 'w') as f:
        f.write(template)

def get_field_magnitude(residue):
    filename = f'FE-HEH{residue}/ElecField_cutoff_summary.dat'
    try:
        with open(filename, 'r') as f:
            for _ in range(3):
                next(f)
            first_line = next(f)
            return float(first_line.split()[1])
    except (FileNotFoundError, IndexError, ValueError):
        return None

def read_output():
    residues = read_linear_sequence()
    for residue in residues:
        magnitude = get_field_magnitude(residue)
        if magnitude is not None:
            print(f"Residue {residue}: {magnitude:.6f} MV/cm")
        else:
            print(f"Residue {residue}: Error reading magnitude")

def main():
    parser = argparse.ArgumentParser(description='Electric field analysis tool')
    parser.add_argument('mode', choices=['write-input', 'read-output'],
                        help='Operation mode: write-input generates config file, read-output processes results')
    
    args = parser.parse_args()
    
    if args.mode == 'write-input':
        generate_config()
    else:
        read_output()

if __name__ == '__main__':
    main()
