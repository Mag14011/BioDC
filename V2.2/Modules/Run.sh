#!/bin/bash

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
       python script.py num_tot=5 num_s=2 num_t=3 seq=TSTST num_sets=1E6 free_eng_opt=false \
           pdbname=min.pdb seqids=1280,1274,1268,1262,1256

    2. With target D and parallel processing:
       python script.py num_tot=5 num_s=2 num_t=3 seq=TSTST num_sets=1E6 free_eng_opt=false \
           pdbname=min.pdb seqids=1280,1274,1268,1262,1256 d_target=2.05E-4 \
           d_tolerance=0.1 max_attempts=1E6 num_processes=4

Notes:
    - The script automatically scales diffusion coefficients by the square of the average heme spacing
    - Results are sorted by diffusion coefficient (largest to smallest)
    - Progress updates show acceptance rate and estimated time remaining
    - When using d_target, the acceptance rate may be low if the target is rare in the par


python ParameterExploration.py \
	num_tot=5 \
	num_s=2 \
	num_t=3 \
	seq=TSTST \
	num_sets=1000 \
	free_eng_opt=false \
	d_target="2.053061e-04" \
	d_tolerance=0.1 \
	max_attempts=1000000

#python ParameterExploration.py num_tot=5 num_s=2 num_t=3 seq=TSTST num_sets=1 free_eng_opt=true 

