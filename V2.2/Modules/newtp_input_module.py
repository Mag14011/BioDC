import math
import re
import os

def find_file(directory, pattern):
    for filename in os.listdir(directory):
        if re.match(pattern, filename, re.IGNORECASE):
            return os.path.join(directory, filename)
    raise FileNotFoundError(f"No file matching pattern '{pattern}' found in directory '{directory}'")

def read_hda_values(filename):
    hda_values = []
    with open(filename, 'r') as file:
        for line in file:
            match = re.search(r'Hda =\s+(\d+\.\d+)\s+meV', line)
            if match:
                hda_meV = float(match.group(1))
                hda_J = hda_meV * 1.60218e-22  # Convert meV to J
                hda_values.append(hda_J)
    return hda_values

def read_dg_values(filename):
    dg_values = []
    with open(filename, 'r') as file:
        for line in file:
            match = re.search(r'DG =\s+([-]?\d+\.\d+)\s+eV', line)
            if match:
                dg_eV = float(match.group(1))
                dg_values.append(dg_eV)  # Keep in eV
    return dg_values

def read_lambda_values(filename):
    lambda_values = []
    with open(filename, 'r') as file:
        for line in file:
            if 'Reorg. Eng.' in line:
                lambda_eV = float(line.split('=')[1].strip())
                lambda_values.append(lambda_eV)  # Keep in eV
    return lambda_values

def read_r_values(filename):
    r_values = []
    with open(filename, 'r') as file:
        for line in file:
            match = re.search(r'r\s*=\s*([-]?\d+\.?\d*)\s*Å', line)
            if match:
                r_angstrom = float(match.group(1))
                r_values.append(r_angstrom)
    return r_values

def calculate_prefactor_calccoup(Hda, lambda_rxn, T):
    h_bar = 1.054571817e-34  # Reduced Planck's constant in J·s
    k_b = 1.380649e-23  # Boltzmann constant in J/K
    lambda_rxn_J = lambda_rxn * 1.60218e-19  # Convert eV to J
    return (2 * math.pi * (Hda)**2) / (h_bar * math.sqrt(4 * math.pi * lambda_rxn_J * k_b * T))

def calculate_prefactor_moser_dutton(r):
    A0 = 1e13  # s^-1
    r0 = 3.6  # Å
    y = 0.6  # Å^-1
    return A0 * math.exp(-y * (r - r0))

def calculate_k_et_components(Hda, delta_G, lambda_rxn, T, r, calc_method):
    k_b = 1.380649e-23  # Boltzmann constant in J/K
    delta_G_J = delta_G * 1.60218e-19  # Convert eV to J
    lambda_rxn_J = lambda_rxn * 1.60218e-19  # Convert eV to J

    if calc_method == 'CalcCoup':
        prefactor = calculate_prefactor_calccoup(Hda, lambda_rxn, T)
    elif calc_method == 'MoserDuttonRuler':
        prefactor = calculate_prefactor_moser_dutton(r)
    else:
        raise ValueError("Invalid calculation method. Use 'CalcCoup' or 'MoserDuttonRuler'.")

    activation_energy = (delta_G_J + lambda_rxn_J)**2 / (4 * lambda_rxn_J)
    exponent = -activation_energy / (k_b * T)
    exponential_term = math.exp(exponent)
    k_et = prefactor * exponential_term

    return prefactor, activation_energy, exponential_term, k_et

def write_output_file(filename, prefactors, dg_values, lambda_values):
    N = len(prefactors)
    with open(filename, 'w') as file:
        file.write(f"{N}\n")
        for g, dg, lam in zip(prefactors, dg_values, lambda_values):
            file.write(f"{g*1E-12:<15.3f}{dg:<15.3f}{lam:<15.3f} ! g, DeltaG, Lambda\n")
        file.write("2.00            2.00                          ! g in,out\n")
        file.write("0.0             0.0                           ! DeltaG in,out\n")
        file.write("0.4             0.4                           ! Lambda in,out\n")
        file.write(f"H1 to H{N}\n")

def print_results_table(results, calc_method):
    print(f"Calculation method: {calc_method}")
    print(f"{'Index':^5} | {'Hda (eV)':^12} | {'Prefactor (s⁻¹)':^20} | {'ΔG (eV)':^10} | {'λ (eV)':^10} | {'E_a (eV)':^10} | {'Exponential':^15} | {'k_ET (s⁻¹)':^15}")
    print('-' * 117)
    for result in results:
        index, hda, prefactor, dg, lambda_rxn, activation_energy, exponential_term, k_et = result
        hda_eV = hda / 1.60218e-19  # Convert J to eV
        print(f"{index:^5} | {hda_eV:^12.6f} | {prefactor:^20.2e} | {dg:^10.3f} | {lambda_rxn:^10.3f} | {activation_energy/1.60218e-19:^10.3f} | {exponential_term:^15.2e} | {k_et:^15.2e}")

def process_data(directory, calc_method, T=298):
    hda_file = find_file(f"{directory}/EE", r'hda\.txt')
    dg_file = find_file(f"{directory}/EE", r'dg\.txt')
    lambda_file = find_file(f"{directory}/EE", r'lambda\.txt')
    r_file = find_file(f"{directory}/EE", r'r\.txt')
    output_file = os.path.join(f"{directory}/RCP", f'output_{calc_method}.txt')

    hda_values = read_hda_values(hda_file)
    dg_values = read_dg_values(dg_file)
    lambda_values = read_lambda_values(lambda_file)
    r_values = read_r_values(r_file)

    results = []
    for i, (hda, dg, lambda_rxn, r) in enumerate(zip(hda_values, dg_values, lambda_values, r_values), start=1):
        prefactor, activation_energy, exponential_term, k_et = calculate_k_et_components(hda, dg, lambda_rxn, T, r, calc_method)
        results.append((i, hda, prefactor, dg, lambda_rxn, activation_energy, exponential_term, k_et))

    write_output_file(output_file, [r[2] for r in results], dg_values, lambda_values)
    print_results_table(results, calc_method)
    print(f"\nOutput file '{output_file}' has been created successfully.")

    return results
