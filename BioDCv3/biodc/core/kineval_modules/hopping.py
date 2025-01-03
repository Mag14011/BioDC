import sympy
import numpy as np
import math
from string import ascii_uppercase

# Define j as a global symbol for the flux
j = sympy.symbols('J')

def generate_site_names(n_sites):
    """Generate names for sites (A, B, C, ...). Falls back to numbered names if more than 26 sites."""
    if n_sites <= 26:
        return ascii_uppercase[:n_sites]
    return [str(i) for i in range(1, n_sites + 1)]

def build_jji(start, stop, ps, kforward, kbackward):
    """Creates the base equation for a hop from site `start` to `stop`."""
    pj, pi = ps[stop], ps[start]
    kji = kforward[stop+start]
    kij = kbackward[start+stop]
    return sympy.Eq(j, kji*pi*(1-pj)-kij*pj*(1-pi))

def is_effectively_real(value, tolerance=1e-10):
    """Check if a complex number is effectively real."""
    if isinstance(value, complex):
        return abs(value.imag) < tolerance
    return True

def safe_float_conversion(value, tolerance=1e-10):
    """Safely convert a possibly complex number to float."""
    if isinstance(value, complex):
        if abs(value.imag) < tolerance:
            return float(value.real)
        raise ValueError(f"Value has significant imaginary component: {value}")
    return float(value)

def is_valid_population(value, tolerance=1e-6):
    """
    Check if a population value is physically valid within tolerance.
    
    Parameters:
    value: float - population value to check
    tolerance: float - acceptable deviation from [0,1] bounds
    
    Returns:
    bool - True if population is valid within tolerance
    """
    return -tolerance <= value <= 1 + tolerance

def solve_flux(kfor, kback, verbose=False, population_tolerance=1e-6):
    """
    Solve for maximum flux in a linear chain of n sites.
    
    Parameters:
    kfor: list of forward rates [k21, k32, k43, ...]
    kback: list of backward rates [k12, k23, k34, ...]
    verbose: whether to print progress information
    population_tolerance: tolerance for population bounds (default 1e-6)
    
    Returns:
    list of populations for each site plus the flux [P1, P2, ..., Pn, J]
    """
    # Number of sites is one more than number of rates
    n_sites = len(kfor) + 1
    
    if verbose:
        print(f"  Setting up equations for {n_sites} sites...")
    
    # Generate site names
    hops = generate_site_names(n_sites)
    
    # Set input/output rates (large values for maximum flux calculation)
    input_rate = 1e15
    output_rate = 1e15
    
    # Create rate pairs
    rates = list(zip(kfor, kback))
    
    # Initialize symbol dictionaries
    kforward = {}
    kbackward = {}
    ps = {}
    replacements = {}
    
    # Define symbols
    kin, kout = sympy.symbols(r'k_\text{in} k_\text{out}')
    replacements[kin] = input_rate
    replacements[kout] = output_rate
    
    # Create rate symbols and their replacements
    for start, stop, rate in zip(hops[:-1], hops[1:], rates):
        kf = sympy.symbols(f'k_{stop}{start}')
        kb = sympy.symbols(f'k_{start}{stop}')
        key_f = f'{stop}{start}'
        key_b = f'{start}{stop}'
        kforward[key_f] = kf
        kbackward[key_b] = kb
        replacements[kf] = rate[0]
        replacements[kb] = rate[1]
    
    # Create population symbols
    for hop in hops:
        ps[hop] = sympy.symbols(f'P_{hop}')
    
    try:
        if verbose:
            print("  Building system of equations...")
            
        # First and last flux relations
        initial = sympy.Eq(j, kin*(1-ps[hops[0]]))
        final = sympy.Eq(j, kout*ps[hops[-1]])
        
        # Solve system of equations
        population_expressions = []
        substitutes = sympy.solve(initial, ps[hops[0]])
        population_expressions.append(substitutes[0])
        
        if verbose:
            print("  Solving equations sequentially...")
            
        # Build and solve equations for each hop
        for i, (start, stop) in enumerate(zip(hops[:-1], hops[1:]), 1):
            if verbose:
                print(f"    Solving for site {i+1}/{n_sites}...")
            base = build_jji(start, stop, ps, kforward, kbackward)
            t = base.subs(ps[start], substitutes[0])
            substitutes = sympy.solve(t, ps[stop])
            population_expressions.append(substitutes[0])
        
        if verbose:
            print("  Building characteristic polynomial...")
            
        # Substitute final population
        t = final.subs(ps[hops[-1]], substitutes[0])
        
        # Build and solve characteristic polynomial
        jpol = sympy.poly(t.lhs * t.rhs.args[1].args[0] - t.rhs.args[0] * t.rhs.args[2], j)
        base_expression = jpol.args[0].subs(replacements)
        coefficients = sympy.poly(base_expression).coeffs()
        
        if verbose:
            print("  Finding roots of characteristic polynomial...")
            
        roots = np.roots(coefficients)
        
        # Find physically meaningful solution
        if verbose:
            print(f"  Found {len(roots)} roots, checking for physical solutions...")
            
        real_roots = []
        for root in roots:
            if is_effectively_real(root):
                root_real = root.real
                if root_real > 0:
                    real_roots.append(root_real)
        
        if not real_roots:
            raise ValueError("No positive real roots found. The system may be physically impossible.")
        
        if verbose:
            print(f"  Found {len(real_roots)} positive real roots, checking population constraints...")
        
        # Try each real root
        for i, root in enumerate(real_roots, 1):
            if verbose:
                print(f"  Testing solution {i}/{len(real_roots)}...")
            try:
                populations = []
                valid_solution = True
                for site, population_expression in enumerate(population_expressions, 1):
                    base_expression = population_expression.subs(replacements).subs(j, root)
                    pop_value = safe_float_conversion(complex(base_expression))
                    if not is_valid_population(pop_value, population_tolerance):
                        if verbose:
                            print(f"    Rejected: Site {site} population = {pop_value:.6f} (out of bounds)")
                        valid_solution = False
                        break
                    populations.append(pop_value)
                
                if valid_solution:
                    # If we get here, we found a valid solution
                    if verbose:
                        print("  Found valid solution!")
                        print(f"    Flux: {root/1e6:.3f} × 10⁶ s⁻¹")
                    return populations + [root]
                
            except (ValueError, TypeError) as e:
                if verbose:
                    print(f"    Rejected: {str(e)}")
                continue  # Try next root if this one fails
        
        raise ValueError("No physically meaningful solution found. All solutions lead to invalid populations.")
        
    except Exception as e:
        print(f"\nDebug information:")
        print(f"Number of hemes: {n_sites}")
        print(f"Forward rates (s⁻¹):")
        for i, rate in enumerate(kfor, 1):
            print(f"  k{i}->{i+1}: {rate:.2e}")
        print(f"Backward rates (s⁻¹):")
        for i, rate in enumerate(kback, 1):
            print(f"  k{i+1}->{i}: {rate:.2e}")
        print(f"Error occurred: {str(e)}")
        raise ValueError("Failed to solve the electron transfer system. Check the rates and system parameters.")

if __name__ == "__main__":
    # Example usage
    kfor = [1e6, 1e6, 1e6]  # Example forward rates
    kback = [1e5, 1e5, 1e5]  # Example backward rates
    result = solve_flux(kfor, kback)
    print("Populations:", result[:-1])
    print("Flux:", result[-1])
