import sympy
import numpy as np
from string import ascii_uppercase
from typing import List, Tuple, Dict, Union
from functools import lru_cache

# Define j as a global symbol for the flux
j = sympy.Symbol('J')

# Define input/output rate symbols globally
k_in = sympy.Symbol('k_in')
k_out = sympy.Symbol('k_out')

@lru_cache(maxsize=128)
def generate_site_names(n_sites: int) -> List[str]:
    """
    Generate names for sites (A, B, C, ...). Falls back to numbered names if more than 26 sites.
    Now cached for repeated calls with same n_sites.
    """
    if n_sites <= 26:
        return list(ascii_uppercase[:n_sites])
    return [str(i) for i in range(1, n_sites + 1)]

def build_jji(start: str, stop: str, ps: Dict[str, sympy.Symbol], 
              kforward: Dict[str, sympy.Symbol], kbackward: Dict[str, sympy.Symbol]) -> sympy.Eq:
    """Creates the base equation for a hop from site `start` to `stop`."""
    pj, pi = ps[stop], ps[start]
    kji = kforward[stop+start]
    kij = kbackward[start+stop]
    return sympy.Eq(j, kji*pi*(1-pj) - kij*pj*(1-pi))

def is_effectively_real(value: complex, tolerance: float=1e-10) -> bool:
    """Check if a complex number is effectively real."""
    return abs(getattr(value, 'imag', 0)) < tolerance

def safe_float_conversion(value: Union[complex, float], tolerance: float=1e-10) -> float:
    """Safely convert a possibly complex number to float."""
    if isinstance(value, complex):
        if abs(value.imag) < tolerance:
            return float(value.real)
        raise ValueError(f"Value has significant imaginary component: {value}")
    return float(value)

class HoppingSystem:
    """Class to encapsulate the hopping system state and calculations."""
    
    def __init__(self, kfor: List[float], kback: List[float], 
                 input_rate: float=1e15, output_rate: float=1e15):
        """Initialize the hopping system with given rates."""
        self.kfor = kfor
        self.kback = kback
        self.n_sites = len(kfor) + 1
        self.input_rate = input_rate
        self.output_rate = output_rate
        self.hops = generate_site_names(self.n_sites)
        
        # Pre-calculate commonly used values
        self.setup_rate_dictionaries()
        self.ps = {hop: sympy.Symbol(f'P_{hop}') for hop in self.hops}
        
    def setup_rate_dictionaries(self) -> None:
        """Set up forward and backward rate dictionaries."""
        self.kforward = {}
        self.kbackward = {}
        self.replacements = {
            k_in: self.input_rate,  # Using global symbols instead of creating new ones
            k_out: self.output_rate
        }
        
        for start, stop, (rate_f, rate_b) in zip(
            self.hops[:-1], self.hops[1:], zip(self.kfor, self.kback)):
            kf = sympy.Symbol(f'k_{stop}{start}')
            kb = sympy.Symbol(f'k_{start}{stop}')
            self.kforward[f'{stop}{start}'] = kf
            self.kbackward[f'{start}{stop}'] = kb
            self.replacements[kf] = rate_f
            self.replacements[kb] = rate_b

    def solve_system(self, verbose: bool=True) -> List[float]:
        """Solve the hopping system for populations and flux."""
        try:
            if verbose:
                print(f"  Setting up equations for {self.n_sites} sites...")
            
            # First and last flux relations using global symbols
            initial = sympy.Eq(j, k_in*(1-self.ps[self.hops[0]]))
            final = sympy.Eq(j, k_out*self.ps[self.hops[-1]])
            
            # Solve system sequentially
            population_expressions = []
            substitutes = sympy.solve(initial, self.ps[self.hops[0]])
            if not substitutes:
                raise ValueError("Could not solve initial equation")
            population_expressions.append(substitutes[0])
            
            # Build and solve equations for each hop
            for i, (start, stop) in enumerate(zip(self.hops[:-1], self.hops[1:]), 1):
                if verbose:
                    print(f"    Solving for site {i+1}/{self.n_sites}...")
                base = build_jji(start, stop, self.ps, self.kforward, self.kbackward)
                t = base.subs(self.ps[start], substitutes[0])
                substitutes = sympy.solve(t, self.ps[stop])
                if not substitutes:
                    raise ValueError(f"Could not solve equation for site {i+1}")
                population_expressions.append(substitutes[0])
            
            # Build and solve characteristic polynomial
            t = final.subs(self.ps[self.hops[-1]], substitutes[0])
            jpol = sympy.poly(t.lhs * t.rhs.args[1].args[0] - t.rhs.args[0] * t.rhs.args[2], j)
            base_expression = jpol.args[0].subs(self.replacements)
            coefficients = sympy.poly(base_expression).coeffs()
            
            # Find physically meaningful solution
            roots = np.roots(coefficients)
            real_roots = [root.real for root in roots if is_effectively_real(root) and root.real > 0]
            
            if not real_roots:
                raise ValueError("No positive real roots found")
            
            # Try each real root
            for root in real_roots:
                try:
                    populations = []
                    for population_expression in population_expressions:
                        base_expression = population_expression.subs(self.replacements).subs(j, root)
                        pop_value = safe_float_conversion(complex(base_expression))
                        if not (0 <= pop_value <= 1):
                            break
                        populations.append(pop_value)
                    else:
                        if verbose:
                            print(f"  Found valid solution! Flux: {root/1e6:.3f} × 10⁶ s⁻¹")
                        return populations + [root]
                except (ValueError, TypeError) as e:
                    if verbose:
                        print(f"    Rejected solution: {str(e)}")
                    continue
            
            raise ValueError("No physically meaningful solution found")
            
        except Exception as e:
            raise ValueError(f"Error solving system: {str(e)}")

def solve_flux(kfor: List[float], kback: List[float], verbose: bool=True) -> List[float]:
    """Wrapper function for backward compatibility."""
    try:
        system = HoppingSystem(kfor, kback)
        return system.solve_system(verbose)
    except Exception as e:
        print(f"\nDebug information:")
        print(f"Number of sites: {len(kfor) + 1}")
        print(f"Forward rates (s⁻¹):")
        for i, rate in enumerate(kfor, 1):
            print(f"  k{i}->{i+1}: {rate:.2e}")
        print(f"Backward rates (s⁻¹):")
        for i, rate in enumerate(kback, 1):
            print(f"  k{i+1}->{i}: {rate:.2e}")
        print(f"Error occurred: {str(e)}")
        raise ValueError("Failed to solve the electron transfer system.") from e

if __name__ == "__main__":
    # Example usage with the rates from your error message
    kfor = [1.91e9, 9.80e7, 7.40e8, 3.25e8]
    kback = [1.50e10, 3.78e7, 1.19e10, 2.74e7]
    try:
        result = solve_flux(kfor, kback)
        print("Populations:", result[:-1])
        print("Flux:", result[-1])
    except ValueError as e:
        print(f"Error: {str(e)}")
