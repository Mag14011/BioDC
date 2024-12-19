import sympy
import numpy as np
from string import ascii_uppercase
from typing import List, Tuple, Dict, Union, Optional
from functools import lru_cache
import warnings
from scipy.optimize import minimize

warnings.filterwarnings('ignore', category=np.ComplexWarning)

# Global symbols
j = sympy.Symbol('J')

@lru_cache(maxsize=128)
def generate_site_names(n_sites: int) -> List[str]:
    return list(ascii_uppercase[:n_sites]) if n_sites <= 26 else [str(i) for i in range(1, n_sites + 1)]

class HoppingSystem:
    def __init__(self, kfor: List[float], kback: List[float], 
                 input_rate: float=1e15, output_rate: float=1e15):
        self.kfor = np.array(kfor, dtype=np.float64)
        self.kback = np.array(kback, dtype=np.float64)
        self.n_sites = len(kfor) + 1
        self.input_rate = input_rate
        self.output_rate = output_rate
        self.hops = generate_site_names(self.n_sites)
        self.setup_symbols()

    def setup_symbols(self):
        """Setup symbolic variables for the system."""
        # Create population symbols
        self.pop_symbols = [sympy.Symbol(f'p{i}') for i in range(self.n_sites)]
        
        # Create rate symbols and their values
        self.rate_subs = {}
        
        # Input and output rates
        kin = sympy.Symbol('k_in')
        kout = sympy.Symbol('k_out')
        self.rate_subs[kin] = self.input_rate
        self.rate_subs[kout] = self.output_rate
        
        # Forward and backward rates
        self.kf_symbols = []
        self.kb_symbols = []
        for i in range(self.n_sites - 1):
            kf = sympy.Symbol(f'kf{i}')
            kb = sympy.Symbol(f'kb{i}')
            self.kf_symbols.append(kf)
            self.kb_symbols.append(kb)
            self.rate_subs[kf] = self.kfor[i]
            self.rate_subs[kb] = self.kback[i]

    def build_equations(self) -> List[sympy.Eq]:
        """Build the system of equations."""
        equations = []
        
        # First site equation (input)
        equations.append(
            sympy.Eq(j, self.rate_subs[sympy.Symbol('k_in')] * (1 - self.pop_symbols[0]))
        )
        
        # Internal site equations
        for i in range(self.n_sites - 1):
            kf = self.rate_subs[self.kf_symbols[i]]
            kb = self.rate_subs[self.kb_symbols[i]]
            p1 = self.pop_symbols[i]
            p2 = self.pop_symbols[i + 1]
            
            equations.append(
                sympy.Eq(j, kf * p1 * (1 - p2) - kb * p2 * (1 - p1))
            )
        
        # Last site equation (output)
        equations.append(
            sympy.Eq(j, self.rate_subs[sympy.Symbol('k_out')] * self.pop_symbols[-1])
        )
        
        return equations

    def solve_symbolic(self, verbose: bool=True) -> Optional[float]:
        """Get symbolic solution and convert to numerical."""
        if verbose:
            print("Attempting symbolic solution...")
        
        try:
            # Build equations
            equations = self.build_equations()
            
            # Convert to polynomial system
            poly_system = []
            for eq in equations:
                # Move everything to LHS
                lhs = eq.lhs - eq.rhs
                # Multiply by denominators to clear fractions
                poly = sympy.together(lhs).as_numer_denom()[0]
                poly_system.append(poly)
            
            # Solve the system symbolically
            all_vars = self.pop_symbols + [j]
            solutions = sympy.solve(poly_system, all_vars)
            
            if verbose:
                print(f"Found {len(solutions)} symbolic solutions")
            
            # Find valid numerical solutions
            valid_solutions = []
            for sol in solutions:
                try:
                    # Convert to float and check validity
                    nums = [complex(x).real for x in sol]
                    flux = nums[-1]
                    pops = nums[:-1]
                    
                    if flux > 0 and all(0 <= p <= 1 for p in pops):
                        valid_solutions.append((pops, flux))
                except:
                    continue
            
            if valid_solutions:
                # Return the solution with highest flux
                return max(valid_solutions, key=lambda x: x[1])
            
            return None
            
        except Exception as e:
            if verbose:
                print(f"Symbolic solution failed: {str(e)}")
            return None

    def solve_numeric(self, initial_guess: Optional[Tuple[List[float], float]]=None,
                     verbose: bool=True) -> List[float]:
        """Solve numerically with optional initial guess."""
        if verbose:
            print("Attempting numerical solution...")
        
        # Define objective function
        def objective(x):
            flux = x[-1]
            pops = x[:-1]
            
            if flux <= 0 or any(p < 0 or p > 1 for p in pops):
                return 1e20
            
            total_error = 0
            
            # Input flux error
            input_flux = self.input_rate * (1 - pops[0])
            total_error += (input_flux - flux)**2
            
            # Internal fluxes
            for i in range(self.n_sites - 1):
                kf = self.kfor[i]
                kb = self.kback[i]
                p1, p2 = pops[i], pops[i+1]
                
                internal_flux = kf * p1 * (1 - p2) - kb * p2 * (1 - p1)
                total_error += (internal_flux - flux)**2
            
            # Output flux error
            output_flux = self.output_rate * pops[-1]
            total_error += (output_flux - flux)**2
            
            return total_error

        # Set initial guess
        if initial_guess is None:
            x0 = np.zeros(self.n_sites + 1)
            x0[:-1] = 0.5  # Initial populations
            x0[-1] = min(self.kfor.min(), self.kback.min()) / 10  # Initial flux
        else:
            x0 = np.concatenate([initial_guess[0], [initial_guess[1]]])

        # Try multiple optimization methods
        methods = ['Nelder-Mead', 'Powell', 'BFGS']
        best_result = None
        best_error = float('inf')

        for method in methods:
            try:
                result = minimize(
                    objective,
                    x0,
                    method=method,
                    options={'maxiter': 10000}
                )
                
                if result.fun < best_error:
                    best_result = result
                    best_error = result.fun
                    
            except:
                continue

        if best_result is not None and best_error < 1e-6:
            return best_result.x.tolist()
        
        raise ValueError("Could not find valid numerical solution")

    def solve_system(self, verbose: bool=True) -> List[float]:
        """Main solving routine using hybrid approach."""
        try:
            if verbose:
                print(f"Solving system with {self.n_sites} sites...")
            
            # Try symbolic solution first
            symbolic_sol = self.solve_symbolic(verbose)
            
            # Try numeric solution with or without symbolic initial guess
            try:
                numeric_sol = self.solve_numeric(symbolic_sol, verbose)
                if verbose:
                    print(f"Found solution with flux: {numeric_sol[-1]/1e6:.3f} × 10⁶ s⁻¹")
                return numeric_sol
            except:
                if verbose:
                    print("Trying numeric solution without initial guess...")
                numeric_sol = self.solve_numeric(None, verbose)
                if verbose:
                    print(f"Found solution with flux: {numeric_sol[-1]/1e6:.3f} × 10⁶ s⁻¹")
                return numeric_sol
                
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
    # Test with the problematic rates
    kfor = [1.13e8, 5.05e8, 4.85e8, 3.00e8]
    kback = [1.66e7, 1.88e10, 3.80e6, 5.36e9]
    
    try:
        result = solve_flux(kfor, kback)
        print("\nResults:")
        print("Populations:", [f"{p:.6f}" for p in result[:-1]])
        print(f"Flux: {result[-1]/1e6:.3f} × 10⁶ s⁻¹")
    except ValueError as e:
        print(f"Error: {str(e)}")
