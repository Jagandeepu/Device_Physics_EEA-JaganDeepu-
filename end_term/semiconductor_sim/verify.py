# verify.py
# Analytical calculations for a 1D Silicon p-n junction

import math

# --- Physical Constants for Silicon at 300 K ---
q = 1.602e-19        # Elementary charge (Coulombs)
kB = 1.3806e-23      # Boltzmann constant (J/K)
T = 300.0            # Temperature (Kelvin)
ni = 1.0e10          # Intrinsic carrier concentration (cm^-3)
eps_0 = 8.854e-14    # Vacuum permittivity (F/cm)
eps_r = 11.7         # Relative permittivity of Silicon
eps_s = eps_r * eps_0

def calc_built_in_potential(Na, Nd):
    """Calculates the built-in potential (V_bi) of a p-n junction."""
    Vt = (kB * T) / q  # Thermal voltage
    Vbi = Vt * math.log((Na * Nd) / (ni**2))
    return Vbi

def calc_depletion_width(Na, Nd, Vbi):
    """Calculates the total depletion width (W) in centimeters."""
    # W = sqrt( (2 * eps_s / q) * (1/Na + 1/Nd) * Vbi )
    term1 = (2 * eps_s) / q
    term2 = (1.0 / Na) + (1.0 / Nd)
    W = math.sqrt(term1 * term2 * Vbi)
    return W

def main():
    # Define your chosen doping profiles
    Na = 1e16  # Acceptor concentration on p-side (cm^-3)
    Nd = 1e16  # Donor concentration on n-side (cm^-3)
    
    print("--- Analytical Theory Targets: 1D p-n Junction ---")
    print(f"p-side Doping (Na): {Na:.1e} cm^-3")
    print(f"n-side Doping (Nd): {Nd:.1e} cm^-3")
    print("-" * 50)
    
    # Calculate values
    Vbi = calc_built_in_potential(Na, Nd)
    W_cm = calc_depletion_width(Na, Nd, Vbi)
    
    # Convert width to micrometers for easier reading/meshing
    W_um = W_cm * 1e4
    
    print(f"Built-in Potential (V_bi) : {Vbi:.4f} Volts")
    print(f"Total Depletion Width (W) : {W_um:.4f} micrometers")
    print("-" * 50)
    print("Use these values to verify your DEVSIM simulation output.")

if __name__ == "__main__":
    main()
