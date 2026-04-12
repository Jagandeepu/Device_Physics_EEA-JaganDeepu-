# visualize.py
import os
import csv
import math
import matplotlib.pyplot as plt

def read_csv(filepath):
    data = {}
    with open(filepath, 'r') as f:
        reader = csv.reader(f)
        headers = next(reader)
        for h in headers:
            data[h] = []
        for row in reader:
            for h, val in zip(headers, row):
                data[h].append(float(val))
    return data

def main():
    results_dir = os.path.join(os.getcwd(), "results")
    
    # Load Data
    iv_data = read_csv(os.path.join(results_dir, "iv_data.csv"))
    sp_eq = read_csv(os.path.join(results_dir, "spatial_eq.csv"))
    sp_fwd = read_csv(os.path.join(results_dir, "spatial_fwd03.csv"))
    sp_rev = read_csv(os.path.join(results_dir, "spatial_rev03.csv"))
    
    x_um = [x * 1e4 for x in sp_eq["x_cm"]]

    # --- Plot 1: Electrostatic Potential (3 Curves) ---
    plt.figure(figsize=(8, 6))
    plt.plot(x_um, sp_eq["Potential_V"], 'k-', label="Equilibrium (0.0 V)")
    plt.plot(x_um, sp_fwd["Potential_V"], 'b--', label="Forward Bias (+0.3 V)")
    plt.plot(x_um, sp_rev["Potential_V"], 'r--', label="Reverse Bias (-0.3 V)")
    plt.title("Plot 1: Electrostatic Potential")
    plt.xlabel("Position x (µm)")
    plt.ylabel("Potential Φ (V)")
    plt.legend()
    plt.grid(True)
    plt.savefig(os.path.join(results_dir, "plot1_potential.png"))
    plt.close()

    # --- Plot 2: Carrier Densities at Equilibrium ---
    plt.figure(figsize=(8, 6))
    plt.semilogy(x_um, sp_eq["Electrons_cm3"], 'b-', label="Electrons (n)", linewidth=2)
    plt.semilogy(x_um, sp_eq["Holes_cm3"], 'r-', label="Holes (p)", linewidth=2)
    # Marking Depletion Edges (Calculated as W/2 ~ 0.215 um from verify.py)
    plt.axvline(x=-0.215, color='k', linestyle=':', label="Depletion Edge (p-side)")
    plt.axvline(x=0.215, color='k', linestyle=':', label="Depletion Edge (n-side)")
    plt.title("Plot 2: Carrier Densities at Equilibrium")
    plt.xlabel("Position x (µm)")
    plt.ylabel("Carrier Density (cm⁻³)")
    plt.legend()
    plt.grid(True, which="both", ls="--", alpha=0.7)
    plt.savefig(os.path.join(results_dir, "plot2_carriers.png"))
    plt.close()

    # --- Plot 3: I-V Curve (Linear & Log Scales) ---
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Linear Subplot (Forward bias only)
    v_fwd = [v for v in iv_data["Voltage_V"] if v >= 0]
    i_fwd = [i for v, i in zip(iv_data["Voltage_V"], iv_data["Current_Density_A_cm2"]) if v >= 0]
    ax1.plot(v_fwd, i_fwd, 'b.-')
    ax1.set_title("Forward Bias (Linear Scale)")
    ax1.set_xlabel("Applied Voltage V")
    ax1.set_ylabel("Current I (A/cm²)")
    ax1.grid(True)

    # Log Subplot (Full Range, absolute current)
    i_abs = [abs(i) if i != 0 else 1e-18 for i in iv_data["Current_Density_A_cm2"]]
    ax2.semilogy(iv_data["Voltage_V"], i_abs, 'r.-')
    ax2.set_title("Full Range (Log Scale)")
    ax2.set_xlabel("Applied Voltage V")
    ax2.set_ylabel("|Current I| (A/cm²)")
    ax2.grid(True, which="both", ls="--")
    
    plt.tight_layout()
    plt.savefig(os.path.join(results_dir, "plot3_iv_curve.png"))
    plt.close()

    # --- Plot 4: Theory vs Simulation (Shockley Overlay) ---
    plt.figure(figsize=(8, 6))
    Vt = 0.02585  # Thermal voltage at 300K
    
    target_v = 0.4
    # Find the index of the voltage closest to 0.4 (handles floating point rounding)
    idx_04 = min(range(len(iv_data["Voltage_V"])), key=lambda i: abs(iv_data["Voltage_V"][i] - target_v))
    actual_v = iv_data["Voltage_V"][idx_04]

    I_sim_04 = iv_data["Current_Density_A_cm2"][idx_04]
    Is_analytical = I_sim_04 / (math.exp(actual_v / Vt) - 1)
    
    # Calculate Shockley Equation points
    shockley_I = [Is_analytical * (math.exp(v / Vt) - 1) for v in v_fwd]
    
    plt.semilogy(v_fwd, i_fwd, 'bo', label="DEVSIM Simulation", alpha=0.6)
    plt.semilogy(v_fwd, shockley_I, 'k--', label=f"Shockley Eq (Is={Is_analytical:.1e})", linewidth=2)
    plt.title("Plot 4: Theory vs Simulation")
    plt.xlabel("Applied Voltage V")
    plt.ylabel("Current I (A/cm²)")
    plt.legend()
    plt.grid(True, which="both", ls="--")
    plt.savefig(os.path.join(results_dir, "plot4_theory.png"))
    plt.close()

    print("Success! The 4 checklist-compliant plots are in the 'results' folder.")

if __name__ == "__main__":
    main()