# visualize_2d.py
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

def compute_and_plot_current(filename, title, output_name):
    filepath = os.path.join(os.getcwd(), "results_2d", filename)
    if not os.path.exists(filepath):
        print(f"Missing {filename}")
        return

    df = pd.read_csv(filepath)
    
    # 1. Create a uniform grid for plotting
    xi = np.linspace(df['x'].min(), df['x'].max(), 100)
    yi = np.linspace(df['y'].min(), df['y'].max(), 50)
    X, Y = np.meshgrid(xi, yi)

    # Interpolate DEVSIM node data onto our uniform grid
    V = griddata((df['x'], df['y']), df['Potential'], (X, Y), method='linear')
    N = griddata((df['x'], df['y']), df['Electrons'], (X, Y), method='linear')
    P = griddata((df['x'], df['y']), df['Holes'], (X, Y), method='linear')

    # 2. Calculate Spatial Gradients (dx, dy)
    dx = xi[1] - xi[0]
    dy = yi[1] - yi[0]
    
    dV_dy, dV_dx = np.gradient(V, dy, dx)
    dN_dy, dN_dx = np.gradient(N, dy, dx)
    dP_dy, dP_dx = np.gradient(P, dy, dx)

    # Electric Field: E = -∇V
    Ex, Ey = -dV_dx, -dV_dy

    # Physical Constants
    q = 1.602e-19
    Vt = 0.02585
    mu_n = 1000 # Approx e- mobility
    mu_p = 400  # Approx h+ mobility

    # 3. Calculate Current Densities (Drift + Diffusion)
    # J_n = q * mu_n * (n*E + Vt * ∇n)
    Jnx = q * mu_n * (N * Ex + Vt * dN_dx)
    Jny = q * mu_n * (N * Ey + Vt * dN_dy)
    
    # J_p = q * mu_p * (p*E - Vt * ∇p)
    Jpx = q * mu_p * (P * Ex - Vt * dP_dx)
    Jpy = q * mu_p * (P * Ey - Vt * dP_dy)

    # Total Current
    Jx = Jnx + Jpx
    Jy = Jny + Jpy

    # Magnitude for coloring the stream lines
    magnitude = np.sqrt(Jx**2 + Jy**2)
    log_mag = np.log10(np.clip(magnitude, 1e-15, None)) # Clip to avoid log(0)

    # 4. Plotting
    plt.figure(figsize=(10, 4))
    
    # Draw the contact locations (Anode left, Cathode right)
    plt.plot([0.0, 0.4], [0.5, 0.5], 'r-', linewidth=6, label="Anode (p-type)")
    plt.plot([1.6, 2.0], [0.5, 0.5], 'b-', linewidth=6, label="Cathode (n-type)")
    plt.axvline(x=1.0, color='k', linestyle='--', alpha=0.5, label="Metallurgical Junction")

    # Streamplot (Vector Field)
    strm = plt.streamplot(X, Y, Jx, Jy, color=log_mag, cmap='plasma', 
                          linewidth=1.5, density=1.2, arrowsize=1.5)
    
    plt.colorbar(strm.lines, label="Log10( Current Density Magnitude )")
    
    plt.title(f"2D Current Vector Field: {title}")
    plt.xlabel("Position X (µm)")
    plt.ylabel("Position Y (µm) [Depth]")
    plt.legend(loc="lower left")
    plt.xlim(0, 2)
    plt.ylim(0, 0.52)
    
    out_path = os.path.join(os.getcwd(), "results_2d", output_name)
    plt.savefig(out_path, bbox_inches='tight', dpi=200)
    plt.close()
    print(f"Saved: {output_name}")

if __name__ == "__main__":
    os.makedirs("results_2d", exist_ok=True)
    compute_and_plot_current("forward_06.csv", "Forward Bias (+0.6V)", "bonus_2d_forward.png")
    compute_and_plot_current("reverse_06.csv", "Reverse Bias (-0.6V)", "bonus_2d_reverse.png")