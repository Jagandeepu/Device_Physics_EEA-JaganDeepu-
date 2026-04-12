# simulate.py
import devsim
from devsim.python_packages.simple_physics import *
import os
import csv
from device_setup import create_1d_pn_mesh, apply_doping

def setup_and_solve_equilibrium(device, region):
    print("--- Setting up Physics (300 K) ---")
    SetSiliconParameters(device, region, 300)
    CreateSolution(device, region, "Potential")
    CreateSiliconPotentialOnly(device, region)
    for contact in devsim.get_contact_list(device=device):
        devsim.set_parameter(device=device, name=GetContactBiasName(contact), value=0.0)
        CreateSiliconPotentialOnlyContact(device, region, contact)
    devsim.solve(type="dc", absolute_error=1.0, relative_error=1e-5, maximum_iterations=50)

def setup_and_solve_drift_diffusion(device, region):
    CreateSolution(device, region, "Electrons")
    CreateSolution(device, region, "Holes")
    devsim.set_node_values(device=device, region=region, name="Electrons", init_from="IntrinsicElectrons")
    devsim.set_node_values(device=device, region=region, name="Holes", init_from="IntrinsicHoles")
    CreateSiliconDriftDiffusion(device, region)
    for contact in devsim.get_contact_list(device=device):
        CreateSiliconDriftDiffusionAtContact(device, region, contact)
    devsim.solve(type="dc", absolute_error=1e10, relative_error=1e-5, maximum_iterations=50)

def save_spatial_data(device, filename):
    """Helper function to save internal profiles at specific voltages."""
    x = devsim.get_node_model_values(device=device, region="silicon", name="x")
    potential = devsim.get_node_model_values(device=device, region="silicon", name="Potential")
    electrons = devsim.get_node_model_values(device=device, region="silicon", name="Electrons")
    holes = devsim.get_node_model_values(device=device, region="silicon", name="Holes")
    
    filepath = os.path.join(os.getcwd(), "results", filename)
    with open(filepath, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["x_cm", "Potential_V", "Electrons_cm3", "Holes_cm3"])
        for i in range(len(x)):
            writer.writerow([x[i], potential[i], electrons[i], holes[i]])

def run_full_sweep(device):
    results_dir = os.path.join(os.getcwd(), "results")
    os.makedirs(results_dir, exist_ok=True)
    
    voltage_list = []
    current_list = []
    
    # 1. Save Equilibrium (0.0V)
    save_spatial_data(device, "spatial_eq.csv")
    
    # 2. Forward Sweep (0.0 to 0.8V)
    print("--- Sweeping Forward Bias ---")
    v = 0.0
    while v <= 0.8 + 1e-9:
        devsim.set_parameter(device=device, name=GetContactBiasName("top"), value=v)
        try:
            devsim.solve(type="dc", absolute_error=1e10, relative_error=1e-4, maximum_iterations=50)
        except devsim.error:
            break
            
        i_elec = devsim.get_contact_current(device=device, contact="top", equation="ElectronContinuityEquation")
        i_hole = devsim.get_contact_current(device=device, contact="top", equation="HoleContinuityEquation")
        voltage_list.append(v)
        current_list.append(i_elec + i_hole)
        
        # Save Forward 0.3V state
        if abs(v - 0.3) < 1e-5:
            save_spatial_data(device, "spatial_fwd03.csv")
        v += 0.05

    # 3. Reset and Reverse Sweep (0.0 to -0.3V)
    print("--- Sweeping Reverse Bias ---")
    devsim.set_parameter(device=device, name=GetContactBiasName("top"), value=0.0)
    devsim.solve(type="dc", absolute_error=1e10, relative_error=1e-5, maximum_iterations=50)
    
    v = -0.05
    while v >= -0.3 - 1e-9:
        devsim.set_parameter(device=device, name=GetContactBiasName("top"), value=v)
        devsim.solve(type="dc", absolute_error=1e10, relative_error=1e-4, maximum_iterations=50)
        
        i_elec = devsim.get_contact_current(device=device, contact="top", equation="ElectronContinuityEquation")
        i_hole = devsim.get_contact_current(device=device, contact="top", equation="HoleContinuityEquation")
        voltage_list.append(v)
        current_list.append(i_elec + i_hole)
        
        # Save Reverse -0.3V state
        if abs(v - (-0.3)) < 1e-5:
            save_spatial_data(device, "spatial_rev03.csv")
        v -= 0.05

    # Save full IV curve
    with open(os.path.join(results_dir, "iv_data.csv"), mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Voltage_V", "Current_Density_A_cm2"])
        # Sort by voltage so the plot draws correctly
        sorted_data = sorted(zip(voltage_list, current_list))
        for volt, curr in sorted_data:
            writer.writerow([volt, curr])
            
    print("Simulation Complete. Data exported for visualization.")

if __name__ == "__main__":
    dev, reg = create_1d_pn_mesh()
    apply_doping(dev, reg)
    setup_and_solve_equilibrium(dev, reg)
    setup_and_solve_drift_diffusion(dev, reg)
    run_full_sweep(dev)