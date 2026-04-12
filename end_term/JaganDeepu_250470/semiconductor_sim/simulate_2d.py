# simulate_2d.py
import devsim
from devsim.python_packages.simple_physics import *
import os
import csv

def build_lateral_2d_diode():
    device = "diode_2d"
    region = "silicon"
    devsim.create_2d_mesh(mesh="mesh2d")
    
    # X-axis: 2.0 um wide. Junction is exactly in the middle at x = 1.0
    devsim.add_2d_mesh_line(mesh="mesh2d", dir="x", pos=0.0, ps=0.05)
    devsim.add_2d_mesh_line(mesh="mesh2d", dir="x", pos=1.0, ps=0.01) # Refined at junction
    devsim.add_2d_mesh_line(mesh="mesh2d", dir="x", pos=2.0, ps=0.05)
    
    # Y-axis: 0.5 um deep.
    devsim.add_2d_mesh_line(mesh="mesh2d", dir="y", pos=0.0, ps=0.05)
    devsim.add_2d_mesh_line(mesh="mesh2d", dir="y", pos=0.5, ps=0.05)
    
    devsim.add_2d_region(mesh="mesh2d", material="Si", region=region)
    
    # Contacts are placed on the TOP surface (y=0.5)
    # Anode on the left (x: 0 to 0.4)
    # Anode on the left (x: 0 to 0.4)
    # Anode on the left (x: 0 to 0.4)
    devsim.add_2d_contact(mesh="mesh2d", name="anode", material="metal", 
                          region="silicon", xl=0.0, xh=0.4, yl=0.5, yh=0.5)
                          
    # Cathode on the right (x: 1.6 to 2.0)
    devsim.add_2d_contact(mesh="mesh2d", name="cathode", material="metal", 
                          region="silicon", xl=1.6, xh=2.0, yl=0.5, yh=0.5)
    
    devsim.finalize_mesh(mesh="mesh2d")
    devsim.create_device(mesh="mesh2d", device=device)
    
    # Apply Step Doping: Left side is p-type, Right side is n-type
    devsim.set_parameter(device=device, region=region, name="Na", value=1e16)
    devsim.set_parameter(device=device, region=region, name="Nd", value=1e16)
    devsim.node_model(device=device, region=region, name="NetDoping", equation="ifelse(x < 1.0, -Na, Nd)")
    
    return device, region

def save_2d_state(device, filename):
    x = devsim.get_node_model_values(device=device, region="silicon", name="x")
    y = devsim.get_node_model_values(device=device, region="silicon", name="y")
    pot = devsim.get_node_model_values(device=device, region="silicon", name="Potential")
    elec = devsim.get_node_model_values(device=device, region="silicon", name="Electrons")
    hole = devsim.get_node_model_values(device=device, region="silicon", name="Holes")
    
    os.makedirs("results_2d", exist_ok=True)
    filepath = os.path.join(os.getcwd(), "results_2d", filename)
    with open(filepath, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["x", "y", "Potential", "Electrons", "Holes"])
        for i in range(len(x)):
            writer.writerow([x[i], y[i], pot[i], elec[i], hole[i]])

def run_2d_simulation():
    dev, reg = build_lateral_2d_diode()
    
    print("--- 2D Equilibrium ---")
    SetSiliconParameters(dev, reg, 300)
    CreateSolution(dev, reg, "Potential")
    CreateSiliconPotentialOnly(dev, reg)
    for contact in devsim.get_contact_list(device=dev):
        devsim.set_parameter(device=dev, name=GetContactBiasName(contact), value=0.0)
        CreateSiliconPotentialOnlyContact(dev, reg, contact)
    devsim.solve(type="dc", absolute_error=1.0, relative_error=1e-5, maximum_iterations=50)
    
    CreateSolution(dev, reg, "Electrons")
    CreateSolution(dev, reg, "Holes")
    devsim.set_node_values(device=dev, region=reg, name="Electrons", init_from="IntrinsicElectrons")
    devsim.set_node_values(device=dev, region=reg, name="Holes", init_from="IntrinsicHoles")
    CreateSiliconDriftDiffusion(dev, reg)
    for contact in devsim.get_contact_list(device=dev):
        CreateSiliconDriftDiffusionAtContact(dev, reg, contact)
    devsim.solve(type="dc", absolute_error=1e10, relative_error=1e-5, maximum_iterations=50)
    
    print("--- 2D Forward Bias (+0.6V) ---")
    devsim.set_parameter(device=dev, name=GetContactBiasName("anode"), value=0.6)
    devsim.solve(type="dc", absolute_error=1e10, relative_error=1e-4, maximum_iterations=50)
    save_2d_state(dev, "forward_06.csv")
    
    print("--- 2D Reverse Bias (-0.6V) ---")
    devsim.set_parameter(device=dev, name=GetContactBiasName("anode"), value=-0.6)
    devsim.solve(type="dc", absolute_error=1e10, relative_error=1e-4, maximum_iterations=50)
    save_2d_state(dev, "reverse_06.csv")
    print("2D Simulation Complete!")

if __name__ == "__main__":
    run_2d_simulation()