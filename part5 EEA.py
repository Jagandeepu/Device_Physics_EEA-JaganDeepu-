import devsim
import math

device = "Diode"
region = "Silicon"

# --- 1. Mesh Construction ---
devsim.create_1d_mesh(mesh="mesh1")
devsim.add_1d_mesh_line(mesh="mesh1", pos=0.0, ps=1e-5, tag="left")
devsim.add_1d_mesh_line(mesh="mesh1", pos=50e-4, ps=1e-5, tag="right")
devsim.add_1d_contact(mesh="mesh1", name="anode", tag="left", material="metal")
devsim.add_1d_contact(mesh="mesh1", name="cathode", tag="right", material="metal")
devsim.add_1d_region(mesh="mesh1", material="Si", region=region, tag1="left", tag2="right")
devsim.finalize_mesh(mesh="mesh1")
devsim.create_device(mesh="mesh1", device=device)

# --- 2. Parameters ---
ni, Vt, conc = 1.0e10, 0.02585, 1.0e17
eps_si, q_val = 1.03594e-12, 1.602e-19
phi_p = -Vt * math.log(conc / ni) 
phi_n = Vt * math.log(conc / ni)  

devsim.set_parameter(name="ni", value=ni)
devsim.set_parameter(name="Vt", value=Vt)
devsim.set_parameter(name="eps_si", value=eps_si)
devsim.set_parameter(name="q", value=q_val)

# --- 3. Step Function Doping ---
devsim.node_model(device=device, region=region, name="NetDoping", 
                  equation=f"ifelse(x < 25e-4, -{conc}, {conc})")

# --- 4. Initialization ---
devsim.node_solution(device=device, region=region, name="Potential")
x_coords = devsim.get_node_model_values(device=device, region=region, name="x")
init_pot = [phi_p if x < 25e-4 else phi_n for x in x_coords]
devsim.set_node_values(device=device, region=region, name="Potential", values=init_pot)
devsim.edge_from_node_model(device=device, region=region, node_model="Potential")

# --- 5. Physics Models ---
devsim.node_model(device=device, region=region, name="n", equation="ni * exp(Potential / Vt)")
devsim.node_model(device=device, region=region, name="p", equation="ni * exp(-Potential / Vt)")
devsim.node_model(device=device, region=region, name="ChargeDensity", equation="q * (p - n + NetDoping)")
devsim.node_model(device=device, region=region, name="ChargeDensity:Potential", 
                  equation="-q * (ni/Vt * exp(Potential/Vt) + ni/Vt * exp(-Potential/Vt))")

devsim.edge_model(device=device, region=region, name="ElectricField", 
                  equation="(Potential@n0 - Potential@n1)*EdgeInverseLength")
devsim.edge_model(device=device, region=region, name="PotentialEdgeFlux", 
                  equation="eps_si * ElectricField")

# --- 6. Equation ---
devsim.equation(device=device, region=region, name="PotentialEquation", 
                variable_name="Potential", node_model="ChargeDensity", 
                edge_model="PotentialEdgeFlux", variable_update="log_damp")

# --- 7. Initial Boundary Conditions (Equilibrium) ---
for contact, pot in [("anode", phi_p), ("cathode", phi_n)]:
    devsim.contact_node_model(device=device, contact=contact, name=f"{contact}_bc", equation=f"Potential - {pot}")
    devsim.contact_node_model(device=device, contact=contact, name=f"{contact}_bc:Potential", equation="1")
    devsim.contact_equation(device=device, contact=contact, name="PotentialEquation", node_model=f"{contact}_bc")

# --- 8. Bias Sweep Loop ---
# Voltages: 0, Forward (0.1, 0.3, 0.5), Reverse (-0.1, -0.3, -0.5)
bias_list = [0.0, 0.1, 0.3, 0.5, -0.1, -0.3, -0.5]

print("\n--- Starting PN Junction Bias Sweep ---")
print("Bias (V) | Status")
print("------------------")

for v in bias_list:
    # Update Anode Potential: Equilibrium Potential + Bias
    devsim.contact_node_model(device=device, contact="anode", name="anode_bc", 
                              equation=f"Potential - ({phi_p + v})")
    
    try:
        devsim.solve(type="dc", solver_type="direct", absolute_error=1e-12, 
                     relative_error=1e-12, maximum_iterations=50)
        
        # Calculate Excess Carrier Concentration: delta_n = n(v) - n(equilibrium)
        # Note: We define n_intrinsic_p as n at equilibrium on p-side for delta_n analysis
        devsim.node_model(device=device, region=region, name="delta_n", 
                          equation=f"n - (ni * exp(({phi_p if v >= 0 else phi_p})/Vt))")
        
        # Save VTU file for each bias
        fname = f"diode_{v}V.vtu".replace("-", "neg")
        devsim.write_devices(file=fname, type="vtk")
        print(f"{v:7.1f}  | Success -> {fname}")
        
    except Exception as e:
        print(f"{v:7.1f}  | Failed: {e}")

print("------------------")
print("All files saved. Process complete.")
