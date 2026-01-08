import devsim
import math

device = "Resistor"
region = "Silicon"

# --- 1. Mesh Construction ---
devsim.create_1d_mesh(mesh="mesh1")
devsim.add_1d_mesh_line(mesh="mesh1", pos=0.0, ps=5e-6, tag="top")
devsim.add_1d_mesh_line(mesh="mesh1", pos=50e-4, ps=5e-6, tag="bot")
devsim.add_1d_contact(mesh="mesh1", name="top", tag="top", material="metal")
devsim.add_1d_contact(mesh="mesh1", name="bot", tag="bot", material="metal")
devsim.add_1d_region(mesh="mesh1", material="Si", region=region, tag1="top", tag2="bot")
devsim.finalize_mesh(mesh="mesh1")
devsim.create_device(mesh="mesh1", device=device)

# --- 2. Parameters (Fonstad Fundamentals) ---
ni, Vt, c1 = 1.0e10, 0.02585, 1.0e17
eps_si, q_val = 1.03594e-12, 1.602e-19
phi_equil = Vt * math.log(c1 / ni)

devsim.set_parameter(name="ni", value=ni)
devsim.set_parameter(name="Vt", value=Vt)
devsim.set_parameter(name="eps_si", value=eps_si)
devsim.set_parameter(name="q", value=q_val)

# --- 3. Solution Variable & Initialization ---
devsim.node_solution(device=device, region=region, name="Potential")
# Initialize with equilibrium potential
devsim.set_node_values(device=device, region=region, name="Potential", values=[phi_equil] * 1001)
devsim.edge_from_node_model(device=device, region=region, node_model="Potential")

# --- 4. Physics Models ---
devsim.node_model(device=device, region=region, name="NetDoping", equation=f"{c1}")
devsim.node_model(device=device, region=region, name="n", equation="ni * exp(Potential / Vt)")
devsim.node_model(device=device, region=region, name="p", equation="ni * exp(-Potential / Vt)")
devsim.node_model(device=device, region=region, name="ChargeDensity", equation="q * (p - n + NetDoping)")

# Jacobian derivative - CRITICAL to stop "Singular Matrix"
devsim.node_model(device=device, region=region, name="ChargeDensity:Potential", 
                  equation="-q * (ni/Vt * exp(Potential/Vt) + ni/Vt * exp(-Potential/Vt))")

devsim.edge_model(device=device, region=region, name="ElectricField", 
                  equation="(Potential@n0 - Potential@n1)*EdgeInverseLength")
devsim.edge_model(device=device, region=region, name="PotentialEdgeFlux", 
                  equation="eps_si * ElectricField")

# --- 5. Potential Equation ---
devsim.equation(device=device, region=region, name="PotentialEquation", 
                variable_name="Potential", node_model="ChargeDensity", 
                edge_model="PotentialEdgeFlux", variable_update="positive")

# --- 6. Boundary Conditions (0.3V Bias) ---
# Top Contact
devsim.contact_node_model(device=device, contact="top", name="top_bc", equation=f"Potential - {phi_equil + 0.3}")
devsim.contact_node_model(device=device, contact="top", name="top_bc:Potential", equation="1")
devsim.contact_equation(device=device, contact="top", name="PotentialEquation", node_model="top_bc")

# Bot Contact
devsim.contact_node_model(device=device, contact="bot", name="bot_bc", equation=f"Potential - {phi_equil}")
devsim.contact_node_model(device=device, contact="bot", name="bot_bc:Potential", equation="1")
devsim.contact_equation(device=device, contact="bot", name="PotentialEquation", node_model="bot_bc")

# --- 7. Solve ---
devsim.solve(type="dc", solver_type="direct", absolute_error=1e-12, relative_error=1e-12, maximum_iterations=50)

# --- 8. Post-Processing ---
devsim.node_model(device=device, region=region, name="delta_n", equation="n - NetDoping")
devsim.write_devices(file="part3_final.vtu", type="vtk")

print("SUCCESS: Result saved to part3_final.vtu")
