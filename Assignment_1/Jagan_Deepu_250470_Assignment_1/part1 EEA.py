
import devsim

device = "Resistor"
region = "Silicon"

# --- 1. Mesh Creation ---
# 1D mesh (0 to 50um) with 100 node points
# Spacing (ps) = 50e-4 / 100 = 5e-7
devsim.create_1d_mesh(mesh="mesh1")
devsim.add_1d_mesh_line(mesh="mesh1", pos=0.0, ps=5e-6, tag="top")
devsim.add_1d_mesh_line(mesh="mesh1", pos=50e-4, ps=5e-6, tag="bot")
devsim.add_1d_contact(mesh="mesh1", name="top", tag="top", material="metal")
devsim.add_1d_contact(mesh="mesh1", name="bot", tag="bot", material="metal")
devsim.add_1d_region(mesh="mesh1", material="Si", region=region, tag1="top", tag2="bot")
devsim.finalize_mesh(mesh="mesh1")
devsim.create_device(mesh="mesh1", device=device)

# --- 2. Parameters ---
ni = 1.0e10        
Vt = 0.02585       
c1 = 1.0e17        # Magnitude of Donor concentration
c2 = 0.0           # Magnitude of Acceptor concentration (per prompt example)

devsim.set_parameter(device=device, region=region, name="ni", value=ni)
devsim.set_parameter(device=device, region=region, name="Vt", value=Vt)

# --- 3. Net Doping Profile ---
# Following your prompt: Net doping = (c1 * Ndn) - (c2 * Nap)
# For a uniform resistor (Fonstad 6.2), Ndn = 1 across the whole device
devsim.node_model(device=device, region=region, name="NDn", equation=f"{c1}")
devsim.node_model(device=device, region=region, name="NAp", equation=f"{c2}")
devsim.node_model(device=device, region=region, name="NetDoping", equation="NDn - NAp")

# --- 4. Potential Solution (Zero Bias) ---
devsim.node_solution(device=device, region=region, name="Potential")
devsim.node_model(device=device, region=region, name="n", equation="ni * exp(Potential / Vt)")
devsim.node_model(device=device, region=region, name="p", equation="ni * exp(-Potential / Vt)")

for c in ("top", "bot"):
    devsim.contact_node_model(device=device, contact=c, name=f"{c}_bc", 
                              equation="Potential - Vt * log(NetDoping / ni)")
    devsim.contact_node_model(device=device, contact=c, name=f"{c}_bc:Potential", equation="1")
    devsim.contact_equation(device=device, contact=c, name="PotentialEquation", node_model=f"{c}_bc")

# --- 5. Initial Guess and Solve ---
devsim.node_model(device=device, region=region, name="Equil", equation="Vt * log(NetDoping / ni)")
vals = devsim.get_node_model_values(device=device, region=region, name="Equil")
devsim.set_node_values(device=device, region=region, name="Potential", values=vals)

# --- 6. Excess Carrier Models ---
devsim.node_model(device=device, region=region, name="n0", equation="ni * exp(Vt * log(NetDoping / ni) / Vt)")
devsim.node_model(device=device, region=region, name="p0", equation="ni * exp(-Vt * log(NetDoping / ni) / Vt)")
devsim.node_model(device=device, region=region, name="delta_n", equation="n - n0")
devsim.node_model(device=device, region=region, name="delta_p", equation="p - p0")

devsim.write_devices(file="part1_resistor_final.vtu", type="vtk")
print("SUCCESS: Resistor profile generated according to Fonstad 6.2.")
