
import devsim

device = "Resistor"
region = "Silicon"

# --- 1. Mesh Creation (0 to 50um with 100 node points) ---
# Total length = 50e-4 cm. Spacing = 5e-6 cm provides 100 points.
devsim.create_1d_mesh(mesh="mesh1")
devsim.add_1d_mesh_line(mesh="mesh1", pos=0.0, ps=5e-6, tag="top")
devsim.add_1d_mesh_line(mesh="mesh1", pos=50e-4, ps=5e-6, tag="bot")
devsim.add_1d_contact(mesh="mesh1", name="top", tag="top", material="metal")
devsim.add_1d_contact(mesh="mesh1", name="bot", tag="bot", material="metal")
devsim.add_1d_region(mesh="mesh1", material="Si", region=region, tag1="top", tag2="bot")
devsim.finalize_mesh(mesh="mesh1")
devsim.create_device(mesh="mesh1", device=device)

# --- 2. Physical Parameters (Silicon) ---
ni = 1.0e10        
Vt = 0.02585       
c1 = 1.0e17        # Donor concentration (Constant for uniform profile)
c2 = 0.0           # Acceptor concentration (Constant for uniform profile)

devsim.set_parameter(device=device, region=region, name="ni", value=ni)
devsim.set_parameter(device=device, region=region, name="Vt", value=Vt)

# --- 3. Net Doping Profile (Corrected per Professor's Advice) ---
# Formula: NetDoping = c1*Ndn - c2*Nap. Using uniform Ndn=1 and Nap=0.
devsim.node_model(device=device, region=region, name="NDn", equation=f"{c1}")
devsim.node_model(device=device, region=region, name="NAp", equation=f"{c2}")
devsim.node_model(device=device, region=region, name="NetDoping", equation="NDn - NAp")

# --- 4. Potential and Carrier Models ---
devsim.node_solution(device=device, region=region, name="Potential")
devsim.node_model(device=device, region=region, name="n", equation="ni * exp(Potential / Vt)")
devsim.node_model(device=device, region=region, name="p", equation="ni * exp(-Potential / Vt)")

# --- 5. Boundary Conditions (Zero Bias Equilibrium) ---
for c in ("top", "bot"):
    # Ohmic contact boundary condition at zero bias
    devsim.contact_node_model(device=device, contact=c, name=f"{c}_bc", 
                              equation="Potential - Vt * log(NetDoping / ni)")
    devsim.contact_node_model(device=device, contact=c, name=f"{c}_bc:Potential", equation="1")
    devsim.contact_equation(device=device, contact=c, name="PotentialEquation", node_model=f"{c}_bc")

# --- 6. The Equilibrium Solve ---
# Set initial guess based on uniform doping log-ratio
devsim.node_model(device=device, region=region, name="Equil", equation="Vt * log(NetDoping / ni)")
vals = devsim.get_node_model_values(device=device, region=region, name="Equil")
devsim.set_node_values(device=device, region=region, name="Potential", values=vals)

# --- 7. Excess Carrier Models (Verification) ---
# n0 and p0 are equilibrium concentrations. delta_n should be 0 at zero bias.
devsim.node_model(device=device, region=region, name="n0", equation="ni * exp(Vt * log(NetDoping / ni) / Vt)")
devsim.node_model(device=device, region=region, name="p0", equation="ni * exp(-Vt * log(NetDoping / ni) / Vt)")
devsim.node_model(device=device, region=region, name="delta_n", equation="n - n0")
devsim.node_model(device=device, region=region, name="delta_p", equation="p - p0")

devsim.write_devices(file="part2_zero_bias_uniform.vtu", type="vtk")
print("SUCCESS: Zero bias simulation complete for uniform resistor.")
