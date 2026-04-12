# device_setup.py
# Defines the geometry, non-uniform mesh, and doping profile for the 1D p-n junction.

import devsim

def create_1d_pn_mesh():
    """
    Creates a 1D non-uniform mesh for a silicon p-n junction.
    Returns the device and region names for use in the simulator.
    """
    device_name = "pn_junction"
    region_name = "silicon"
    
    # 1. Initialize the 1D mesh
    devsim.create_1d_mesh(mesh="mesh1")
    
    # 2. Add mesh lines (Non-uniform spacing)
    # Total length: 2.0 um (from x = -1.0 to x = 1.0)
    # The analytical depletion width is ~0.43 um (-0.215 to +0.215)
    
    # Coarse spacing at the p-type contact (100 nm)
    devsim.add_1d_mesh_line(mesh="mesh1", pos=-1.0, ps=0.1, tag="top")
    
    # Transition to fine spacing at the edge of the depletion region (10 nm)
    devsim.add_1d_mesh_line(mesh="mesh1", pos=-0.25, ps=0.01, tag="p_depletion")
    
    # Extremely fine spacing exactly at the metallurgical junction (5 nm)
    devsim.add_1d_mesh_line(mesh="mesh1", pos=0.0, ps=0.005, tag="junction")
    
    # Transition back to 10 nm at the n-side depletion edge
    devsim.add_1d_mesh_line(mesh="mesh1", pos=0.25, ps=0.01, tag="n_depletion")
    
    # Coarse spacing at the n-type contact (100 nm)
    devsim.add_1d_mesh_line(mesh="mesh1", pos=1.0, ps=0.1, tag="bottom")
    
    # 3. Add regions and contacts
    devsim.add_1d_region(mesh="mesh1", material="Si", region=region_name, tag1="top", tag2="bottom")
    devsim.add_1d_contact(mesh="mesh1", name="top", tag="top", material="metal")
    devsim.add_1d_contact(mesh="mesh1", name="bottom", tag="bottom", material="metal")
    
    # 4. Finalize the mesh
    devsim.finalize_mesh(mesh="mesh1")
    devsim.create_device(mesh="mesh1", device=device_name)
    
    print(f"[{device_name}] Non-uniform mesh successfully generated.")
    
    return device_name, region_name

def apply_doping(device, region):
    """
    Applies the step-junction doping profile.
    p-side (x < 0): Acceptors (Na) = 1e16 cm^-3
    n-side (x >= 0): Donors (Nd) = 1e16 cm^-3
    """
    # Define spatial parameters
    devsim.set_parameter(device=device, region=region, name="Na", value=1e16)
    devsim.set_parameter(device=device, region=region, name="Nd", value=1e16)
    
    # Use a step function: if x < 0, doping is -Na (p-type). If x > 0, doping is +Nd (n-type).
    # DEVSIM uses 'NetDoping' as the standard node model variable.
    devsim.node_model(device=device, region=region, name="NetDoping", 
                      equation="ifelse(x < 0, -Na, Nd)")
    
    print(f"[{device}] Doping profile applied: Na = 1e16, Nd = 1e16.")

if __name__ == "__main__":
    # Test the setup script
    dev, reg = create_1d_pn_mesh()
    apply_doping(dev, reg)