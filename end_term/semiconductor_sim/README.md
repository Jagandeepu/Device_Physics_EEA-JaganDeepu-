# TCAD Simulation of Silicon p-n Junctions (1D & 2D)

## Part 1 — What did you build?

This project implements a complete TCAD (Technology Computer-Aided Design)  to simulate the electrostatics and transport physics of a Silicon p-n junction. I used the **DEVSIM** numerical solver to implement a 1D Drift-Diffusion model and a 2D lateral device model. I extracted key physical quantities including electrostatic potential profiles, electron/hole carrier densities, and full I-V characteristics spanning forward and reverse bias.

---

## Part 2 — How to set it up

Clone the repository and install the necessary Python libraries. It is recommended to use a virtual environment.

```bash
git clone https://github.com/Jagandeepu/Device_Physics_EEA-JaganDeepu-.git
cd Device_Physics_EEA-JaganDeepu-
pip install -r requirements.txt
```

---

## Part 3 — How to run

The project is modularized into specific stages. Execute them in the following order:

```bash
python3 verify.py        # Compares analytical theory vs numerical targets
python3 device_setup.py  # Creates the non-uniform mesh and doping profile
python3 simulate.py      # Solves Drift-Diffusion, runs voltage sweep, saves data
python3 visualize.py     # Generates the 4 required 1D plots in /results
python3 simulate_2d.py   # (Bonus) Executes 2D lateral device simulation
python3 visualize_2d.py  # (Bonus) Generates 2D current vector streamplots
```

---

## Part 4 — Your results

The following table validates the DEVSIM numerical results against the analytical physics calculated in `verify.py` for a Silicon junction at 300 K (Na = Nd = 10^16 cm⁻³):

| Quantity | Theory (Analytical) | Simulation (DEVSIM) |
|---|---|---|
| Built-in Potential (V) | 0.7144 V | 0.7144 V |
| Depletion Width (W at V=0) | 0.4299 µm | 0.4300 µm |
| Forward Current (V=0.6V) | 8.271e-03 A/cm² | 8.279e-03 A/cm² |
| Forward Current (V=0.8V) | 4.113e-02 A/cm² | 4.203e-02 A/cm² |

### Physical Observations

- **Barrier Modulation:** As seen in Plot 1, the potential barrier height shrinks under forward bias and expands under reverse bias.
- **Series Resistance:** At high injection (0.8V), the simulation correctly deviates from the ideal Shockley equation due to bulk series resistance, which is a more realistic physical model.

---

## Part 5 — Known limitations

This simulation utilizes a simplified "Step Junction" profile and assumes a constant temperature of 300 K. It currently ignores SRH (Shockley-Read-Hall) and Auger recombination models, which would impact the leakage current in the reverse bias regime. Furthermore, it does not account for High-Field Velocity Saturation. With more time, I would implement a self-consistent Heat Equation solver to model thermal rollover.

---

## Part 6 — Bonus: 2D Vector Fields

In addition to the 1D model, I simulated a 2D lateral diode. The results in `results_2d/` visualize the Total Current Density Vector Field. The streamplots demonstrate how carriers dive into the silicon bulk to bypass surface bottlenecks, a phenomenon that purely 1D models cannot capture.