# Compact In-Line Solid-State Substation (CSSS)

Welcome to the **Compact In-Line Solid-State Substation (CSSS)** project. This repository contains the architecture, documentation, and simulation models for a novel single-stage, non-isolated Solid-State Transformer (SST).

## Project Overview
The CSSS aims to revolutionize power distribution by converting transmission-level AC voltage (e.g., 138kV) to standard distribution feeder voltages through a compact, in-line solid-state substation. By utilizing advanced 10kV Silicon Carbide (SiC) MOSFETs operating at high switching frequencies (>10kHz) and a non-isolated autotransformer with an amorphous magnetic core, this design minimizes size and conversion losses while eliminating the need for galvanic isolation.

The project models serve as a virtual prototype to simulate, analyze, and validate the performance, control strategies (such as voltage regulation, power factor correction, and fault current limiting), and operational limits of this architecture.

## Repository Structure

- `docs/`: Contains project documentation.
  - `PRD.md`: The Product Requirements Document outlining the specifications and goals of the CSSS model.
  - `RESEARCH.md`: A summary of updated research on Solid-State Transformers and 10kV Silicon Carbide (SiC) technology.
- `simulink/`: Contains the Simulink models for the CSSS.
  - `main_sst_model.slx`: The primary Simulink model for the system.
  - `components/`: Sub-components used within the main model.
- `scripts/`: Contains Octave/MATLAB scripts for managing the simulation.
  - `set_parameters.m`: Script to define and load the simulation parameters (voltages, frequencies, component properties) into the workspace.
  - `run_simulation.m`: Script to initialize parameters and execute the Simulink model.

## Usage

To run the simulations, you will need MATLAB/Simulink (R2024a or later) with the Simscape Electrical™ toolbox, or a compatible environment like GNU Octave for executing the setup scripts.

1.  **Initialize Parameters**: Run `set_parameters.m` to load the base parameters for the 138kV system and high-frequency SiC components into your workspace.
2.  **Run Simulation**: Execute `run_simulation.m` to trigger the simulation and begin analysis.

Please refer to `docs/PRD.md` and `docs/RESEARCH.md` for deeper insights into the design requirements and the underlying technologies driving the CSSS architecture.
