This repository contains a MATLAB/Simulink project for a Compact In-Line Solid-State Substation (CSSS) designed to convert 138kV transmission to 34.5kV distribution.
​System Configuration
​High-Side: 138kV AC, 3-Phase, 60 Hz
​Low-Side: 34.5kV AC, 3-Phase, 60 Hz
​Power Rating: 20 MVA nominal
​Transformer: High-frequency non-isolated autotransformer with an amorphous alloy core (Bsat = 1.56 T)
​Converter: Single-stage direct AC-AC converter built with Silicon Carbide (SiC) switches
​Development Strategy
​Binary .slx files are difficult to version. We use programmatic MATLAB scripts (.m files) in /scripts to generate, parameterize, and test system models.
​Continuous Improvement Checklist
​The Coder Agent must execute exactly one unchecked step per daily run, mark it as completed, and submit a PR.

​[ ] STEP 6: DIRECT AC-AC CONVERTER MODELING - Draft an initial Simscape implementation script for the direct AC-AC matrix converter topology, mapping the 10kV SiC switches to the appropriate electrical blocks.

​[ ] STEP 7: AMORPHOUS CORE SATURATION FINE-TUNING - Expand set_parameters.m to include detailed piece-wise linear saturation curve data points for the amorphous core, accurately reflecting the B-H curve up to Bsat = 1.56 T.

​[ ] STEP 8: ADVANCED FAULT SCENARIO PARAMETERS - Introduce new fault simulation parameters in the scripts folder to model asymmetric grid faults (e.g., single-line-to-ground) and their impact on the non-isolated substation.

​[ ] STEP 9: VOLTAGE REGULATION PI TUNING - Create a control initialization script init_vreg_controls.m to define the PI controller gains (Kp, Ki) and anti-windup limits for the output voltage regulation loop.

​[ ] STEP 10: HARMONIC FILTER DESIGN - Add LCL filter sizing calculations to the parameters script, ensuring the direct AC-AC converter meets the <3% THD requirement at the 34.5kV output.

​Completed Iterations Log
​This section is managed by the Architect Agent. Completed checklists are archived here to maintain a history of the system's evolution.

​Iteration v1.1.0 (Continuous Improvement)
​Advanced optimization, component validation, and documentation auditing were performed to enhance simulation robustness.
​[x] STEP 1: SOLVER OPTIMIZATION - Analyze scripts/set_parameters.m and write optimizations for numerical stability and solver tolerances to handle 20kHz switching.
​[x] STEP 2: AUTOTRANSFORMER VALIDATION - Expand tests/test_parameters.m or create a validation script for the non-linear amorphous core magnetizing inductance.
​[x] STEP 3: CONVERTER STRESS ANALYSIS - Review the switching frequencies and write a script calculating theoretical switching/conduction losses for the SiC stage.
​[x] STEP 4: PROTECTION & FCL AUDIT - Implement basic overcurrent threshold logic or parameter definitions in the scripts folder to prepare for fault current limiting.
​[x] STEP 5: DOCUMENTATION REFINE - Review all scripts for proper comment structures, unit mappings, and verify system compliance against the README.md PRD.

​Iteration v1.0.0 (Baseline Setup)
​[x] Initial structure established.
​[x] Simulation parameter specs defined.
​[x] Project README and PRD compiled.