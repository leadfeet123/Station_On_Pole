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

​[x] STEP 1: SOLVER OPTIMIZATION - Analyze scripts/set_parameters.m and write optimizations for numerical stability and solver tolerances to handle 20kHz switching.

​[x] STEP 2: AUTOTRANSFORMER VALIDATION - Expand tests/test_parameters.m or create a validation script for the non-linear amorphous core magnetizing inductance.

​[ ] STEP 3: CONVERTER STRESS ANALYSIS - Review the switching frequencies and write a script calculating theoretical switching/conduction losses for the SiC stage.

​[ ] STEP 4: PROTECTION & FCL AUDIT - Implement basic overcurrent threshold logic or parameter definitions in the scripts folder to prepare for fault current limiting.

​[ ] STEP 5: DOCUMENTATION REFINE - Review all scripts for proper comment structures, unit mappings, and verify system compliance against the README.md PRD.

​Completed Iterations Log
​This section is managed by the Architect Agent. Completed checklists are archived here to maintain a history of the system's evolution.

​Iteration v1.0.0 (Baseline Setup)
​[x] Initial structure established.
​[x] Simulation parameter specs defined.
​[x] Project README and PRD compiled.