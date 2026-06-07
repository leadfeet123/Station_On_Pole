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

​- [ ] STEP 6: SIMSCAPE TOPOLOGY - Draft initial parameterization and connectivity scripts for the specific direct AC-AC matrix converter topology.

​- [ ] STEP 7: LCL FILTER DESIGN - Implement and define parameter scripts for the grid-side LCL filter to attenuate 20kHz switching harmonics.

​- [ ] STEP 8: SATURATION CURVE TUNING - Fine-tune the amorphous core autotransformer saturation curves based on the newly added validation scripts.

​- [ ] STEP 9: FAULT EXPANSION - Expand fault simulation parameters in the protection script to cover line-to-ground and three-phase faults.

​- [ ] STEP 10: PI CONTROLS - Enhance the PI voltage regulator loop controls by defining dynamic gains in set_parameters.m.

​Completed Iterations Log
​This section is managed by the Architect Agent. Completed checklists are archived here to maintain a history of the system's evolution.

​Iteration v1.1.0 (Refinement and Optimization)
​Accomplished baseline solver optimization, autotransformer validation, SiC converter loss analysis, FCL parameters auditing, and documentation refinement.
​[x] STEP 1: SOLVER OPTIMIZATION - Analyze scripts/set_parameters.m and write optimizations for numerical stability and solver tolerances to handle 20kHz switching.
​[x] STEP 2: AUTOTRANSFORMER VALIDATION - Expand tests/test_parameters.m or create a validation script for the non-linear amorphous core magnetizing inductance.
​[x] STEP 3: CONVERTER STRESS ANALYSIS - Review the switching frequencies and write a script calculating theoretical switching/conduction losses for the SiC stage.
​[x] STEP 4: PROTECTION & FCL AUDIT - Implement basic overcurrent threshold logic or parameter definitions in the scripts folder to prepare for fault current limiting.
​[x] STEP 5: DOCUMENTATION REFINE - Review all scripts for proper comment structures, unit mappings, and verify system compliance against the README.md PRD.

​Iteration v1.0.0 (Baseline Setup)
​[x] Initial structure established.
​[x] Simulation parameter specs defined.
​[x] Project README and PRD compiled.