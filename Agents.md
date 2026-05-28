This is a MATLAB/Simulink project for a Compact In-Line Solid-State Substation (CSSS).
​Environment: MATLAB R2024a+ / Simscape Electrical.
​Goal: Build a single-stage, non-isolated direct AC-AC converter with an amorphous core autotransformer (138kV to 34.5kV).
​Development Strategy: Since binary .slx files are hard to modify directly via Git, we use MATLAB script files (.m files) to programmatically build, parameterize, and modify Simulink blocks (using add_block, set_param, and Simulink.Bus scripts) alongside standard Simulink library models.
​Directory Structure:
​/models for Simulink .slx models.
​/scripts for .m parameter and programmatic model-building scripts.
​/data for .sldd data dictionaries.