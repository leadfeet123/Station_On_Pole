% set_parameters.m
% Script to define base parameters for the Compact In-Line Solid-State Substation (CSSS)

disp('Loading parameters for CSSS simulation...');

% Grid Parameters
V_in_rms = 138000; % Input Voltage (V) - 138kV Transmission
f_grid = 60;       % Grid Frequency (Hz)

% SST Parameters
f_sw = 20000;      % Switching Frequency (Hz) - 20kHz for compact magnetics
% Note: Using 10kV SiC MOSFETs is assumed in the model design

% Solver Optimizations for 20kHz switching
solver_MaxStep = 1 / (100 * f_sw); % High resolution for 20kHz switching
solver_RelTol = 1e-4;              % Tighter relative tolerance
solver_AbsTol = 1e-5;              % Tighter absolute tolerance
solver_Type = 'ode23tb';           % Stiff solver suited for power electronics

% Autotransformer Amorphous Core Parameters
B_sat = 1.56;        % Saturation flux density (T)
L_mag = 0.05;        % Nominal magnetizing inductance (H) for validation

% Backward-compatible aliases used by alternate scripts/branches
T_s = 1e-6;        % Sample time (s) - recommended for high-freq switching
RelTol = solver_RelTol;
AbsTol = solver_AbsTol;
MaxStep = solver_MaxStep;

solver_type = solver_Type;
rel_tol = solver_RelTol;
abs_tol = solver_AbsTol;
max_step = solver_MaxStep;
solver_max_step = solver_MaxStep;
solver_rel_tol = solver_RelTol;
solver_abs_tol = solver_AbsTol;
Ts = solver_MaxStep;

disp('Parameters loaded successfully.');
