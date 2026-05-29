% set_parameters.m
% Script to define base parameters for the Compact In-Line Solid-State Substation (CSSS)

disp('Loading parameters for CSSS simulation...');

% Grid Parameters
V_in_rms = 138000; % Input Voltage (V) - 138kV Transmission
f_grid = 60;       % Grid Frequency (Hz)

% SST Parameters
f_sw = 20000;      % Switching Frequency (Hz) - 20kHz for compact magnetics
% Note: Using 10kV SiC MOSFETs is assumed in the model design

% Solver Optimization Parameters
solver_type = 'ode23tb';      % Stiff solver for power electronics
rel_tol = 1e-4;               % Relative tolerance for 20kHz switching
abs_tol = 1e-5;               % Absolute tolerance
max_step = 1 / (f_sw * 100);  % Max step size (100 points per switching cycle)

disp('Parameters loaded successfully.');
