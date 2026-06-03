% set_parameters.m
% Script to define base parameters for the Compact In-Line Solid-State Substation (CSSS)

disp('Loading parameters for CSSS simulation...');

% Grid Parameters
V_in_rms = 138000; % Input Voltage (V) - 138kV Transmission
f_grid = 60;       % Grid Frequency (Hz)

% SST Parameters
f_sw = 20000;      % Switching Frequency (Hz) - >10kHz for compact magnetics
% Note: Using 10kV SiC MOSFETs is assumed in the model design

% Solver Optimizations for 20kHz switching
solver_MaxStep = 1 / (100 * f_sw); % High resolution for 20kHz switching
solver_RelTol = 1e-4;              % Tighter relative tolerance
solver_AbsTol = 1e-5;              % Tighter absolute tolerance
solver_Type = 'ode23tb';           % Stiff solver suited for power electronics

disp('Parameters loaded successfully.');
