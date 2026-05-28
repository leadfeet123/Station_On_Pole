% set_parameters.m
% Script to define base parameters for the Compact In-Line Solid-State Substation (CSSS)

disp('Loading parameters for CSSS simulation...');

% Grid Parameters
V_in_rms = 138000; % Input Voltage (V) - 138kV Transmission
f_grid = 60;       % Grid Frequency (Hz)

% SST Parameters
f_sw = 20000;      % Switching Frequency (Hz) - >10kHz for compact magnetics, updated to 20kHz
% Note: Using 10kV SiC MOSFETs is assumed in the model design

% Solver Optimization Parameters for 20kHz switching
% Using rule of thumb: 100 points per switching period for stability
MaxStep = 1 / (f_sw * 100);
RelTol = 1e-4;
AbsTol = 1e-4;

disp('Parameters loaded successfully.');
