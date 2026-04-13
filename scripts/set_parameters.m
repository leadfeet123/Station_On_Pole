% set_parameters.m
% Script to define base parameters for the Compact In-Line Solid-State Substation (CSSS)

disp('Loading parameters for CSSS simulation...');

% Grid Parameters
V_in_rms = 138000; % Input Voltage (V) - 138kV Transmission
f_grid = 60;       % Grid Frequency (Hz)

% SST Parameters
f_sw = 15000;      % Switching Frequency (Hz) - >10kHz for compact magnetics
% Note: Using 10kV SiC MOSFETs is assumed in the model design

disp('Parameters loaded successfully.');
