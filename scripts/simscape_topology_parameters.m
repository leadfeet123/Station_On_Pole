% simscape_topology_parameters.m
% Initial parameterization and connectivity scripts for the specific direct AC-AC matrix converter topology.

disp('Loading Simscape Topology Parameters for Matrix Converter...');

mc_num_phases_in = 3;
mc_num_phases_out = 3;
mc_efficiency_target = 0.98;
mc_switching_freq = 20000; % 20kHz

% Aliases for backward compatibility
num_input_phases = mc_num_phases_in;
num_output_phases = mc_num_phases_out;

disp('Simscape Topology Parameters loaded successfully.');
