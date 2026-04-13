% run_simulation.m
% Script to initialize and run the CSSS simulation

disp('Initializing simulation environment...');

% Load parameters
run('set_parameters.m');

% In a real MATLAB environment, this would start the Simulink model:
% sim('simulink/main_sst_model.slx');
disp('Starting Simulink model execution: simulink/main_sst_model.slx');
disp('Simulation complete.');
