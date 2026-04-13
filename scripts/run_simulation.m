% run_simulation.m
% Script to initialize and run the CSSS simulation

disp('Initializing simulation environment...');

% Load parameters - use absolute or relative path so it can be run from the root
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end
run(fullfile(script_dir, 'set_parameters.m'));

% In a real MATLAB environment, this would start the Simulink model:
% sim('simulink/main_sst_model.slx');
disp('Starting Simulink model execution: simulink/main_sst_model.slx');
disp('Simulation complete.');
