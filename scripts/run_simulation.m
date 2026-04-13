% run_simulation.m
% Script to initialize and run the CSSS simulation

disp('Initializing simulation environment...');

% Load parameters - use absolute or relative path so it can be run from the root
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end
run(fullfile(script_dir, 'set_parameters.m'));

% To run the Simulink model, uncomment the line below in MATLAB/Simulink:
% sim(fullfile(fileparts(script_dir), 'simulink', 'main_sst_model.slx'));
disp('Parameters loaded. To run the simulation, open MATLAB/Simulink and call sim() on simulink/main_sst_model.slx.');
