% test_parameters.m
% Validation script for the CSSS simulation parameters

disp('Running tests on parameters...');

% Load parameters - use absolute or relative path so it can be run from the root
script_dir = fileparts(mfilename('fullpath'));
% Go one level up to reach the root folder
root_dir = fileparts(script_dir);

if isempty(root_dir)
    root_dir = pwd;
end

run(fullfile(root_dir, 'scripts', 'set_parameters.m'));

% Verify critical parameters
assert(exist('f_sw', 'var') == 1, 'f_sw is not defined');
assert(f_sw == 20000, 'f_sw should be 20000');

assert(exist('solver_MaxStep', 'var') == 1, 'solver_MaxStep is not defined');
assert(solver_MaxStep == 1 / (100 * f_sw), 'solver_MaxStep is incorrect');

assert(exist('solver_RelTol', 'var') == 1, 'solver_RelTol is not defined');
assert(solver_RelTol == 1e-4, 'solver_RelTol is incorrect');

assert(exist('solver_AbsTol', 'var') == 1, 'solver_AbsTol is not defined');
assert(solver_AbsTol == 1e-5, 'solver_AbsTol is incorrect');

assert(exist('solver_Type', 'var') == 1, 'solver_Type is not defined');
assert(strcmp(solver_Type, 'ode23tb'), 'solver_Type is incorrect');

disp('All tests passed successfully!');
