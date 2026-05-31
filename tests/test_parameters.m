% test_parameters.m
% Test script to validate parameter definitions

disp('Running tests for parameter definitions...');

% Run parameter setup script
% Calculate absolute path for set_parameters.m
[test_dir, ~, ~] = fileparts(mfilename('fullpath'));
if isempty(test_dir)
    test_dir = pwd;
end
script_dir = fullfile(fileparts(test_dir), 'scripts');
run(fullfile(script_dir, 'set_parameters.m'));

% Test variables exist and have correct values
assert(exist('V_in_rms', 'var') == 1, 'V_in_rms is not defined');
assert(V_in_rms == 138000, 'V_in_rms should be 138000');

assert(exist('f_grid', 'var') == 1, 'f_grid is not defined');
assert(f_grid == 60, 'f_grid should be 60');

assert(exist('f_sw', 'var') == 1, 'f_sw is not defined');
assert(f_sw == 15000, 'f_sw should be 15000');

% Test solver parameters
assert(exist('solver_type', 'var') == 1, 'solver_type is not defined');
assert(strcmp(solver_type, 'ode23tb'), 'solver_type should be ode23tb');

assert(exist('solver_max_step', 'var') == 1, 'solver_max_step is not defined');
assert(solver_max_step <= 1/(20 * f_sw), 'solver_max_step should be small enough to capture 20kHz switching');

assert(exist('solver_rel_tol', 'var') == 1, 'solver_rel_tol is not defined');
assert(solver_rel_tol <= 1e-4, 'solver_rel_tol should be adequately tight for stability');

disp('All parameter tests passed successfully.');
