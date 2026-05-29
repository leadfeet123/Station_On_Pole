% test_parameters.m
% Test script to validate that set_parameters.m loads correctly and parameters are set as expected.

disp('Running parameter validation tests...');

% Run the setup script
base_dir = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(base_dir, 'scripts', 'set_parameters.m'));

% Verify critical parameters
assert(exist('f_sw', 'var') == 1, 'f_sw is not defined');
assert(f_sw == 20000, 'f_sw should be 20000');

assert(exist('solver_type', 'var') == 1, 'solver_type is not defined');
assert(strcmp(solver_type, 'ode23tb'), 'solver_type should be ode23tb');

assert(exist('rel_tol', 'var') == 1, 'rel_tol is not defined');
assert(rel_tol == 1e-4, 'rel_tol should be 1e-4');

assert(exist('abs_tol', 'var') == 1, 'abs_tol is not defined');
assert(abs_tol == 1e-5, 'abs_tol should be 1e-5');

assert(exist('max_step', 'var') == 1, 'max_step is not defined');
expected_max_step = 1 / (20000 * 100);
assert(abs(max_step - expected_max_step) < 1e-10, 'max_step value is incorrect');

disp('All parameter validation tests passed!');
