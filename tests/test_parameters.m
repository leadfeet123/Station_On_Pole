% tests/test_parameters.m
disp('Running test_parameters.m...');
script_dir = fullfile(fileparts(pwd), 'scripts');
% Check if we are running from project root
if exist('scripts/set_parameters.m', 'file')
    script_dir = fullfile(pwd, 'scripts');
end
run(fullfile(script_dir, 'set_parameters.m'));

assert(f_sw == 20000, 'f_sw should be 20000');
assert(exist('Ts', 'var') == 1, 'Ts should be defined');
assert(strcmp(solver_type, 'ode23tb'), 'solver_type should be ode23tb');
assert(rel_tol == 1e-4, 'rel_tol should be 1e-4');
assert(abs_tol == 1e-5, 'abs_tol should be 1e-5');
disp('All tests passed.');
