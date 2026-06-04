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
assert(exist('V_in_rms', 'var') == 1, 'V_in_rms is not defined');
assert(V_in_rms == 138000, 'V_in_rms should be 138000');

assert(exist('f_grid', 'var') == 1, 'f_grid is not defined');
assert(f_grid == 60, 'f_grid should be 60');

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

assert(exist('solver_type', 'var') == 1, 'solver_type is not defined');
assert(strcmp(solver_type, 'ode23tb'), 'solver_type should be ode23tb');

assert(exist('solver_max_step', 'var') == 1, 'solver_max_step is not defined');
assert(solver_max_step == 1 / (100 * f_sw), 'solver_max_step value is incorrect');

assert(exist('solver_rel_tol', 'var') == 1, 'solver_rel_tol is not defined');
assert(solver_rel_tol == 1e-4, 'solver_rel_tol should be 1e-4');

assert(exist('solver_abs_tol', 'var') == 1, 'solver_abs_tol is not defined');
assert(solver_abs_tol == 1e-5, 'solver_abs_tol should be 1e-5');

assert(exist('rel_tol', 'var') == 1, 'rel_tol is not defined');
assert(rel_tol == 1e-4, 'rel_tol should be 1e-4');

assert(exist('abs_tol', 'var') == 1, 'abs_tol is not defined');
assert(abs_tol == 1e-5, 'abs_tol should be 1e-5');

assert(exist('max_step', 'var') == 1, 'max_step is not defined');
expected_max_step = 1 / (20000 * 100);
assert(abs(max_step - expected_max_step) < 1e-10, 'max_step value is incorrect');

assert(exist('Ts', 'var') == 1, 'Ts should be defined');
assert(Ts == expected_max_step, 'Ts value is incorrect');

assert(exist('B_sat', 'var') == 1, 'B_sat should be defined');
assert(B_sat == 1.56, 'B_sat value is incorrect');

assert(exist('L_mag', 'var') == 1, 'L_mag should be defined');
assert(L_mag == 0.05, 'L_mag value is incorrect');

assert(exist('S_rated_va', 'var') == 1, 'S_rated_va should be defined');
assert(S_rated_va == 20e6, 'S_rated_va value is incorrect');

assert(exist('V_out_ll_rms', 'var') == 1, 'V_out_ll_rms should be defined');
assert(V_out_ll_rms == 34500, 'V_out_ll_rms value is incorrect');

assert(exist('I_rated_rms', 'var') == 1, 'I_rated_rms should be defined');
assert(I_rated_rms > 0, 'I_rated_rms should be positive');

assert(exist('pf_target', 'var') == 1, 'pf_target should be defined');
assert(pf_target == 0.99, 'pf_target should be 0.99');

assert(exist('fcl_limit_factor', 'var') == 1, 'fcl_limit_factor should be defined');
assert(fcl_limit_factor == 1.5, 'fcl_limit_factor value is incorrect');

assert(exist('fcl_trip_threshold_rms', 'var') == 1, 'fcl_trip_threshold_rms should be defined');
assert(fcl_trip_threshold_rms == fcl_limit_factor * I_rated_rms, 'fcl_trip_threshold_rms value is incorrect');

assert(exist('fcl_release_threshold_rms', 'var') == 1, 'fcl_release_threshold_rms should be defined');
assert(fcl_release_threshold_rms < fcl_trip_threshold_rms, 'fcl_release_threshold_rms should be below trip threshold');

assert(exist('fcl_response_ms', 'var') == 1, 'fcl_response_ms should be defined');
assert(fcl_response_ms == 2, 'fcl_response_ms should be 2 ms');

assert(exist('fcl_detection_samples', 'var') == 1, 'fcl_detection_samples should be defined');
assert(fcl_detection_samples == ceil((fcl_response_ms / 1000) / T_s), 'fcl_detection_samples value is incorrect');

assert(exist('fcl_current_limit_a', 'var') == 1, 'fcl_current_limit_a should be defined');
assert(fcl_current_limit_a == fcl_current_limit_rms, 'fcl_current_limit_a alias mismatch');

disp('All parameter validation tests passed!');
