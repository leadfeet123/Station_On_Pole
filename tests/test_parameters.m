% test_parameters.m
% Script to test parameter generation and solver optimizations

% Find the scripts directory relative to the tests directory
script_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'scripts');
% Construct absolute path for Octave run compatibility
run(fullfile(script_dir, 'set_parameters.m'));

% Assert f_sw is 20000
assert(exist('f_sw', 'var') == 1, 'f_sw is not defined');
assert(f_sw == 20000, 'f_sw is not equal to 20000');

% Assert solver optimization parameters exist
assert(exist('MaxStep', 'var') == 1, 'MaxStep is not defined');
assert(exist('RelTol', 'var') == 1, 'RelTol is not defined');
assert(exist('AbsTol', 'var') == 1, 'AbsTol is not defined');

% Verify MaxStep is small enough for 20kHz
% Rule of thumb: at least 100 points per switching period -> MaxStep <= 1 / (20000 * 100)
assert(MaxStep <= 1 / (20000 * 100), 'MaxStep is too large for 20kHz switching stability');

disp('All parameter tests passed successfully.');
