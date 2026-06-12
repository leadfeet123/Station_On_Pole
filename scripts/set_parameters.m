% set_parameters.m
% Script to define base parameters for the Compact In-Line Solid-State Substation (CSSS)

disp('Loading parameters for CSSS simulation...');

% Grid Parameters
V_in_rms = 138000; % Input Voltage (V) - 138kV Transmission
f_grid = 60;       % Grid Frequency (Hz)

% SST Parameters
f_sw = 20000;      % Switching Frequency (Hz) - 20kHz for compact magnetics
% Note: Using 10kV SiC MOSFETs is assumed in the model design

% Unit mapping for key variables:
%   *_rms: RMS electrical quantities (SI base units)
%   *_va:  Apparent power (VA)
%   *_ms:  Time in milliseconds
%   *_s:   Time in seconds
%   *_pu:  Per-unit ratios

% Solver Optimizations for 20kHz switching
solver_MaxStep = 1 / (100 * f_sw); % High resolution for 20kHz switching
solver_RelTol = 1e-4;              % Tighter relative tolerance
solver_AbsTol = 1e-5;              % Tighter absolute tolerance
solver_Type = 'ode23tb';           % Stiff solver suited for power electronics

% Autotransformer Amorphous Core Parameters
B_sat = 1.56;        % Saturation flux density (T)
L_mag = 0.05;        % Nominal magnetizing inductance (H) for validation

% Sample time (s) - recommended for high-freq switching
T_s = 1e-6;

% Protection & Fault Current Limiting (FCL) Parameters
S_rated_va = 20e6;        % Nominal three-phase apparent power rating (VA)
V_out_ll_rms = 34500;     % Nominal low-side line-to-line RMS voltage (V)
I_rated_rms = S_rated_va / (sqrt(3) * V_out_ll_rms); % Rated line current (A RMS)
pf_target = 0.99;         % Input-side target power factor (FR-04)

fcl_limit_factor = 1.5;   % Fault current limit multiplier (pu of rated current)
fcl_current_limit_rms = fcl_limit_factor * I_rated_rms; % Current limit threshold (A RMS)
fcl_trip_threshold_rms = fcl_current_limit_rms;         % Overcurrent trip threshold (A RMS)
fcl_release_factor = 0.95;                              % Reset threshold as fraction of trip level
fcl_release_threshold_rms = fcl_release_factor * fcl_trip_threshold_rms;
fcl_response_ms = 2;                                    % Target limit response time (ms)
fcl_response_s = fcl_response_ms / 1000;                % Response time in seconds
fcl_detection_samples = ceil(fcl_response_s / T_s);     % Sample window for overcurrent confirmation

% Backward-compatible aliases used by alternate scripts/branches
RelTol = solver_RelTol;
AbsTol = solver_AbsTol;
MaxStep = solver_MaxStep;

solver_type = solver_Type;
rel_tol = solver_RelTol;
abs_tol = solver_AbsTol;
max_step = solver_MaxStep;
solver_max_step = solver_MaxStep;
solver_rel_tol = solver_RelTol;
solver_abs_tol = solver_AbsTol;
Ts = solver_MaxStep;

fcl_limit_pu = fcl_limit_factor;
fcl_current_limit_a = fcl_current_limit_rms;
fcl_trip_threshold_a = fcl_trip_threshold_rms;
fcl_release_threshold_a = fcl_release_threshold_rms;
pf_setpoint = pf_target;

run(fullfile(fileparts(mfilename('fullpath')), 'simscape_topology_parameters.m'));
run(fullfile(fileparts(mfilename('fullpath')), 'design_lcl_filter.m'));

disp('Parameters loaded successfully.');
