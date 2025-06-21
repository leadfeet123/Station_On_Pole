% initialize_CSSS_parameters.m
%
% Description:
% This script initializes all necessary parameters for the Compact Solid-State
% Substation (CSSS) simulation. Parameters are organized into structures and
% saved to 'CSSS_simulation_parameters.mat' for use by Simulink models
% and other scripts.
%
% As per NFR-01: All key component parameters shall be user-configurable
% from this single script.
%
% Version: 1.0
% Date:    YYYY-MM-DD % To be filled

disp('Initializing CSSS simulation parameters...');

% --- 1. Simulation Control Parameters ---
SimControl.StartTime = 0.0; % s
SimControl.StopTime = 0.2;  % s (short duration for initial testing)
SimControl.MaxStep = 1e-6;  % s (Max step size for variable step solver, suitable for power electronics)
SimControl.SolverType = 'ode23tb'; % Stiffness-aware solver, good for power electronics
SimControl.RelTol = 1e-4;
SimControl.AbsTol = 'auto'; % Or specify if needed, e.g., 1e-5 for currents, 1e-3 for voltages
SimControl.SimscapeLocalSolverSampleTime = SimControl.MaxStep / 10; % e.g. 1e-7 s, for Simscape Electrical solver

% --- 2. Grid Parameters (Input: 138kV) ---
Grid.V_LL_RMS = 138e3;    % V (Line-to-Line RMS)
Grid.V_LN_RMS = Grid.V_LL_RMS / sqrt(3); % V (Line-to-Neutral RMS)
Grid.Frequency_Hz = 60;   % Hz
Grid.Omega_rad_s = 2 * pi * Grid.Frequency_Hz; % rad/s
Grid.SourceInductance_H = 1e-3; % H (Example source inductance)
Grid.SourceResistance_Ohm = 0.1; % Ohm (Example source resistance)
% For three-phase source blocks in Simulink:
Grid.PhaseAngles_deg = [0, -120, 120]; % [deg] for phases A, B, C

% --- 3. Distribution Load Parameters ---
Load.V_LL_RMS_Nominal_kV = 11; % kV (Configurable, e.g., 11, 13.8, 25, 34.5)
Load.V_LL_RMS_Nominal = Load.V_LL_RMS_Nominal_kV * 1e3; % V
Load.V_LN_RMS_Nominal = Load.V_LL_RMS_Nominal / sqrt(3); % V
Load.PowerRating_MVA = 10;    % MVA (as per Model Specifications)
Load.PowerRating_VA = Load.PowerRating_MVA * 1e6; % VA
Load.PowerFactor = 0.95; % lagging (example)
Load.ActivePower_W = Load.PowerRating_VA * Load.PowerFactor; % W (Rated P)
Load.ReactivePower_VAR = Load.PowerRating_VA * sin(acos(Load.PowerFactor)); % VAR (Rated Q)
% For Simscape RLC Load block (assuming 3-phase balanced):
Load.Resistance_Ohm_per_phase = (Load.V_LN_RMS_Nominal^2) / (Load.ActivePower_W / 3);
Load.Inductance_H_per_phase = (Load.V_LN_RMS_Nominal^2) / (Load.ReactivePower_VAR / 3) / Grid.Omega_rad_s; % For inductive load
Load.Capacitance_F_per_phase = 0; % Assuming inductive load

% --- 4. AC-AC Converter Parameters ---
Converter.SwitchingFrequency_Hz = 20e3; % Hz (Configurable: 10kHz-50kHz)
Converter.Ts_sw = 1 / Converter.SwitchingFrequency_Hz; % s (Switching period)

% Input LCL Filter (example values, per phase)
Converter.L1_in_H = 0.5e-3;  % H
Converter.L2_in_H = 0.2e-3;  % H
Converter.C_in_F  = 5e-6;    % F
Converter.R_damp_in_Ohm = 0.5; % Ohm (Damping resistor for LCL filter)

% Output LC Filter (example values, per phase)
Converter.L_out_H = 0.3e-3;  % H
Converter.C_out_F = 10e-6;   % F
Converter.R_damp_out_Ohm = 0.3; % Ohm (Damping resistor for LC filter)

% Topology specific (e.g. if it's a DAB, DC link voltage might be a parameter)
Converter.Intermediate_DC_Link_Voltage_V = 2500; % V (Example, if topology has an internal DC link)

% --- 5. SiC MOSFET Parameters (`SimMOSFET` struct for component model) ---
SimMOSFET.Ron_Ohm = 0.012;      % Ohm (On-state resistance)
SimMOSFET.Vth_V = 2.2;          % V (Gate threshold voltage)
SimMOSFET.Cgs_F = 1.2e-9;       % F (Gate-Source capacitance)
SimMOSFET.Cgd_F = 0.3e-9;       % F (Gate-Drain capacitance)
SimMOSFET.Cds_F = 0.6e-9;       % F (Drain-Source capacitance)
SimMOSFET.Diode_Ron_Ohm = 0.015; % Ohm (Body diode on-resistance)
SimMOSFET.Diode_Vf_V = 0.9;     % V (Body diode forward voltage)
% Note: Actual switching times are emergent from these + gate drive.

% --- 6. High-Frequency Autotransformer Parameters (`SimAutotransformer` struct) ---
% Based on 138kV input, potential HF link voltage, and output voltage (e.g. 11kV)
% These are conceptual values for the HF stage.
% Assume a certain HF link RMS voltage, e.g., 2000V for primary side of HF transformer
SimAutotransformer.HF_Link_Voltage_RMS_Primary_Target = 2000; % V
SimAutotransformer.HF_Link_Voltage_RMS_Secondary_Target = 400; % V (Example step-down before unfolding to distribution)

% Autotransformer turns ratio calculation is complex for HF AC link.
% For now, placeholder values for N1 (series part) and N2 (common part).
% Voltage_HV_P1_to_Tap / N1 = Voltage_Tap_to_Common / N2
% Total primary voltage (HV_P1 to Common) related to (N1+N2)
% Output voltage (Tap to Common) related to N2
SimAutotransformer.N1_turns = 80; % Turns in the series winding part (non-common)
SimAutotransformer.N2_turns = 20; % Turns in the common winding part
% Effective primary turns = N1 + N2 = 100
% Effective secondary turns = N2 = 20
% Voltage ratio (Primary_Total / Secondary_Output) approx = (N1+N2)/N2 = 100/20 = 5:1

SimAutotransformer.R1_Ohm = 5e-3;    % Resistance of series winding (N1)
SimAutotransformer.Llk1_H = 2e-6;    % Leakage inductance of series winding (N1)
SimAutotransformer.R2_Ohm = 1e-3;    % Resistance of common winding (N2)
SimAutotransformer.Llk2_H = 0.5e-6;  % Leakage inductance of common winding (N2)

SimAutotransformer.CoreArea_m2 = 0.0025; % m^2 (Core cross-sectional area)
SimAutotransformer.CorePathLength_m = 0.35; % m (Mean magnetic path length)
% B-H Curve for Amorphous Alloy (Example: [H (A/m), B (T)])
SimAutotransformer.BH_Data = [0, 0; 20, 0.5; 50, 1.0; 100, 1.3; 200, 1.45; 500, 1.55; 1000, 1.58; 5000, 1.62];
% Steinmetz Coefficients for core loss (example for an amorphous material at specific freq)
SimAutotransformer.Kh_Steinmetz = 50.0; % Hysteresis coefficient
SimAutotransformer.alpha_Steinmetz = 1.6; % Hysteresis exponent
SimAutotransformer.beta_Steinmetz = 1.8;  % Hysteresis exponent
SimAutotransformer.Kc_Steinmetz = 0.5;  % Eddy current coefficient (classical losses)
% Note: Steinmetz params are frequency dependent. Values here are illustrative.

% --- 7. Control System Parameters (`SimControlParams` struct) ---
SimControlParams.ControlLoopSampleTime_s = Converter.Ts_sw / 2; % e.g., 25e-6 s for 20kHz switching

% Output Voltage Regulator (Outer Loop, controlling V_out_RMS or V_out_dq)
SimControlParams.V_out_RMS_Setpoint_V = Load.V_LN_RMS_Nominal; % V (Per-phase Line-to-Neutral RMS)
SimControlParams.VoltReg.Kp = 0.8;
SimControlParams.VoltReg.Ki = 15;
SimControlParams.VoltReg.Limit_Upper = 1.5; % Output limit for controller (e.g. current ref magnitude)
SimControlParams.VoltReg.Limit_Lower = -1.5;

% Input Power Factor Correction / Current Controller (Inner Loop, controlling I_in_dq)
SimControlParams.PF_Setpoint = 1.0; % Target unity power factor
SimControlParams.PF_Reg.I_d_ref_A = 50; % Placeholder, actual would come from outer voltage loop power demand
SimControlParams.PF_Reg.I_q_ref_A = 0;  % For unity PF, q-component of current ref is zero
SimControlParams.PF_Reg.Kp = 1.2;
SimControlParams.PF_Reg.Ki = 20;
SimControlParams.PF_Reg.Limit_Upper = 1.2; % Output limit for controller (e.g. modulation index)
SimControlParams.PF_Reg.Limit_Lower = -1.2;

% Fault Current Limiting (FCL)
SimControlParams.FCL.CurrentLimit_RMS_Factor = 1.5; % e.g., 1.5x rated current
SimControlParams.FCL.RatedCurrent_RMS_A = (Load.PowerRating_VA / 3) / Load.V_LN_RMS_Nominal; % A (per phase)
SimControlParams.FCL.CurrentLimit_RMS_A = SimControlParams.FCL.RatedCurrent_RMS_A * SimControlParams.FCL.CurrentLimit_RMS_Factor;
SimControlParams.FCL.ResponseTime_ms = 2; % ms (Target from FR-05)

% PLL Parameters for Grid Synchronization
SimControlParams.PLL.Kp = 100;
SimControlParams.PLL.Ki = 2000;
SimControlParams.PLL.Bandwidth_Hz = 10; % Approximate bandwidth of PLL

% --- 8. Save Parameters to MAT-file ---
% Group all structs into a single struct for saving, if desired, or save them individually.
AllParams.SimControl = SimControl;
AllParams.Grid = Grid;
AllParams.Load = Load;
AllParams.Converter = Converter;
AllParams.SimMOSFET = SimMOSFET;
AllParams.SimAutotransformer = SimAutotransformer;
AllParams.SimControlParams = SimControlParams;

save('CSSS_simulation_parameters.mat', 'AllParams', 'SimControl', 'Grid', 'Load', 'Converter', 'SimMOSFET', 'SimAutotransformer', 'SimControlParams');
disp('CSSS simulation parameters initialized and saved to CSSS_simulation_parameters.mat');

% --- End of Script ---
