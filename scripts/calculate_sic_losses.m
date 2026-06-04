% calculate_sic_losses.m
% Script to calculate theoretical switching and conduction losses for the SiC stage

disp('Calculating theoretical switching and conduction losses for SiC stage...');

% Ensure parameters are loaded
if ~exist('V_in_rms', 'var') || ~exist('f_sw', 'var')
    % Resolve path dynamically
    script_dir = fileparts(mfilename('fullpath'));
    run(fullfile(script_dir, 'set_parameters.m'));
end

% SiC MOSFET properties (assumed 10kV class)
R_ds_on = 0.05;         % On-state resistance (Ohms)
E_on = 0.015;           % Turn-on energy per pulse (J) at nominal conditions
E_off = 0.012;          % Turn-off energy per pulse (J) at nominal conditions
V_nom = 138000;         % Nominal voltage for Energy numbers
I_nom = 100;            % Nominal current for Energy numbers

% System Power
S_nominal = 20e6;       % 20 MVA nominal
I_rms = S_nominal / (sqrt(3) * V_in_rms);
I_peak = I_rms * sqrt(2);

% Theoretical Conduction Losses per switch
P_cond = I_rms^2 * R_ds_on;

% Theoretical Switching Losses per switch (scaled by actual current)
% Average current over a half sine wave is (2/pi) * I_peak
I_avg_sw = (2/pi) * I_peak;

% Scale energy by ratio of actual current to nominal current (simplified linear model)
E_on_actual = E_on * (I_avg_sw / I_nom);
E_off_actual = E_off * (I_avg_sw / I_nom);

% Switching losses
P_sw = (E_on_actual + E_off_actual) * f_sw;

% Total losses per switch
P_total_switch = P_cond + P_sw;

disp(['Conduction losses per switch: ', num2str(P_cond), ' W']);
disp(['Switching losses per switch: ', num2str(P_sw), ' W']);
disp(['Total theoretical losses per switch: ', num2str(P_total_switch), ' W']);
