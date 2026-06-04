% analyze_converter_stress.m
% Estimate theoretical SiC converter switching and conduction losses.

disp('Running converter stress analysis...');

script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end
run(fullfile(script_dir, 'set_parameters.m'));

% User-overridable assumptions
if ~exist('S_nominal_va', 'var')
    S_nominal_va = 10e6; % VA
end
if ~exist('V_ll_nominal', 'var')
    V_ll_nominal = 34.5e3; % V
end
if ~exist('power_factor', 'var')
    power_factor = 0.99; % -
end
if ~exist('modulation_index', 'var')
    modulation_index = 0.9; % -
end
if ~exist('num_switches', 'var')
    num_switches = 12; % 3-phase, 2-level equivalent estimate
end
if ~exist('V_device', 'var')
    V_device = 10e3; % V
end
if ~exist('E_sw_ref_j', 'var')
    E_sw_ref_j = 0.15; % J at reference conditions, per event
end
if ~exist('I_ref_a', 'var')
    I_ref_a = 300; % A reference current for E_sw_ref_j
end
if ~exist('alpha_sw', 'var')
    alpha_sw = 1.1; % switching energy current scaling exponent
end
if ~exist('R_on', 'var')
    R_on = 0.012; % Ohm, from set_parameters.m
end

% Derived electrical quantities
I_line_rms = S_nominal_va / (sqrt(3) * V_ll_nominal);
I_switch_rms = I_line_rms / sqrt(2); % each switch conducts about half-cycle
I_switch_avg = (2 / pi) * I_line_rms / 2; % rough sinusoidal average while on

% Loss calculations
E_sw_scaled = E_sw_ref_j * (max(I_line_rms, 1e-9) / I_ref_a)^alpha_sw;
P_sw_per_switch = 2 * E_sw_scaled * f_sw; % turn-on + turn-off
P_cond_per_switch = (I_switch_rms^2) * R_on;

P_sw_total = num_switches * P_sw_per_switch;
P_cond_total = num_switches * P_cond_per_switch;
P_total = P_sw_total + P_cond_total;

% Normalized converter stress indicators
apparent_output_power = sqrt(3) * V_ll_nominal * I_line_rms;
loss_fraction = P_total / max(apparent_output_power, 1e-9);
estimated_efficiency = max(0, 1 - loss_fraction);

fprintf('\n--- Converter Stress Summary ---\n');
fprintf('Switching frequency f_sw:           %.0f Hz\n', f_sw);
fprintf('Estimated line current (RMS):       %.2f A\n', I_line_rms);
fprintf('Per-switch switching loss:          %.2f W\n', P_sw_per_switch);
fprintf('Per-switch conduction loss:         %.2f W\n', P_cond_per_switch);
fprintf('Total switching loss:               %.2f W\n', P_sw_total);
fprintf('Total conduction loss:              %.2f W\n', P_cond_total);
fprintf('Total converter device loss:        %.2f W\n', P_total);
fprintf('Estimated efficiency (device-only): %.4f pu\n', estimated_efficiency);
fprintf('Assumed power factor:               %.2f\n', power_factor);
fprintf('Assumed modulation index:           %.2f\n', modulation_index);
fprintf('--------------------------------\n\n');

StressResults = struct( ...
    'f_sw_hz', f_sw, ...
    'S_nominal_va', S_nominal_va, ...
    'V_ll_nominal_v', V_ll_nominal, ...
    'line_current_rms_a', I_line_rms, ...
    'switching_loss_per_switch_w', P_sw_per_switch, ...
    'conduction_loss_per_switch_w', P_cond_per_switch, ...
    'switching_loss_total_w', P_sw_total, ...
    'conduction_loss_total_w', P_cond_total, ...
    'total_loss_w', P_total, ...
    'estimated_efficiency_pu', estimated_efficiency, ...
    'num_switches', num_switches);

save(fullfile(script_dir, 'converter_stress_results.mat'), 'StressResults');
disp('Saved results to scripts/converter_stress_results.mat');
