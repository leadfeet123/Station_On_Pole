% design_lcl_filter.m
% Parameter script for the grid-side LCL filter to attenuate 20kHz switching harmonics.

disp('Designing LCL filter...');

% Ensure parameters are loaded
if ~exist('V_in_rms', 'var') || ~exist('f_sw', 'var')
    script_dir = fileparts(mfilename('fullpath'));
    run(fullfile(script_dir, 'set_parameters.m'));
end

% Grid-side LCL filter parameters (Attenuate 20kHz switching frequency)
L_f1 = 1e-3; % Converter-side inductor (H)
L_f2 = 0.5e-3; % Grid-side inductor (H)
C_f = 10e-6; % Filter capacitor (F)
R_d = 0.5; % Damping resistor (Ohms)

disp('LCL filter parameters defined.');
