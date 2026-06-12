% design_lcl_filter.m
% Implements parameter scripts for the grid-side LCL filter to attenuate 20kHz switching harmonics.

disp('Loading LCL Filter Parameters...');

% Base definitions for 20kHz LCL filter
lcl_L1 = 1.5e-3; % Converter-side inductor (H)
lcl_L2 = 0.5e-3; % Grid-side inductor (H)
lcl_Cf = 20e-6;  % Filter capacitor (F)
lcl_Rd = 0.5;    % Damping resistor (Ohms)

disp('LCL Filter Parameters loaded successfully.');
