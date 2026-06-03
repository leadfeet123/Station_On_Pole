% post_process_CSSS_results.m
%
% Description:
% This script loads and processes the simulation output data from the CSSS
% simulation ('CSSS_simulation_output.mat'). It generates various plots to
% visualize key electrical quantities and performance metrics.
%
% Version: 1.0
% Date:    YYYY-MM-DD % To be filled

disp('Starting CSSS post-processing script...');

% --- Configuration ---
outputDataFile = 'CSSS_simulation_output.mat'; % Data file from run_CSSS_simulation.m
plotsSubFolder = 'CSSS_simulation_plots';     % Folder to save plots

% --- 1. Load Simulation Data ---
if exist(outputDataFile, 'file')
    disp(['Loading simulation output from ' outputDataFile '...']);
    load(outputDataFile); % Loads simOut, AllParams, SimControl, Grid, etc.
    disp('Simulation data loaded successfully.');
else
    disp(['Error: Simulation output file (' outputDataFile ') not found.']);
    disp('Please run the simulation first using run_CSSS_simulation.m.');
    return;
end

% Check if simOut exists (it should if the file loaded correctly)
if ~exist('simOut', 'var') || isempty(simOut)
    disp('Error: simOut variable not found or is empty in the loaded data.');
    disp('Ensure the simulation ran correctly and logged data to simOut.');
    return;
end

% Create a subfolder for plots if it doesn't exist
if ~exist(plotsSubFolder, 'dir')
    mkdir(plotsSubFolder);
    disp(['Created folder: ./' plotsSubFolder]);
end

% --- 2. Extract Data from simOut ---
disp('Extracting data from simOut object...');
try
    % Time vector
    timeVector_s = simOut.tout;
    timeVector_ms = timeVector_s * 1000; % Time in milliseconds for better plot resolution

    % Helper function to safely extract data from logsout
    % Assumes simOut.logsout is a Simulink.SimulationData.Dataset object
    function signalStruct = getLoggedSignal(logsoutObj, signalName)
        signalStruct = []; % Default to empty
        try
            if ~isempty(logsoutObj) && logsoutObj.exist(signalName)
                element = logsoutObj.getElement(signalName);
                signalStruct.Time = element.Values.Time;
                signalStruct.Data = squeeze(element.Values.Data); % Squeeze to remove singleton dimensions
                 % Ensure data is samples x dimensions (e.g., N x 3 for 3-phase)
                if size(signalStruct.Data, 1) < size(signalStruct.Data, 2) && size(signalStruct.Data, 1) ~= 1
                     signalStruct.Data = signalStruct.Data';
                end
            else
                disp(['Warning: Logged signal "' signalName '" not found in simOut.logsout.']);
            end
        catch ME_extract
            disp(['Warning: Error extracting signal "' signalName '": ' ME_extract.message]);
        end
    end

    % Extract key signals (names should match logging in CSSS_MainModel.xml or actual model)
    % These are examples; actual logged signal names might differ.
    % Assuming signals are logged via a "To Workspace" block named 'simout_CSSS' which creates logsout
    % or directly via signal logging. The XML used a ToFile block "SimulationDataLogger"
    % and a scope "KeyWaveformsScope". We'll assume signals are also available via simOut.logsout for this script.

    if isprop(simOut, 'logsout') && isa(simOut.logsout, 'Simulink.SimulationData.Dataset') && simOut.logsout.numElements > 0
        logs = simOut.logsout;

        % Input (Grid Side) Voltages and Currents
        V_in_abc_struct = getLoggedSignal(logs, 'V_in_abc'); % Placeholder name
        I_in_abc_struct = getLoggedSignal(logs, 'I_in_abc'); % Placeholder name

        % Output (Load Side) Voltages and Currents
        V_out_abc_struct = getLoggedSignal(logs, 'V_out_abc'); % Placeholder name
        I_out_abc_struct = getLoggedSignal(logs, 'I_out_abc'); % Placeholder name

        % Control System Signals (examples)
        V_out_RMS_actual_struct = getLoggedSignal(logs, 'V_out_RMS_actual'); % Placeholder
        Input_PF_actual_struct = getLoggedSignal(logs, 'Input_PF_actual');   % Placeholder
        FCL_Status_struct = getLoggedSignal(logs, 'FCL_Status');             % Placeholder

        % If specific signals were logged via ToFile block 'CSSS_simulation_data.mat',
        % you might need to load that MAT file separately if not included in simOut.
        % load('CSSS_simulation_data.mat'); % if logged_data struct exists there.

    else
        disp('Warning: simOut.logsout is empty or not found. Cannot extract detailed telemetry for plotting.');
        disp('Ensure signals are correctly logged in the Simulink model (e.g., using scopes with logging, ToWorkspace, or signal logging).');
        % Attempt to use timeVector_s for plots if available, but data will be missing.
    end

catch ME_data_extract
    disp(['Error extracting data from simOut: ' ME_data_extract.message]);
    disp('Check the names of logged signals and the structure of simOut.');
    % return; % Optional: stop script if critical data is missing
end

% --- 3. Generate Plots ---
disp('Generating plots...');

% Plot 1: Input Three-Phase Voltages (Grid Side)
if exist('V_in_abc_struct','var') && ~isempty(V_in_abc_struct) && size(V_in_abc_struct.Data,2) == 3
    figure('Name', 'Input Voltages (Grid Side)');
    plot(V_in_abc_struct.Time*1000, V_in_abc_struct.Data(:,1), 'DisplayName', 'Va_in'); hold on;
    plot(V_in_abc_struct.Time*1000, V_in_abc_struct.Data(:,2), 'DisplayName', 'Vb_in');
    plot(V_in_abc_struct.Time*1000, V_in_abc_struct.Data(:,3), 'DisplayName', 'Vc_in');
    title('Input Three-Phase Voltages (Grid Side)');
    xlabel('Time (ms)'); ylabel('Voltage (V)');
    legend show; grid on; xlim([max(0, timeVector_ms(end)-50) timeVector_ms(end)]); % Show last 50ms or full if shorter
    saveas(gcf, fullfile(plotsSubFolder, 'input_voltages_abc.png'));
    disp('Saved input_voltages_abc.png');
else
    disp('Skipping Input Voltages Plot: Data not available or incorrect dimensions.');
end

% Plot 2: Input Three-Phase Currents (Grid Side)
if exist('I_in_abc_struct','var') && ~isempty(I_in_abc_struct) && size(I_in_abc_struct.Data,2) == 3
    figure('Name', 'Input Currents (Grid Side)');
    plot(I_in_abc_struct.Time*1000, I_in_abc_struct.Data(:,1), 'DisplayName', 'Ia_in'); hold on;
    plot(I_in_abc_struct.Time*1000, I_in_abc_struct.Data(:,2), 'DisplayName', 'Ib_in');
    plot(I_in_abc_struct.Time*1000, I_in_abc_struct.Data(:,3), 'DisplayName', 'Ic_in');
    title('Input Three-Phase Currents (Grid Side)');
    xlabel('Time (ms)'); ylabel('Current (A)');
    legend show; grid on; xlim([max(0, timeVector_ms(end)-50) timeVector_ms(end)]);
    saveas(gcf, fullfile(plotsSubFolder, 'input_currents_abc.png'));
    disp('Saved input_currents_abc.png');
else
    disp('Skipping Input Currents Plot: Data not available or incorrect dimensions.');
end

% Plot 3: Output Three-Phase Voltages (Load Side)
if exist('V_out_abc_struct','var') && ~isempty(V_out_abc_struct) && size(V_out_abc_struct.Data,2) == 3
    figure('Name', 'Output Voltages (Load Side)');
    plot(V_out_abc_struct.Time*1000, V_out_abc_struct.Data(:,1), 'DisplayName', 'Va_out'); hold on;
    plot(V_out_abc_struct.Time*1000, V_out_abc_struct.Data(:,2), 'DisplayName', 'Vb_out');
    plot(V_out_abc_struct.Time*1000, V_out_abc_struct.Data(:,3), 'DisplayName', 'Vc_out');
    % Add setpoint line if available
    if exist('SimControlParams','var') && isfield(SimControlParams, 'V_out_RMS_Setpoint_V')
        yline(SimControlParams.V_out_RMS_Setpoint_V * sqrt(2), '--r', 'Vpeak Setpoint', 'LabelVerticalAlignment','bottom');
        yline(-SimControlParams.V_out_RMS_Setpoint_V * sqrt(2), '--r','-Vpeak Setpoint', 'LabelVerticalAlignment','top');
    end
    title('Output Three-Phase Voltages (Load Side)');
    xlabel('Time (ms)'); ylabel('Voltage (V)');
    legend show; grid on; xlim([max(0, timeVector_ms(end)-50) timeVector_ms(end)]);
    saveas(gcf, fullfile(plotsSubFolder, 'output_voltages_abc.png'));
    disp('Saved output_voltages_abc.png');
else
    disp('Skipping Output Voltages Plot: Data not available or incorrect dimensions.');
end

% Plot 4: Output Three-Phase Currents (Load Side)
if exist('I_out_abc_struct','var') && ~isempty(I_out_abc_struct) && size(I_out_abc_struct.Data,2) == 3
    figure('Name', 'Output Currents (Load Side)');
    plot(I_out_abc_struct.Time*1000, I_out_abc_struct.Data(:,1), 'DisplayName', 'Ia_out'); hold on;
    plot(I_out_abc_struct.Time*1000, I_out_abc_struct.Data(:,2), 'DisplayName', 'Ib_out');
    plot(I_out_abc_struct.Time*1000, I_out_abc_struct.Data(:,3), 'DisplayName', 'Ic_out');
    title('Output Three-Phase Currents (Load Side)');
    xlabel('Time (ms)'); ylabel('Current (A)');
    legend show; grid on; xlim([max(0, timeVector_ms(end)-50) timeVector_ms(end)]);
    saveas(gcf, fullfile(plotsSubFolder, 'output_currents_abc.png'));
    disp('Saved output_currents_abc.png');
else
    disp('Skipping Output Currents Plot: Data not available or incorrect dimensions.');
end

% Plot 5: Output Voltage RMS vs Setpoint (FR-03)
if exist('V_out_RMS_actual_struct','var') && ~isempty(V_out_RMS_actual_struct) && exist('SimControlParams','var')
    figure('Name', 'Output Voltage RMS vs Setpoint');
    plot(V_out_RMS_actual_struct.Time, V_out_RMS_actual_struct.Data, 'DisplayName', 'V_out RMS Actual'); hold on;
    yline(SimControlParams.V_out_RMS_Setpoint_V, '--r', 'V_out RMS Setpoint');
    yline(SimControlParams.V_out_RMS_Setpoint_V * 1.02, ':k', '+2% Limit');
    yline(SimControlParams.V_out_RMS_Setpoint_V * 0.98, ':k', '-2% Limit');
    title('Output Voltage RMS Regulation (FR-03)');
    xlabel('Time (s)'); ylabel('RMS Voltage (V_LN)');
    legend show; grid on;
    ylim([SimControlParams.V_out_RMS_Setpoint_V*0.95 SimControlParams.V_out_RMS_Setpoint_V*1.05]); % Zoom around setpoint
    saveas(gcf, fullfile(plotsSubFolder, 'output_voltage_rms.png'));
    disp('Saved output_voltage_rms.png');
else
    disp('Skipping Output Voltage RMS Plot: Data or Setpoint not available.');
end

% Plot 6: Input Power Factor (FR-04)
if exist('Input_PF_actual_struct','var') && ~isempty(Input_PF_actual_struct) && exist('SimControlParams','var')
    figure('Name', 'Input Power Factor');
    plot(Input_PF_actual_struct.Time, Input_PF_actual_struct.Data, 'DisplayName', 'Input PF Actual'); hold on;
    yline(SimControlParams.PF_Setpoint, '--r', 'PF Setpoint (0.99 Target)');
    title('Input Power Factor Correction (FR-04)');
    xlabel('Time (s)'); ylabel('Power Factor');
    legend show; grid on; ylim([0.9 1.05]); % Zoom around target PF
    saveas(gcf, fullfile(plotsSubFolder, 'input_power_factor.png'));
    disp('Saved input_power_factor.png');
else
    disp('Skipping Input Power Factor Plot: Data or Setpoint not available.');
end

% --- Advanced Plots (Placeholders - require more complex calculations) ---
% Plot 7: Harmonic Distortion (THD) - FR-06
% Requires FFT analysis of V_out_abc or I_in_abc.
% Example: if exist('V_out_abc_struct','var') && ~isempty(V_out_abc_struct)
%   % Select a steady-state portion of the signal
%   % Perform FFT: thd_V_out = thd(V_out_abc_struct.Data(:,1), Grid.Frequency_Hz, num_harmonics_to_consider);
%   disp('THD analysis placeholder: Requires Signal Processing Toolbox and steady-state data selection.');
% end

% Plot 8: Efficiency Analysis - FR-07
% Requires logging of input power and output power (or losses).
% P_in = V_in_abc' * I_in_abc (instantaneous, then average)
% P_out = V_out_abc' * I_out_abc (instantaneous, then average)
% Efficiency = P_out_avg / P_in_avg
% disp('Efficiency analysis placeholder: Requires calculation of average input/output power from waveforms.');


% --- 4. Generate Summary Report (Optional) ---
disp('--- CSSS Simulation Summary ---');
if exist('V_out_RMS_actual_struct','var') && ~isempty(V_out_RMS_actual_struct) && exist('SimControlParams','var')
    % Analyze last part of simulation for steady state
    settling_time_estimate = 0.1; % seconds, assume system settles after this time
    idx_steady_state = find(V_out_RMS_actual_struct.Time >= max(V_out_RMS_actual_struct.Time) - settling_time_estimate);
    if ~isempty(idx_steady_state)
        avg_V_out_RMS = mean(V_out_RMS_actual_struct.Data(idx_steady_state));
        regulation_error_percent = ((avg_V_out_RMS - SimControlParams.V_out_RMS_Setpoint_V) / SimControlParams.V_out_RMS_Setpoint_V) * 100;
        disp(['Average Output Voltage (L-N RMS, steady state): ' num2str(avg_V_out_RMS, '%.2f') ' V']);
        disp(['Setpoint Output Voltage (L-N RMS):              ' num2str(SimControlParams.V_out_RMS_Setpoint_V, '%.2f') ' V']);
        disp(['Voltage Regulation Error:                       ' num2str(regulation_error_percent, '%.3f') ' %']);
    else
        disp('Could not determine steady-state output voltage for summary.');
    end
end
if exist('Input_PF_actual_struct','var') && ~isempty(Input_PF_actual_struct)
    idx_steady_state_pf = find(Input_PF_actual_struct.Time >= max(Input_PF_actual_struct.Time) - settling_time_estimate);
     if ~isempty(idx_steady_state_pf)
        avg_Input_PF = mean(Input_PF_actual_struct.Data(idx_steady_state_pf));
        disp(['Average Input Power Factor (steady state):      ' num2str(avg_Input_PF, '%.4f')]);
     else
        disp('Could not determine steady-state input power factor for summary.');
     end
end
disp('-----------------------------');

disp('CSSS post-processing script finished.');
disp(['Plots saved in: ./' plotsSubFolder]);

% Close all figures (optional)
% close all;
