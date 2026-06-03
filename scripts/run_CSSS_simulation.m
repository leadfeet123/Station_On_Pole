% run_CSSS_simulation.m
%
% Description:
% This script runs the Compact Solid-State Substation (CSSS) Simulink simulation.
% It loads parameters, configures, and executes the main model 'CSSS_MainModel.slx'.
% Simulation results are saved to 'CSSS_simulation_output.mat'.
%
% Version: 1.0
% Date:    YYYY-MM-DD % To be filled

disp('Starting CSSS simulation run script...');

% --- Configuration ---
mainModelName = 'CSSS_MainModel'; % Name of the main Simulink model file (without .slx)
paramScriptName = 'initialize_CSSS_parameters';
paramDataFile = 'CSSS_simulation_parameters.mat';
outputDataFile = 'CSSS_simulation_output.mat';

% --- 1. Initialize and Load Parameters ---
if ~exist(paramDataFile, 'file')
    disp(['Parameter file ' paramDataFile ' not found. Running ' paramScriptName '.m ...']);
    try
        run(paramScriptName); % This script should save paramDataFile
        if ~exist(paramDataFile, 'file')
            error([paramScriptName '.m did not create ' paramDataFile '. Cannot proceed.']);
        end
        disp([paramDataFile ' created by ' paramScriptName '.m.']);
    catch ME_init
        disp(['Error running ' paramScriptName '.m: ' ME_init.message]);
        disp('Please ensure initialize_CSSS_parameters.m can run successfully.');
        return; % Stop execution
    end
end

disp(['Loading parameters from ' paramDataFile ' ...']);
load(paramDataFile); % Loads AllParams, SimControl, Grid, Load, Converter, etc.
disp('Parameters loaded successfully.');

% --- 2. Check for Model File ---
modelFileWithExt = [mainModelName '.slx'];
if ~exist(modelFileWithExt, 'file')
    modelFileWithExt_mdl = [mainModelName '.mdl']; % Check for older format
    if ~exist(modelFileWithExt_mdl, 'file')
        disp(['Error: Main Simulink model file (' modelFileWithExt ' or ' modelFileWithExt_mdl ') not found.']);
        disp('Please ensure the model exists in the MATLAB path or current directory.');
        return; % Stop execution
    else
        modelFileWithExt = modelFileWithExt_mdl;
    end
end
disp(['Found main model: ' modelFileWithExt]);

% --- 3. Load and Run Simulation ---
disp(['Loading Simulink model: ' mainModelName]);
try
    % Load the model into memory if not already loaded
    if ~bdIsLoaded(mainModelName)
        load_system(mainModelName);
        disp(['Model ' mainModelName ' loaded.']);
    else
        disp(['Model ' mainModelName ' is already loaded.']);
    end

    % Create a SimulationInput object
    simIn = Simulink.SimulationInput(mainModelName);

    % Make parameter structures available to the model workspace
    % The XML descriptions assume these structs (SimControl, Grid, etc.) are available.
    simIn = simIn.setVariable('SimControl', SimControl);
    simIn = simIn.setVariable('Grid', Grid);
    simIn = simIn.setVariable('Load', Load);
    simIn = simIn.setVariable('Converter', Converter);
    simIn = simIn.setVariable('SimMOSFET', SimMOSFET);
    simIn = simIn.setVariable('SimAutotransformer', SimAutotransformer);
    simIn = simIn.setVariable('SimControlParams', SimControlParams);
    if exist('AllParams', 'var') % If AllParams struct exists, pass it too
       simIn = simIn.setVariable('AllParams', AllParams);
    end
    disp('Parameter structures assigned to model workspace via SimulationInput object.');

    % Configure simulation settings from SimControl struct
    simIn = simIn.setModelParameter('StartTime', num2str(SimControl.StartTime), ...
                                    'StopTime', num2str(SimControl.StopTime), ...
                                    'Solver', SimControl.SolverType, ...
                                    'MaxStep', num2str(SimControl.MaxStep), ...
                                    'RelTol', num2str(SimControl.RelTol), ...
                                    'AbsTol', SimControl.AbsTol); % AbsTol can be 'auto' or numeric string
    disp(['Simulation settings configured: Solver=' SimControl.SolverType ...
          ', MaxStep=' num2str(SimControl.MaxStep) ...
          ', StopTime=' num2str(SimControl.StopTime) 's.']);

    % Set Simscape local solver options if specified
    if isfield(SimControl, 'SimscapeLocalSolverSampleTime') && ~isempty(SimControl.SimscapeLocalSolverSampleTime)
        simIn = simIn.setModelParameter('SimscapeUseLocalSolver','on', ...
                                        'SimscapeLocalSolverChoice','NE_BACKWARD_EULER_ADVANCER', ... % Example solver
                                        'SimscapeLocalSolverSampleTime', num2str(SimControl.SimscapeLocalSolverSampleTime));
        disp(['Simscape local solver enabled with sample time: ' num2str(SimControl.SimscapeLocalSolverSampleTime) 's.']);
    else
        simIn = simIn.setModelParameter('SimscapeUseLocalSolver','off');
    end

    % Enable fast restart for iterative simulations if desired (usually off for single runs)
    % simIn = simIn.setModelParameter('FastRestart','off');

    disp(['Running simulation for ' mainModelName '...']);

    % Run the simulation
    simOut = [];
    simException = [];
    try
        simOut = sim(simIn);
        disp('Simulation finished successfully.');
    catch ME_sim
        simException = ME_sim;
        disp(['Error during simulation run: ' ME_sim.message]);
        fprintf('Error details:\n');
        for k=1:length(ME_sim.stack)
            fprintf('  File: %s, Name: %s, Line: %d\n', ME_sim.stack(k).file, ME_sim.stack(k).name, ME_sim.stack(k).line);
        end
    end

    % --- 4. Save Simulation Output ---
    if ~isempty(simOut)
        disp(['Saving simulation output to ' outputDataFile '...']);
        % Save the Simulink.SimulationOutput object and the parameters used
        save(outputDataFile, 'simOut', 'AllParams', 'SimControl', 'Grid', 'Load', 'Converter', 'SimMOSFET', 'SimAutotransformer', 'SimControlParams');
        disp('Simulation output saved.');

        % Check if any ToFile blocks were used (as described in CSSS_MainModel.xml)
        if exist('CSSS_simulation_data.mat', 'file')
            disp('Data also logged by "ToFile" blocks in the model to CSSS_simulation_data.mat');
        end
    else
        disp('No simulation output object (simOut) generated. This might be due to an error or no logged signals.');
    end

    if ~isempty(simException)
        disp('Simulation completed with errors. Check console output and simException variable in workspace.');
        % Consider rethrowing the error if this script is part of a larger automated process
        % rethrow(simException);
    end

    % --- 5. Close Model (Optional) ---
    % To prevent accidental changes and free up memory.
    % Set second argument to 0 to close without saving any model changes.
    % close_system(mainModelName, 0);
    % disp(['Model ' mainModelName ' closed.']);

catch ME_main
    disp(['An error occurred in run_CSSS_simulation.m: ' ME_main.message]);
    fprintf('Error details:\n');
    for k=1:length(ME_main.stack)
        fprintf('  File: %s, Name: %s, Line: %d\n', ME_main.stack(k).file, ME_main.stack(k).name, ME_main.stack(k).line);
    end
end

disp('run_CSSS_simulation.m script finished.');

% --- Notes ---
% To run this script:
% 1. Ensure MATLAB is in the project's root directory or that the directory
%    is added to the MATLAB path.
% 2. Ensure Simulink and Simscape Electrical toolboxes are installed and licensed.
% 3. Ensure 'initialize_CSSS_parameters.m' and 'CSSS_MainModel.slx' exist.
% 4. Execute `run_CSSS_simulation` in the MATLAB command window.
