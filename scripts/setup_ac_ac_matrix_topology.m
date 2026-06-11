% setup_ac_ac_matrix_topology.m
% Initial parameterization and connectivity for direct AC-AC matrix converter topology.

disp('Setting up Simscape topology parameters...');

MatrixConverter.InputPhases = 3;
MatrixConverter.OutputPhases = 3;
MatrixConverter.SwitchingFreq = 20000;
MatrixConverter.SwitchRon = 0.012;
MatrixConverter.NumSwitches = MatrixConverter.InputPhases * MatrixConverter.OutputPhases;

disp('Matrix converter topology parameters drafted.');
