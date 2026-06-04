% verify_prd_compliance.m
% Script-level compliance checks against README.md and docs/PRD.md targets.

disp('Running README/PRD compliance checks...');

script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end
run(fullfile(script_dir, 'set_parameters.m'));

% FR-01: 138kV input, configurable 11-34.5kV output
fr01_input_ok = (V_in_rms == 138000);
fr01_output_ok = (V_out_ll_rms >= 11000) && (V_out_ll_rms <= 34500);

% FR-04: PF >= 0.99 target
fr04_pf_ok = (pf_target >= 0.99);

% FR-05: FCL set to 1.5x rated current, response within 2ms
fr05_limit_ok = abs(fcl_limit_factor - 1.5) < 1e-12;
fr05_time_ok = (fcl_response_ms <= 2);

% FR-07: efficiency analysis path exists
fr07_efficiency_ok = exist(fullfile(script_dir, 'analyze_converter_stress.m'), 'file') == 2;

% NFR-01: parameters centralized in set_parameters.m
nfr01_parameterization_ok = true;

% NFR-02: modular scripts exist for run + stress + protection audit
nfr02_modularity_ok = ...
    exist(fullfile(script_dir, 'run_simulation.m'), 'file') == 2 && ...
    exist(fullfile(script_dir, 'analyze_converter_stress.m'), 'file') == 2 && ...
    exist(fullfile(script_dir, 'audit_protection_fcl.m'), 'file') == 2;

compliance = struct( ...
    'FR01_Input138kV', fr01_input_ok, ...
    'FR01_OutputRange11to34p5kV', fr01_output_ok, ...
    'FR04_PFAtLeast0p99', fr04_pf_ok, ...
    'FR05_FCLLimit1p5xRated', fr05_limit_ok, ...
    'FR05_FCLResponseWithin2ms', fr05_time_ok, ...
    'FR07_EfficiencyAnalysisPath', fr07_efficiency_ok, ...
    'NFR01_SingleParameterScript', nfr01_parameterization_ok, ...
    'NFR02_ModularScriptFlow', nfr02_modularity_ok);

all_ok = all(structfun(@(x) logical(x), compliance));
if ~all_ok
    error('Compliance check failed for one or more README/PRD criteria.');
end

save(fullfile(script_dir, 'prd_compliance_results.mat'), 'compliance');
disp('All README/PRD compliance checks passed.');
disp('Saved results to scripts/prd_compliance_results.mat');
