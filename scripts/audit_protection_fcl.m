% audit_protection_fcl.m
% Basic overcurrent/FCL parameter audit for CSSS scripts.

disp('Running protection and FCL audit...');

script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end
run(fullfile(script_dir, 'set_parameters.m'));

if ~(fcl_trip_threshold_rms > I_rated_rms)
    error('FCL trip threshold must be above rated current.');
end

if ~(fcl_release_threshold_rms < fcl_trip_threshold_rms)
    error('FCL release threshold must be below trip threshold.');
end

if ~(fcl_detection_samples >= 1)
    error('FCL detection samples must be at least 1.');
end

fprintf('\n--- Protection & FCL Audit Summary ---\n');
fprintf('Rated current (RMS):          %.3f A\n', I_rated_rms);
fprintf('Trip threshold (RMS):         %.3f A\n', fcl_trip_threshold_rms);
fprintf('Release threshold (RMS):      %.3f A\n', fcl_release_threshold_rms);
fprintf('FCL response target:          %.3f ms\n', fcl_response_ms);
fprintf('Detection window (samples):   %d\n', fcl_detection_samples);
fprintf('--------------------------------------\n\n');

FCLAudit = struct( ...
    'rated_current_rms_a', I_rated_rms, ...
    'trip_threshold_rms_a', fcl_trip_threshold_rms, ...
    'release_threshold_rms_a', fcl_release_threshold_rms, ...
    'response_time_ms', fcl_response_ms, ...
    'detection_samples', fcl_detection_samples);

save(fullfile(script_dir, 'fcl_audit_results.mat'), 'FCLAudit');
disp('Saved audit results to scripts/fcl_audit_results.mat');
