function jm_automate_virt_sim_vary_ext_current (varargin)
%% Automate vIRt simulations with varying external currents
% Usage: jm_automate_virt_sim_vary_ext_current (varargin)
% Explanation:
%       This function iterates over a list of parameter files and runs
%       simulations varying external current parameters (Mean Rr/Rp, 
%       Std Rr/Rp, Mean Fm).
%
%       It calls 'jm_run_virt_sim_vary_ext_current' to perform the actual
%       simulation sweeps for each parameter file.
%
%       FORCED PARAMETERS:
%       The following parameters are enforced on all simulations to ensure 
%       consistency, regardless of the input parameter file settings:
%           - SimDur = 30000 ms
%           - RelativeAnalysisStart = 0.5
%           - RelativeTransitionTime = 1.0
%
% Example(s):
%       % Run with default settings (uses virt_moore_params defaults)
%       jm_automate_virt_sim_vary_ext_current;
%
%       % Run with custom current ranges
%       jm_automate_virt_sim_vary_ext_current('VIRtInputCurrents', 100:50:500);
%
% Outputs:
%       This creates the output
%           sprintf('simData_ExtCurrent_%s.mat', fName) 
%           for each parameter file name fName
%       required for jm_postprocess_virt_sim_ExtCurrent_analysis.m
%           and jm_automate_rerun_virt_sim_from_params_intrinsic_vs_driven_freq.m
%           and jm_postprocess_virt_sim_ExtCurrent_analysis.m
%       which produces Figure 6b - 6d and Figure 7b - 7e in Lu et al 2026
%
% Arguments:
%       varargin    - 'SimulationDir': Directory containing the simulation code.
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'ParamFileList': List of parameter files to process.
%                   must be a cell array of character vectors or string array
%                   default == {} (Uses virt_moore_params default)
%                   - 'SimDur': Simulation duration in ms.
%                   must be a positive numeric scalar
%                   default == 30000
%                   - 'IntrinsicProbePeriod': Probe period for intrinsic amp (ms).
%                   must be a positive numeric scalar
%                   default == 200
%                   - 'VIRtInputCurrents': Vector of Rr/Rp mean input currents (pA).
%                   must be a numeric vector
%                   default == 100:30:520
%                   - 'VIRtStdTimes': Vector of Rr/Rp stdTime input currents.
%                   must be a numeric vector
%                   default == 0:5:60
%                   - 'FMNInputCurrents': Vector of Fm mean input currents (pA).
%                   must be a numeric vector
%                   default == 270:20:430
%                   - 'VIRtTauTime': Time constant for noise color (ms).
%                   must be a positive numeric scalar
%                   default == 5
%                   - 'OutDir': Directory to save output .mat files.
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'RngAlgorithm': The random number generator algorithm to use.
%                       Options: 'twister', 'simdTwister', 'combRecursive',
%                                'multFibonacci', 'philox', 'threefry',
%                                'v4', 'v5uniform', 'v5normal'
%                   must be a string scalar or a character vector
%                   default == [] (set in virt_moore_params.m)
%
% Requires:
%       \Shared\Code\Adams_Functions\archive_dependent_scripts.m
%       \Shared\Code\vIRt-Moore\jm_run_virt_sim_vary_ext_current.m
%       \Shared\Code\vIRt-Moore\virt_moore_params.m
%
% Used by:
%       (User automated scripts)

% File History:
% 2025-12-05 Created by Jeff Moore
% 2026-01-09 Updated parameters and file lists.
% 2026-01-13 Refactored by Gemini to match vIRt coding standards.
% 2026-01-15 Updated defaults to use virt_moore_params and enforce specific
%            simulation timing parameters (simDur=30s, relStart=0.5, relTrans=1).
% 2026-01-16 Now saves default parameters used in a file
% 2026-01-16 Now archives dependent scripts
% 2026-01-18 Added 'RngAlgorithm' optional argument

%% Hard-Coded Information For This Script Only
% Forced Parameters (Enforced on all loaded parameter files)
forcedSimDur = 30000;                % ms
forcedRelativeAnalysisStart = 0.5;   % fraction
forcedRelativeTransitionTime = 1.0;  % fraction

validRngAlgorithms = {'twister', 'simdTwister', 'combRecursive', ...
                      'multFibonacci', 'philox', 'threefry', ...
                      'v4', 'v5uniform', 'v5normal'};

%% Default values for optional arguments
simulationDirDefault = pwd;
paramFileListDefault = {};
simDurDefault = 30000;                  % simulation duration
intrinsicProbePeriodDefault = 200;      % probe period
vIRtInputCurrentsDefault = 100:30:520;  % pA
vIRtStdTimesDefault = 0:5:60;           % pA
fmnInputCurrentsDefault = 270:20:430;   % pA
vIRtTauTimeDefault = 5;                 % ms
outDirDefault = pwd;
rngAlgorithmDefault = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = false;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'SimulationDir', simulationDirDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'ParamFileList', paramFileListDefault, ...
    @(x) iscellstr(x) || isstring(x));
addParameter(iP, 'SimDur', simDurDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'IntrinsicProbePeriod', intrinsicProbePeriodDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'VIRtInputCurrents', vIRtInputCurrentsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'VIRtStdTimes', vIRtStdTimesDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'FMNInputCurrents', fmnInputCurrentsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'VIRtTauTime', vIRtTauTimeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'OutDir', outDirDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'RngAlgorithm', rngAlgorithmDefault, ...
    @(x) isempty(x) || any(validatestring(x, validRngAlgorithms)));

% Parse the inputs
parse(iP, varargin{:});
simulationDir = iP.Results.SimulationDir;
paramFileList = iP.Results.ParamFileList;
simDur = iP.Results.SimDur;
intrinsicProbePeriod = iP.Results.IntrinsicProbePeriod;
vIRtInputCurrents = iP.Results.VIRtInputCurrents;
vIRtStdTimes = iP.Results.VIRtStdTimes;
fmnInputCurrents = iP.Results.FMNInputCurrents;
vIRtTauTime = iP.Results.VIRtTauTime;
outDir = iP.Results.OutDir;
rngAlgorithm = iP.Results.RngAlgorithm;

% Ensure output directory exists
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

% Copy this script and its dependencies to the output directory
archive_dependent_scripts(mfilename, 'OutFolder', outDir);

%% Simulation Loop
fprintf('Starting automated vIRt simulation (Vary External Current)...\n');
fprintf('Simulation Directory: %s\n', simulationDir);

% Handle empty paramFileList by creating a placeholder for default generation
if isempty(paramFileList)
    paramFileList = {'virt_moore_params_default'};
end

nFiles = size(paramFileList, 1);

for iFile = 1:nFiles
    % Get current file information
    currentFile = paramFileList{iFile};
    
    fprintf('--- Processing File %d of %d ---\n', iFile, nFiles);
    
    % --- 1. Load or Generate Parameters ---
    if strcmp(currentFile, 'virt_moore_params_default')
        fName = 'virt_moore_params_default';
        fprintf('   > Generating default parameters using virt_moore_params...\n');
        P = virt_moore_params();
    else
        [~, fName] = fileparts(currentFile);
        fprintf('   > Loading parameters from: %s\n', fName);
        loadedData = load(currentFile);
        if isfield(loadedData, 'P')
            P = loadedData.P;
        else
            error('Parameter file must contain a variable named "P".');
        end
    end
    
    % --- 2. Enforce Specific Parameters & Notify ---
    % Check and modify SimDur
    if ~isfield(P, 'simDur') || P.simDur ~= forcedSimDur
        fprintf('     [NOTICE] Changing P.simDur from %g to %g ms.\n', ...
            P.simDur, forcedSimDur);
        P.simDur = forcedSimDur;
    end
    
    % Check and modify RelativeAnalysisStart
    if ~isfield(P.Analysis, 'relativeAnalysisStart') || ...
            P.Analysis.relativeAnalysisStart ~= forcedRelativeAnalysisStart
        fprintf('     [NOTICE] Changing P.Analysis.relativeAnalysisStart from %g to %g.\n', ...
            P.Analysis.relativeAnalysisStart, forcedRelativeAnalysisStart);
        P.Analysis.relativeAnalysisStart = forcedRelativeAnalysisStart;
    end
    
    % Check and modify RelativeTransitionTime
    if ~isfield(P.SquareInput.Pb, 'relativeTransitionTime') || ...
            P.SquareInput.Pb.relativeTransitionTime ~= forcedRelativeTransitionTime
        fprintf('     [NOTICE] Changing P.SquareInput.Pb.relativeTransitionTime from %g to %g.\n', ...
             P.SquareInput.Pb.relativeTransitionTime, forcedRelativeTransitionTime);
        P.SquareInput.Pb.relativeTransitionTime = forcedRelativeTransitionTime;
    end
    
    % Update local simDur variable to match P (to ensure consistency with runner call)
    simDur = P.simDur;
    
    % --- 3. Save Parameter File ---
    % The runner function expects a file path. We save the modified P to a temp file.
    paramFilePath = fullfile(outDir, sprintf('%s_Params.mat', fName));
    save(paramFilePath, 'P');
    
    % Begin computing
    % Calls the runner function to handle the sweeps
    [output, pOrig] = jm_run_virt_sim_vary_ext_current(simulationDir, ...
                                                       paramFilePath, ...
                                                       simDur, ...
                                                       intrinsicProbePeriod, ...
                                                       vIRtInputCurrents, ...
                                                       vIRtStdTimes, ...
                                                       fmnInputCurrents, ...
                                                       vIRtTauTime, ...
                                                       'RngAlgorithm', rngAlgorithm);
    
    % Construct filename for saving final results
    saveFileName = fullfile(outDir, ...
        sprintf('simData_ExtCurrent_%s.mat', fName));
    
    fprintf('     Saving results to: %s\n', saveFileName);
    
    % Save outputs (using -v7.3 for large data support)
    % Note: Variable names in save match original script conventions for consistency
    % We save 'fPath' as empty or the output dir since we generated/modified it
    fPath = outDir; 
    save(saveFileName, 'fPath', 'fName', 'vIRtInputCurrents', ...
         'vIRtStdTimes', 'fmnInputCurrents', 'output', 'pOrig', '-v7.3');
     
    fprintf('File %d complete.\n\n', iFile);
end

fprintf('All automated external current simulations complete.\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%