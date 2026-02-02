function [output, pOrig] = jm_run_virt_sim_vary_ext_current (simulationDir, simulationParamFile, simDur, intrinsicProbePeriod, vIRtInputCurrents, vIRtStdTimes, fmnInputCurrents, vIRtTauTime, varargin)
%% Iterate over parameter space to run multiple vIRt simulations (ExtCurrent)
% Usage: [output, pOrig] = jm_run_virt_sim_vary_ext_current (simulationDir, simulationParamFile, simDur, intrinsicProbePeriod, vIRtInputCurrents, vIRtStdTimes, fmnInputCurrents, vIRtTauTime, varargin)
% Explanation:
%       This function iterates over a defined parameter space to run multiple 
%       simulations of the vIRt model. It systematically varies external 
%       current inputs (Mean and Std for Rr/Rp, Mean for Fm).
%
%       The simulations are run in "Intrinsic" mode (Pb amplitude set to 0).
%
% Outputs:
%       output      - data struct containing model runs.
%                   specified as a structure with fields:
%                       .intrinsicDefault.vIRtInput(i).vIRtStdTime(j).FMNInput(k)
%                           .analysis: Analysis results.
%                           .params: Output parameters.
%                           .effectorVec: Effector position trace.
%       pOrig       - The original parameter structure loaded from file.
%                   specified as a structure
%
% Arguments:
%       simulationDir
%                   - Directory containing the version of simulation code to run.
%                   must be a string scalar or a character vector
%       simulationParamFile
%                   - Path to the base network parameter file (.mat).
%                   must be a string scalar or a character vector
%       simDur      - Simulation duration in milliseconds (ms).
%                   must be a positive numeric scalar
%       intrinsicProbePeriod 
%                   - Period (ms) with which to sample effector outputs.
%                   must be a positive numeric scalar
%       vIRtInputCurrents
%                   - Vector of Rr/Rp mean input currents to test (pA).
%                   must be a numeric vector
%       vIRtStdTimes
%                   - Vector of Rr/Rp stdTime input currents to test.
%                   must be a numeric vector
%       fmnInputCurrents
%                   - Vector of Fm mean input currents to test (pA).
%                   must be a numeric vector
%       vIRtTauTime 
%                   - Time constant for noise (ms).
%                   must be a positive numeric scalar
%       varargin    - 'RngAlgorithm': The random number generator algorithm to use.
%                       Options: 'twister', 'simdTwister', 'combRecursive',
%                                'multFibonacci', 'philox', 'threefry',
%                                'v4', 'v5uniform', 'v5normal'
%                   must be a string scalar or a character vector
%                   default == [] (set in virt_moore_params.m)
%                   - Any other parameter-value pair for virt_moore()
%
% Requires:
%       \Shared\Code\vIRt-Moore\virt_moore.m
%       \Shared\Code\vIRt-Moore\virt_updated_test_name_for_rep.m
%
% Used by:
%       \Shared\Code\vIRt-Moore\jm_automate_virt_sim_vary_ext_current.m

% File History:
% 2025-12-05 Created by Jeff Moore
% 2026-01-13 Refactored by Gemini to match vIRt coding standards.
% 2026-01-16 Fixed virt_moore output sequence
% 2026-01-16 Removed dummy message handler
% 2026-01-16 Refactored to use parfor for parallel execution.
% 2026-01-17 Force toRecompute == true for safety
% 2026-01-18 Added 'RngAlgorithm' optional argument
% 2026-01-18 Now updates testName in parameter structure using virt_updated_test_name_for_rep.m

%% Hard-coded parameters
validRngAlgorithms = {'twister', 'simdTwister', 'combRecursive', ...
                      'multFibonacci', 'philox', 'threefry', ...
                      'v4', 'v5uniform', 'v5normal'};

%% Default values for optional arguments
rngAlgorithmDefault = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options for virt_moore

% Add required inputs to the Input Parser
addRequired(iP, 'simulationDir', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addRequired(iP, 'simulationParamFile', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addRequired(iP, 'simDur', ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addRequired(iP, 'intrinsicProbePeriod', ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addRequired(iP, 'vIRtInputCurrents', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addRequired(iP, 'vIRtStdTimes', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addRequired(iP, 'fmnInputCurrents', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addRequired(iP, 'vIRtTauTime', ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'RngAlgorithm', rngAlgorithmDefault, ...
    @(x) isempty(x) || any(validatestring(x, validRngAlgorithms)));

% Parse the inputs
parse(iP, simulationDir, simulationParamFile, simDur, intrinsicProbePeriod, ...
      vIRtInputCurrents, vIRtStdTimes, fmnInputCurrents, vIRtTauTime, varargin{:});

rngAlgorithm = iP.Results.RngAlgorithm;

% Extract extra arguments to pass down
extraArgs = iP.Unmatched;

%% Setup Environment
% Load parameters
if ~exist(simulationParamFile, 'file')
    error('Parameter file not found: %s', simulationParamFile);
end
loadedData = load(simulationParamFile);
if isfield(loadedData, 'P')
    p = loadedData.P;
else
    error('Parameter file must contain a variable named "P".');
end

% Save original parameters
pOrig = p;

% Update simulation duration
p.simDur = simDur;

% Ensure we return to original directory after the simulations are complete
currentDir = pwd;
cleanupObj = onCleanup(@() cd(currentDir));

% Navigate to simulation directory
cd(simulationDir);

% Convert extraArgs struct to name-value pairs for function calls
extraArgsCell = namedargs2cell(extraArgs);

%% Simulation Loop (Parallelized)
n1 = length(vIRtInputCurrents);
n2 = length(vIRtStdTimes);
n3 = length(fmnInputCurrents);
totalSims = n1 * n2 * n3;

fprintf('Starting Parameter Sweep (%d simulations) using parfor...\n', totalSims);

% Pre-allocate a cell array to store results linearly for parfor
parResults = cell(totalSims, 1);

% Flatten the loop for parallel execution
parfor idx = 1:totalSims
    fprintf('  > Simulation %d of %d...\n', idx, totalSims);

    % Convert linear index back to subscripts (i, j, k)
    [i, j, k] = ind2sub([n1, n2, n3], idx);
    
    % Create a temporary parameter structure for this worker
    pTemp = p; 
    
    % Set Pb input parameters (Intrinsic Mode)
    pTemp.SquareInput.Pb.amplitude = 0;
    pTemp.SquareInput.Pb.period1 = intrinsicProbePeriod;
    pTemp.SquareInput.Pb.perVarPerc1 = 0;
    pTemp.SquareInput.Pb.perVar1 = 0;
    pTemp.SquareInput.Pb.period2 = intrinsicProbePeriod;
    pTemp.SquareInput.Pb.perVarPerc2 = 0;
    pTemp.SquareInput.Pb.perVar2 = 0;
    pTemp.SquareInput.Pb.toRecompute = true;
    
    % Set input current params
    pTemp.ExtCurrent.Rr.mean = vIRtInputCurrents(i);
    pTemp.ExtCurrent.Rp.mean = vIRtInputCurrents(i);
    pTemp.ExtCurrent.Rr.stdTime = vIRtStdTimes(j);
    pTemp.ExtCurrent.Rp.stdTime = vIRtStdTimes(j);
    pTemp.ExtCurrent.Fm.mean = fmnInputCurrents(k);
    pTemp.ExtCurrent.Rr.tauTime = vIRtTauTime; % Set tauTime for Rr
    pTemp.ExtCurrent.Rp.tauTime = vIRtTauTime; % Set tauTime for Rp
    
    % --- Update Test Name for this parameter combination ---
    % 1. Add vIRt Input Mean
    strRrRpMean = replace(sprintf('IExtRrRp-%g', vIRtInputCurrents(i)), '.', 'p');
    pTemp.testName = virt_updated_test_name_for_rep(pTemp.testName, 'IExtRrRp', strRrRpMean);
    
    % 2. Add vIRt Std Time
    strRrRpStd = replace(sprintf('RrRpStdTime-%g', vIRtStdTimes(j)), '.', 'p');
    pTemp.testName = virt_updated_test_name_for_rep(pTemp.testName, 'StdTime', strRrRpStd);

    % 3. Add Fm Input Mean
    strFmMean = replace(sprintf('IExtFm-%g', fmnInputCurrents(k)), '.', 'p');
    pTemp.testName = virt_updated_test_name_for_rep(pTemp.testName, 'IExtFm', strFmMean);
    
    fprintf('  > Running %s\n', pTemp.testName);

    % Run simulation
    % Use flags to disable plotting for efficiency during sweep
    [~, pOut1, ~, analysis1, ~, plotData1] = virt_moore('Params', pTemp, ...
        'RandomizeSeed', false, ...
        'RngAlgorithm', rngAlgorithm, ...
        'SaveOutput', false, ...
        'PlotSampleTraces', false, ...
        'PlotRaster', false, ...
        'PlotJitter', false, ...
        'PlotScatter', false, ...
        'PlotPhaseResponse', false, ...
        'CopyScript', false, ...
        extraArgsCell{:});
    
    % Store Results in a temporary structure
    res = struct();
    res.analysis = analysis1;
    res.params = pOut1;
    res.effectorVec = plotData1.effectorVec;
    
    parResults{idx} = res;
end

% Reconstruct the nested output structure from parallel results
fprintf('Reconstructing output structure...\n');
for idx = 1:totalSims
    [i, j, k] = ind2sub([n1, n2, n3], idx);
    res = parResults{idx};
    
    output.intrinsicDefault.vIRtInput(i).vIRtStdTime(j).FMNInput(k).analysis = res.analysis;
    output.intrinsicDefault.vIRtInput(i).vIRtStdTime(j).FMNInput(k).params = res.params;
    output.intrinsicDefault.vIRtInput(i).vIRtStdTime(j).FMNInput(k).effectorVec = res.effectorVec;
end

fprintf('Sweep complete.\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%