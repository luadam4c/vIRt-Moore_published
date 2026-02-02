function [output, pOrig] = jm_run_virt_sim_replicates (simulationDir, simulationParamFile, simDur, intrinsicProbePeriod, testFrequencies, testPerVar, numRandomSims, varargin)
%% Iterate over parameter space to run multiple vIRt simulations
% Usage: [output, pOrig] = jm_run_virt_sim_replicates (simulationDir, simulationParamFile, simDur, intrinsicProbePeriod, testFrequencies, testPerVar, numRandomSims, varargin)
% Explanation:
%       This function iterates over a defined parameter space to run multiple 
%       simulations of the vIRt model. It systematically varies input 
%       parameters (frequency and period variability) and runs replicates 
%       with randomized seeds.
%
%       It performs the following steps:
%       1. Runs a default "Intrinsic" simulation (no drive).
%       2. Runs N randomized "Intrinsic" replicates.
%       3. Performs a grid sweep over test frequencies and period variability.
%       4. For each grid point, runs a default "Driven" simulation.
%       5. For each grid point, runs N randomized "Driven" replicates.
%
%       REFACTOR NOTE: This function now utilizes 'virt_moore_multiple_reps'
%       to handle the randomized replicates, ensuring consistent seed management
%       and execution logic.
%
% Example(s):
%       [out, p] = jm_run_virt_sim_replicates('C:\Code\vIRt-Moore', 'default_params.mat', ...
%                   5000, 200, [2, 4, 8], [0, 0.1, 0.2], 5);
%
% Outputs:
%       output      - data struct containing model runs.
%                   specified as a structure with fields:
%                       .intrinsicDefault: Results for default intrinsic run.
%                       .intrinsicRep:     Struct array of randomized intrinsic runs.
%                       .drivenDefault:    Nested struct of frequency/variance sweeps.
%                       .drivenRep:        Nested struct of randomized driven runs.
%       pOrig       - The original parameter structure loaded from file.
%                   specified as a structure
%
% Arguments:
%       simulationDir
%                   - Directory containing the version of simulation code to run.
%                   must be a string scalar or a character vector
%       simulationParamFile
%                   - Path to the base network parameter file (.mat) from a GUI run.
%                   must be a string scalar or a character vector
%       simDur      - Simulation duration in milliseconds (ms).
%                   must be a positive numeric scalar
%       intrinsicProbePeriod 
%                   - Period (ms) with which to sample effector outputs during
%                     intrinsic runs.
%                   must be a positive numeric scalar
%       testFrequencies
%                   - Vector of Pb drive frequencies (Hz) to test.
%                   must be a numeric vector
%       testPerVar  - Vector of period variability values (0-1) to test.
%                     Expressed as a fraction of the period.
%                   must be a numeric vector
%       numRandomSims
%                   - Number of simulation re-runs with randomized seeds.
%                   must be a non-negative integer scalar
%       varargin    - 'RngAlgorithm': The random number generator algorithm to use.
%                       Options: 'twister', 'simdTwister', 'combRecursive',
%                                'multFibonacci', 'philox', 'threefry',
%                                'v4', 'v5uniform', 'v5normal'
%                   must be a string scalar or a character vector
%                   default == [] (set in virt_moore_params.m)
%                   - Any other parameter-value pair for virt_moore()
%
% Requires:
%       \Shared\Code\vIRt-Moore\virt_moore_multiple_reps.m
%       \Shared\Code\vIRt-Moore\virt_moore.m
%       \Shared\Code\vIRt-Moore\virt_updated_test_name_for_rep.m
%
% Used by:
%       \Shared\Code\vIRt-Moore\jm_automate_virt_sim_vary_freq_perVar.m
%       \Shared\Code\vIRt-Moore\jm_automate_virt_sim_vary_freq_reps.m
%       \Shared\Code\vIRt-Moore\jm_automate_rerun_virt_sim_from_params_intrinsic_vs_driven_freq.m

% File History:
% 2025-12-01 Created by Jeff Moore
% 2026-01-13 Adapted/Reorganized by Gemini to match vIRt code style.
% 2025-01-13 Refactored to use virt_moore_multiple_reps.m
% 2026-01-16 Fixed virt_moore output sequence
% 2026-01-16 Removed dummy message handler
% 2026-01-16 Refactored to use parfor for Driven Frequency sweep.
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
addRequired(iP, 'testFrequencies', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addRequired(iP, 'testPerVar', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addRequired(iP, 'numRandomSims', ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'nonnegative'}));
addParameter(iP, 'RngAlgorithm', rngAlgorithmDefault, ...
    @(x) isempty(x) || any(validatestring(x, validRngAlgorithms)));

% Parse the inputs
parse(iP, simulationDir, simulationParamFile, simDur, intrinsicProbePeriod, ...
      testFrequencies, testPerVar, numRandomSims, varargin{:});

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

% Common plotting flags to disable graphics during batch runs
batchPlotFlags = {'PlotSampleTraces', false, 'PlotRaster', false, ...
                  'PlotJitter', false, 'PlotScatter', false, ...
                  'PlotPhaseResponse', false, 'ShowFigure', false, ...
                  'SaveOutput', false, 'CopyScript', false};

% Convert extraArgs struct to name-value pairs for function calls
extraArgsCell = namedargs2cell(extraArgs);

%% 1. Intrinsic Simulation
% Run once to calculate intrinsic frequency and amplitude
p.SquareInput.Pb.amplitude = 0;
p.SquareInput.Pb.period1 = intrinsicProbePeriod;
p.SquareInput.Pb.perVarPerc1 = 0;
p.SquareInput.Pb.perVar1 = 0;
p.SquareInput.Pb.period2 = intrinsicProbePeriod;
p.SquareInput.Pb.perVarPerc2 = 0;
p.SquareInput.Pb.perVar2 = 0;

fprintf('Running Intrinsic Simulation (Default)...\n');
% We keep the direct call for the single default run
[~, pOut1, ~, analysis1, ~, plotData1] = virt_moore('Params', p, ...
    'RandomizeSeed', false, ...
    'RngAlgorithm', rngAlgorithm, ...
    batchPlotFlags{:}, ...
    extraArgsCell{:});

output.intrinsicDefault.analysis = analysis1;
output.intrinsicDefault.params = pOut1;
output.intrinsicDefault.data = plotData1;

% --- Refactored: Use virt_moore_multiple_reps for replicates ---
if numRandomSims > 0
    fprintf('Running %d Intrinsic Replicates...\n', numRandomSims);
    
    [~, ~, pAll, ~, analysisAll, ~, plotDataAll] = virt_moore_multiple_reps(...
        'Params', p, ...
        'NRepetitions', numRandomSims, ...
        'RepetitionMode', 'MonteCarlo', ...
        'RandomizeSeed', true, ...
        'UseParpool', true, ...
        'RngAlgorithm', rngAlgorithm, ...
        batchPlotFlags{:}, ...
        extraArgsCell{:});
    
    % Map the results back to the legacy output structure
    % Now we map plotDataAll to .data
    for i = 1:numRandomSims
        output.intrinsicRep(i).analysis = analysisAll(i);
        output.intrinsicRep(i).params = pAll(i);
        output.intrinsicRep(i).data = plotDataAll(i); 
    end
end

%% 2. Driven Simulations (Frequency Sweep - Parallelized)
numFreq = length(testFrequencies);
numVar = length(testPerVar);
totalSteps = numFreq * numVar;

fprintf('Starting Frequency Sweep (%d frequencies, %d variances) using parfor...\n', ...
    numFreq, numVar);

% Pre-allocate cell array to store results from parfor
parResults = cell(totalSteps, 1);

% Flatten loop for parallel execution
parfor idx = 1:totalSteps
    fprintf('  > Simulation %d of %d...\n', idx, totalSteps);

    % Recover indices
    [iPer, iVar] = ind2sub([numFreq, numVar], idx);
    
    % Create temp params struct for this worker
    pTemp = p;
    
    % Restore original amplitude
    pTemp.SquareInput.Pb.amplitude = pOrig.SquareInput.Pb.amplitude;
    
    % Set Period based on Frequency (Period = 1000 ms / Freq)
    currentFreq = testFrequencies(iPer);
    pTemp.SquareInput.Pb.period1 = 1000 / currentFreq;
    pTemp.SquareInput.Pb.period2 = 1000 / currentFreq;
    
    % Set Period Variability (Fraction * Period)
    currentPerVarFrac = testPerVar(iVar);
    pTemp.SquareInput.Pb.perVarPerc1 = currentPerVarFrac * 100;
    pTemp.SquareInput.Pb.perVar1 = currentPerVarFrac * pTemp.SquareInput.Pb.period1;
    pTemp.SquareInput.Pb.perVarPerc2 = currentPerVarFrac * 100;
    pTemp.SquareInput.Pb.perVar2 = currentPerVarFrac * pTemp.SquareInput.Pb.period2;
    pTemp.SquareInput.Pb.toRecompute = true;
    
    fprintf('  > Running Freq: %.2f Hz, Var: %.2f%%\n', ...
        currentFreq, currentPerVarFrac * 100);

    % --- Update Test Name for this parameter combination ---
    % 1. Add Frequency
    strFreq = replace(sprintf('Freq-%gHz', currentFreq), '.', 'p');
    pTemp.testName = virt_updated_test_name_for_rep(pTemp.testName, 'Freq', strFreq);
    
    % 2. Add Variability
    strVar = replace(sprintf('PerVarPerc-%.2f', currentPerVarFrac * 100), '.', 'p');
    pTemp.testName = virt_updated_test_name_for_rep(pTemp.testName, 'PerVar', strVar);

    % --- Run Default Driven Simulation ---
    [~, pOutDriven, ~, analysisDriven, ~, plotDataDriven] = ...
        virt_moore('Params', pTemp, ...
                    'RandomizeSeed', false, ...
                    'RngAlgorithm', rngAlgorithm, ...
                    batchPlotFlags{:}, ...
                    extraArgsCell{:});
    
    % --- Run Replicates (using the helper function) ---
    % Note: virt_moore_multiple_reps runs serially inside this specific parfor worker,
    % but multiple grid points are running in parallel.
    repResultStruct = [];
    if numRandomSims > 0
        % virt_moore_multiple_reps will use pTemp.testName as base and append seed/rep info
        [~, ~, pAll, ~, analysisAll, ~, plotDataAll] = virt_moore_multiple_reps(...
            'Params', pTemp, ...
            'NRepetitions', numRandomSims, ...
            'RepetitionMode', 'MonteCarlo', ...
            'RandomizeSeed', true, ...
            'RngAlgorithm', rngAlgorithm, ...
            batchPlotFlags{:}, ...
            extraArgsCell{:});
        
        % Store temporarily for this index
        for r = 1:numRandomSims
            repResultStruct(r).analysis = analysisAll(r);
            repResultStruct(r).params = pAll(r);
            repResultStruct(r).data = plotDataAll(r);
        end
    end
    
    % Package into a temp struct
    stepRes = struct();
    stepRes.default.analysis = analysisDriven;
    stepRes.default.params = pOutDriven;
    stepRes.default.data = plotDataDriven;
    stepRes.replicates = repResultStruct;
    
    parResults{idx} = stepRes;
end

% Reconstruct the nested output structure
fprintf('Reconstructing output structure...\n');
for idx = 1:totalSteps
    [iPer, iVar] = ind2sub([numFreq, numVar], idx);
    stepRes = parResults{idx};
    
    % Assign Default
    output.drivenDefault.freq(iPer).var(iVar).analysis = stepRes.default.analysis;
    output.drivenDefault.freq(iPer).var(iVar).params = stepRes.default.params;
    output.drivenDefault.freq(iPer).var(iVar).data = stepRes.default.data;
    
    % Assign Replicates
    if numRandomSims > 0
        for r = 1:numRandomSims
            output.drivenRep.freq(iPer).var(iVar).rep(r).analysis = stepRes.replicates(r).analysis;
            output.drivenRep.freq(iPer).var(iVar).rep(r).params = stepRes.replicates(r).params;
            output.drivenRep.freq(iPer).var(iVar).rep(r).data = stepRes.replicates(r).data;
        end
    end
end

fprintf('All simulations complete.\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%