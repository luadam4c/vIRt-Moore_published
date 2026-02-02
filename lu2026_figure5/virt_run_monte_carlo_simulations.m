function [handlesAggregate, aggregatedAnalysis, pAll, analysisAll, ...
            handlesSample, pSample, stateSample, analysisSample, ...
            analysisStringSample, plotDataSample, connSample] = virt_run_monte_carlo_simulations (varargin)
%% Run Monte Carlo simulations of the vIRt model
% Usage: [handlesAggregate, aggregatedAnalysis, pAll, analysisAll, ...
%           handlesSample, pSample, stateSample, analysisSample, ...
%           analysisStringSample, plotDataSample, connSample] = virt_run_monte_carlo_simulations (varargin)
% Explanation:
%       This function executes a two-step Monte Carlo workflow:
%       1. Runs a single representative simulation (`virt_moore.m`).
%       2. Runs the Monte Carlo repetitions (`virt_moore_multiple_reps.m`) with
%          'PlotAggregateOnly' set to true.
%
%       It executes the simulation 30 times (by default) with randomized seeds.
%       By default, figures are not displayed during execution but are saved 
%       to the output directory in .png, .fig, and .epsc formats.
%
% Example(s):
%       % Run with default settings (30 reps, save to pwd subfolders):
%       virt_run_monte_carlo_simulations;
%
%       % Run 30 reps and save to subfolders inside 'C:\MySims':
%       virt_run_monte_carlo_simulations('NRepetitions', 30, 'OutParentDir', 'C:\MySims');
%
% Outputs:
%       Results and figures are saved to the output subdirectories.
%
% Arguments:
%       varargin    - 'NRepetitions': The number of simulations to run.
%                   must be a positive integer scalar
%                   default == 30
%                   - 'ShowFigure': Whether to show figures.
%                   must be a logical scalar
%                   default == false
%                   - 'OutParentDir': The parent directory where sub-folders will be created.
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'CopyScript': Whether to copy the script to the output directory.
%                   must be a logical scalar or empty `[]`. If empty, the script
%                   is copied only when default parameters are used
%                   default == []
%                   - 'SaveOutput': Whether to save figures and data files.
%                   must be a logical scalar
%                   default == true
%                   - 'FigTypes': File formats for saving figures.
%                   must be a cell array of character vectors or string array
%                   default == {'png', 'fig', 'epsc'}
%                   - 'UseParpool': Whether to use parfor for parallel execution.
%                   must be a logical scalar
%                   default == true
%                   - 'RngAlgorithm': The random number generator algorithm to use.
%                       Options: 'twister', 'simdTwister', 'combRecursive',
%                                'multFibonacci', 'philox', 'threefry',
%                                'v4', 'v5uniform', 'v5normal'
%                   must be a string scalar or a character vector
%                   default == [] (set in virt_moore_params.m)
%                   - Any other parameter-value pair for the virt_moore_multiple_reps() function
%                   
% Requires:
%       \Shared\Code\Adams_Functions\archive_dependent_scripts.m
%       \Shared\Code\Adams_Functions\check_dir.m
%       \Shared\Code\vIRt-Moore\virt_moore.m
%       \Shared\Code\vIRt-Moore\virt_moore_multiple_reps.m
%
% Used by:

% File History:
% 2026-01-14 Created by Gemini but updated with desired defaults
% 2026-01-15 Changed nRepetitionsDefault to 50 
% 2026-01-15 Changed nRepetitionsDefault back to 30 
% 2026-01-17 Now uses parfor by default
% 2026-01-18 Added 'RngAlgorithm' optional argument
% 2026-01-18 Modified workflow: Runs virt_moore once with plots, then multiple reps with PlotAggregateOnly.
% 2026-01-18 Fixed input parser to KeepUnmatched=true.

%% Hard-Coded Information For This Script Only
validRngAlgorithms = {'twister', 'simdTwister', 'combRecursive', ...
                      'multFibonacci', 'philox', 'threefry', ...
                      'v4', 'v5uniform', 'v5normal'};

%% Default values for optional arguments
nRepetitionsDefault = 30;
showFigureDefault = false;
outParentDirDefault = pwd;
copyScriptDefault = true;
saveOutputDefault = true;
figTypesDefault = {'png', 'fig', 'epsc'};
useParpoolDefault = true;
rngAlgorithmDefault = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true; % Allowed to pass Params, TestName, etc.

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'NRepetitions', nRepetitionsDefault, @isnumeric);
addParameter(iP, 'ShowFigure', showFigureDefault, ...
    @(x) validateattributes(x, {'logical'}, {'scalar'}));
addParameter(iP, 'OutParentDir', outParentDirDefault, ...
    @(x) ischar(x) || isstring(x));
addParameter(iP, 'CopyScript', copyScriptDefault, ...
    @(x) validateattributes(x, {'logical'}, {'scalar'}));
addParameter(iP, 'SaveOutput', saveOutputDefault, ...
    @(x) validateattributes(x, {'logical'}, {'scalar'}));
addParameter(iP, 'FigTypes', figTypesDefault, ...
    @(x) ischar(x) || isstring(x) || iscell(x));
addParameter(iP, 'UseParpool', useParpoolDefault, ...
    @(x) validateattributes(x, {'logical'}, {'scalar'}));
addParameter(iP, 'RngAlgorithm', rngAlgorithmDefault, ...
    @(x) isempty(x) || any(validatestring(x, validRngAlgorithms)));

% Parse the inputs
parse(iP, varargin{:});
nRepetitions = iP.Results.NRepetitions;
showFigure = iP.Results.ShowFigure;
outParentDir = iP.Results.OutParentDir;
copyScript = iP.Results.CopyScript;
saveOutput = iP.Results.SaveOutput;
figTypes = iP.Results.FigTypes;
useParpool = iP.Results.UseParpool;
rngAlgorithm = iP.Results.RngAlgorithm;

% Keep unmatched arguments for the underlying functions (e.g., Params, TestName)
otherArguments = iP.Unmatched;

%% Preparation
% Manage Output Directory Context
% We temporarily change the current working directory to OutParentDir.
% This allows virt_moore and virt_moore_multiple_reps to default to creating
% their timestamped folders within this parent directory.
if saveOutput
    check_dir(outParentDir); % Ensure parent exists
end

startDir = pwd;
if ~strcmp(startDir, outParentDir) && saveOutput
    cd(outParentDir);
end

% Ensure we return to the original directory when the function finishes or errors
cleanupObj = onCleanup(@() cd(startDir));

% Copy this script and dependencies
if copyScript
    archive_dependent_scripts(mfilename, 'OutFolder', outParentDir);
end

%% Step 1: Run Single Representative Simulation
fprintf('\n=== Step 1: Running single representative simulation ===\n');

% Run virt_moore. 
% We pass 'OutDir' as empty [] so it creates its own 'Test-TimeStamp-Name' 
% folder inside the current directory (which is now OutParentDir).
% We pass the user's 'ShowFigure' preference directly.
[handlesSample, pSample, stateSample, analysisSample, ...
    analysisStringSample, plotDataSample, connSample] = ...
    virt_moore('ShowFigure', showFigure, ...
                'SaveOutput', saveOutput, ...
                'CopyScript', false, ...           % Done above
                'FigTypes', figTypes, ...
                'RngAlgorithm', rngAlgorithm, ...
                otherArguments);

fprintf('Single representative simulation complete.\n');

%% Step 2: Run Monte Carlo Simulations
fprintf('\n=== Step 2: Running Monte Carlo simulations ===\n');
fprintf('Repetitions: %d\n', nRepetitions);

% Run virt_moore_multiple_reps.
% We pass 'OutDir' as empty [] so it creates its own 'MonteCarlo-TimeStamp-Name' 
% folder inside the current directory (which is now OutParentDir).
% We set 'PlotAggregateOnly' to true.
[handlesAggregate, aggregatedAnalysis, pAll, ~, analysisAll, ~, ~] = ...
        virt_moore_multiple_reps(...
            'NRepetitions', nRepetitions, ...
            'RepetitionMode', 'MonteCarlo', ...
            'PlotAggregateOnly', true, ...     % Suppress individual plots
            'ShowFigure', showFigure, ...      % Controls visibility of aggregate plots
            'SaveOutput', saveOutput, ...
            'CopyScript', false, ...           % Done above
            'FigTypes', figTypes, ...
            'UseParpool', useParpool, ...
            'RngAlgorithm', rngAlgorithm, ...
            otherArguments);

fprintf('All Monte Carlo simulations complete.\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%