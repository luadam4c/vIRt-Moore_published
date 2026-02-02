function [handles, aggregatedAnalysis, pAll, stateLast, analysisAll, analysisStringAll, plotDataAll] = ...
            virt_moore_multiple_reps (varargin)
%% Run multiple repetitions of the vIRt model and plot aggregated results
% Usage: [handles, aggregatedAnalysis, pAll, stateLast, analysisAll, analysisStringAll, plotDataAll] = ...
%           virt_moore_multiple_reps (varargin)
% Explanation:
%       This function runs the virt_moore simulation N times. It can operate in
%       two modes: 'MonteCarlo' (default) for simulations with different random
%       seedNumbers, or 'Tuning' to sweep a specific parameter across a range of values.
%       It then aggregates analysis results and generates summary plots.
%       The outputs of the final simulation run are returned to allow for
%       GUI updates.
%
% Outputs:
%       handles             - A structure of handles to all generated figures.
%                           specified as a structure
%       aggregatedAnalysis  - A structure containing the combined data tables.
%                           specified as a structure
%       pAll                - The parameters structures from all simulation runs.
%                           specified as a structure array
%       stateLast           - The state structure from the last simulation run
%                           specified as a structure array
%       analysisAll         - The analysis structures from all runs.
%                           specified as a structure array
%       analysisStringAll   - The analysis summary strings from all runs.
%                           specified as a cell array of string scalars
%       plotDataAll         - The plot data structures from all runs.
%                           specified as a structure array
%
% Arguments:
%       varargin    - 'Params': The base parameters structure for the simulation.
%                   must be a structure
%                   default == Loaded from virt_moore_params.m
%                   - 'Handles': handles to graphics objects
%                   specified as a structure with fields:
%                       fig1, fig2, ax1, ax2, etc.
%                   - 'TestName': The base name for the test run.
%                   must be a string scalar or a character vector
%                   default == '' (managed by virt_moore_params.m)
%                   - 'OutDir': The main output directory for the run.
%                   must be a string scalar or a character vector
%                   default == A timestamped directory in pwd
%                   - 'NRepetitions': The number of simulations to run.
%                   must be a positive integer scalar
%                   default == 10
%                   - 'RepetitionMode': The mode of operation.
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'MonteCarlo' - Vary random seed for each repetition.
%                       'Tuning'     - Vary a specific parameter for each repetition.
%                   default == 'MonteCarlo'
%                   - 'ParamToVary': The path to the parameter to vary in 'Tuning' mode.
%                                       this must start with a '.
%                   must be a string scalar or a character vector
%                   default == ''
%                   - 'RangeToVary': The vector of values to use in 'Tuning' mode.
%                   must be a numeric vector
%                   default == []
%                   - 'StrToVary': A custom string used for figure titles and
%                                  filenames when in 'Tuning' mode. If empty,
%                                  a name is generated from ParamToVary.
%                   must be a string scalar or a character vector
%                   default == ''
%                   - 'RandomizeSeed': Whether to randomize the seed for
%                                       each repetition
%                   must be a logical scalar
%                   default == true
%                   - 'FirstSeedNumber': The seed for the first simulation
%                   must be a numeric scalar or empty `[]`
%                   default == [] (User Input > P.seedNumber > 1)
%                   - 'RngAlgorithm': The random number generator algorithm to use.
%                       Options: 'twister', 'simdTwister', 'combRecursive',
%                                'multFibonacci', 'philox', 'threefry',
%                                'v4', 'v5uniform', 'v5normal'
%                   must be a string scalar or a character vector
%                   default == [] (set in virt_moore_params.m)
%                   - 'MessageHandler': A function handle for displaying progress messages.
%                   must be a function handle
%                   default == [] (prints to command window)
%                   - 'ComputeWhiskAmp': Whether to compute whisk statistics.
%                   must be a logical scalar
%                   default == true
%                   - 'ComputeWhiskInterval': Whether to compute interval statistics.
%                   must be a logical scalar
%                   default == true
%                   - 'SaveOutput': Whether to save figures and data files.
%                   must be a logical scalar
%                   default == true
%                   - 'PlotSampleTraces': Whether to plot sample voltage traces.
%                   must be a logical scalar
%                   default == true
%                   - 'PlotRaster': Whether to plot the spike raster.
%                   must be a logical scalar
%                   default == true
%                   - 'PlotJitter': Whether to plot the whisk log decrement jitter plot.
%                   must be a logical scalar
%                   default == true
%                   - 'PlotScatter': Whether to plot successive whisk amplitude scatter plots.
%                   must be a logical scalar
%                   default == true
%                   - 'PlotPhaseResponse': Whether to plot the phase response curve.
%                   must be a logical scalar
%                   default == true
%                   - 'AlwaysNew': Whether to always create new figures instead of updating.
%                   must be a logical scalar
%                   default == false
%                   - 'ShowFigure': Whether to show figures.
%                   must be a logical scalar
%                   default == true
%                   - 'CopyScript': Whether to copy the script to the output directory.
%                   must be a logical scalar or empty `[]`. If empty, the script
%                   is copied only when default parameters are used
%                   default == []
%                   - 'FigTypes': Figure type(s) for saving.
%                   must be a string scalar, character vector, or a cell array of character vectors
%                   default == {'png'}
%                   - 'UseParpool': Whether to use parfor for parallel execution.
%                   must be a logical scalar
%                   default == false
%
% Requires:
%       \Shared\Code\Adams_Functions\archive_dependent_scripts.m
%       \Shared\Code\Adams_Functions\check_dir.m
%       \Shared\Code\Adams_Functions\combine_tables.m
%       \Shared\Code\Adams_Functions\convert_colors_to_rgb.m
%       \Shared\Code\Adams_Functions\create_time_stamp.m
%       \Shared\Code\Adams_Functions\plot_table.m
%       \Shared\Code\Adams_Functions\force_matrix.m
%       \Shared\Code\Adams_Functions\virt_plot_amplitude_correlation.m
%       \Shared\Code\Adams_Functions\virt_plot_jitter.m
%       \Shared\Code\Adams_Functions\virt_plot_phase_response.m
%       \Shared\Code\Adams_Functions\write_table.m
%       \Shared\Code\vIRt-Moore\virt_moore.m
%       \Shared\Code\vIRt-Moore\virt_moore_params.m
%       \Shared\Code\vIRt-Moore\virt_updated_test_name_for_rep.m
%
% Used by:
%       \Shared\Code\vIRt-Moore\virt_moore_gui.mlapp
%       \Shared\Code\vIRt-Moore\jm_run_virt_sim_replicates.m
%       \Shared\Code\vIRt-Moore\virt_run_monte_carlo_simulations.m
%

% File History:
% 2025-10-02 Created by Gemini.
% 2025-10-03 Modified outputs to return all structures. Added 'SeedNumber'.
% 2025-10-05 Updated basal table structure
% 2025-10-07 Now saves aggregated tables to CSV using write_table.m by Gemini
% 2025-10-07 Added aggregation and plotting of scalar whisk outputs by Gemini
% 2025-10-07 Renamed to virt_moore_multiple_reps.m and added 'Tuning' mode by Gemini.
% 2025-10-08 Changed 'SeedNumber' to 'FirstSeedNumber'
% 2025-10-09 Now passes in StrToVary as an optional argument from the GUI
% 2025-10-09 Now uses the seed number as the grouping variable
% 2025-10-12 Abstracted table combination logic to combine_tables.m
% 2025-10-12 Changed varsToAdd from struct to table for combine_tables.m by Gemini
% 2025-10-15 Updated linked period variability parameters for tuning mode by Gemini
% 2025-10-15 Now plots statistics for each individual simulation as well
% 2025-10-15 Now uses virt_moore_params.m
% 2025-10-15 Changed pooled aggregate statistics figure names to 'pooled'
% 2025-10-15 Added jitter plots for mean statistics across runs by Gemini.
% 2025-10-15 Now aggregates vector outputs as well
% 2025-10-16 Added 'Handles' as an optional argument
% 2025-10-16 Modified by Gemini to store object handles directly and pass them back for plot updates.
% 2025-10-16 Modified by Gemini to add 'FigTypes' as an optional argument for saving plots.
% 2025-10-17 Made 'AlwaysNew' an optional argument, default == false, and pass in from the GUI
% 2025-10-17 Made 'ShowFigure' an optional argument, default == true, and pass in from the GUI
% 2026-01-13 Added 'plotDataAll' output by Gemini
% 2026-01-14 Fixed bug where error bars in PRC plots were not being cleared by Gemini
% 2026-01-14 Added 'PlotAggregateOnly' flag
% 2026-01-14 Default TestName is empty (managed by virt_moore_params.m)
% 2026-01-14 Changed returning stateAll to stateLast to save on memory
% 2026-01-15 Now updates test name with nReps or Tuning info if running from script by Gemini
% 2026-01-16 Added 'UseParpool' optional argument by Gemini
% 2026-01-16 Refactored repetition loops to use common subfunction 'run_single_repetition' by Gemini.
% 2026-01-17 Added 'RngAlgorithm' optional argument by Gemini
% 2026-01-18 Now uses virt_updated_test_name_for_rep.m instead of local subfunction by Gemini
% TODO: Consider linking Rr and Rp parameters like in the GUI

%% Hard-Coded Information For This Script Only
validRepetitionModes = {'MonteCarlo', 'Tuning'};

% Base Figure Titles
figTitleBaseJitterSniffPooled = 'Sniff Start Whisk Log Decrements (Pooled)';
figTitleBaseJitterBasalPooled = 'Basal Respiration Whisk Log Decrements (Pooled)';
figTitleBaseJitterSniffAvg = 'Average Sniff Start Whisk Log Decrements';
figTitleBaseJitterBasalAvg = 'Average Basal Respiration Whisk Log Decrements';
figTitleBaseJitterCorrZScore = 'Basal Respiration Amplitude Correlation Z-Scores';
figTitleBaseScatterSniffPooled = 'Sniff Start Whisk Amplitude Correlations (Pooled)';
figTitleBaseScatterBasalPooled = 'Basal Respiration Whisk Amplitude Correlations (Pooled)';
figTitleBasePRCPooled = 'Whisk Phase Response Curve (Pooled)';
figTitleBaseWhiskScalar = 'Aggregated Whisk Scalar Outputs';

% Base Figure Names (for saving)
figNameBaseJitterSniffPooled = 'sniffstart_whisk_log_decrements_pooled_jitter';
figNameBaseJitterBasalPooled = 'basalresp_whisk_log_decrements_pooled_jitter';
figNameBaseJitterSniffAvg = 'sniffstart_whisk_log_decrements_avg_jitter';
figNameBaseJitterBasalAvg = 'basalresp_whisk_log_decrements_avg_jitter';
figNameBaseJitterBasalAmpCorrZScore = 'basalresp_whisk_corr_zscore_jitter';
figNameBaseScatterSniffPooled = 'sniffstart_whisk_amplitudes_pooled_scatter';
figNameBaseScatterBasalPooled = 'basalresp_whisk_amplitudes_pooled_scatter';
figNameBasePRCPooled = 'phase_response_pooled_scatter';
figNameBaseWhiskScalarSummary = 'aggregated_whisk_scalar_outputs_summary';

% Base File Names for aggregated data tables
fileNameBaseSniffTable = 'aggregated_sniffstart_table';
fileNameBaseBasalTable = 'aggregated_basalresp_table';
fileNameBaseWhiskScalarTable = 'aggregated_whisk_scalar_outputs';
fileNameBaseWhiskVectorTable = 'aggregated_whisk_vector_outputs';

% Output variables to plot
varsToPlot = {'fundamentalFrequency', 'meanFirstPeakAmplitudeToAnalyze'};
varLabels = {'Whisk Fundamental Frequency (Hz)', 'Average First Whisk Amplitude (a.u.)'};

% Default randomization parameters
defaultFirstSeedNumber = 1;         % Default seed is 1 if there is no user input
validRngAlgorithms = {'twister', 'simdTwister', 'combRecursive', ...
                      'multFibonacci', 'philox', 'threefry', ...
                      'v4', 'v5uniform', 'v5normal'};

%% Default values for optional arguments
paramsDefault = [];
handlesDefault = [];            % no handles by default
testNameDefault = '';
dirOutDefault = '';
nRepetitionsDefault = 10;
repetitionModeDefault = 'MonteCarlo';
paramToVaryDefault = '';
rangeToVaryDefault = [];
strToVaryDefault = '';              % Set later
toRandomizeSeedDefault = [];        % Set later
firstSeedNumberDefault = [];
rngAlgorithmDefault = [];           % Set later
messageHandlerDefault = [];
alwaysNewDefault = false;
showFigureDefault = true;
copyScriptDefault = true;
plotSampleTracesDefault = true;
plotRasterDefault = true;
plotJitterDefault = true;
plotScatterDefault = true;
plotPRCDefault = true;
computeWhiskAmpDefault = true;
computeWhiskIntervalDefault = true;
saveOutputDefault = true;
figTypesDefault = {'png'};
plotAggregateOnlyDefault = false;
useParpoolDefault = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set up Input Parser
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Params', paramsDefault, @isstruct);
addParameter(iP, 'Handles', handlesDefault);
addParameter(iP, 'TestName', testNameDefault, @ischar);
addParameter(iP, 'OutDir', dirOutDefault, @ischar);
addParameter(iP, 'NRepetitions', nRepetitionsDefault, @isnumeric);
addParameter(iP, 'RepetitionMode', repetitionModeDefault, ...
    @(x) ~isempty(validatestring(x, validRepetitionModes)));
addParameter(iP, 'ParamToVary', paramToVaryDefault, @ischar);
addParameter(iP, 'RangeToVary', rangeToVaryDefault, @isnumeric);
addParameter(iP, 'StrToVary', strToVaryDefault, @ischar);
addParameter(iP, 'RandomizeSeed', toRandomizeSeedDefault, ...
    @(x) validateattributes(x, {'logical'}, {'scalar'}));
addParameter(iP, 'FirstSeedNumber', firstSeedNumberDefault, ...
    @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
addParameter(iP, 'RngAlgorithm', rngAlgorithmDefault, ...
    @(x) isempty(x) || any(validatestring(x, validRngAlgorithms)));
addParameter(iP, 'MessageHandler', messageHandlerDefault);
addParameter(iP, 'AlwaysNew', alwaysNewDefault, ...
    @(x) validateattributes(x, {'logical'}, {'scalar'}));
addParameter(iP, 'ShowFigure', showFigureDefault, ...
    @(x) validateattributes(x, {'logical'}, {'scalar'}));
addParameter(iP, 'CopyScript', copyScriptDefault, ...
    @(x) validateattributes(x, {'logical'}, {'scalar'}));
addParameter(iP, 'PlotSampleTraces', plotSampleTracesDefault, ...
    @(x) validateattributes(x, {'logical'}, {'scalar'}));
addParameter(iP, 'PlotRaster', plotRasterDefault, ...
    @(x) validateattributes(x, {'logical'}, {'scalar'}));
addParameter(iP, 'PlotJitter', plotJitterDefault, ...
    @(x) validateattributes(x, {'logical'}, {'scalar'}));
addParameter(iP, 'PlotScatter', plotScatterDefault, ...
    @(x) validateattributes(x, {'logical'}, {'scalar'}));
addParameter(iP, 'PlotPhaseResponse', plotPRCDefault, ...
    @(x) validateattributes(x, {'logical'}, {'scalar'}));
addParameter(iP, 'ComputeWhiskAmp', computeWhiskAmpDefault, ...
    @(x) validateattributes(x, {'logical'}, {'scalar'}));
addParameter(iP, 'ComputeWhiskInterval', computeWhiskIntervalDefault, ...
    @(x) validateattributes(x, {'logical'}, {'scalar'}));
addParameter(iP, 'SaveOutput', saveOutputDefault, ...
    @(x) validateattributes(x, {'logical'}, {'scalar'}));
addParameter(iP, 'FigTypes', figTypesDefault, ...
    @(x) ischar(x) || isstring(x) || iscell(x));
addParameter(iP, 'PlotAggregateOnly', plotAggregateOnlyDefault, ...
    @(x) validateattributes(x, {'logical'}, {'scalar'}));
addParameter(iP, 'UseParpool', useParpoolDefault, ...
    @(x) validateattributes(x, {'logical'}, {'scalar'}));


% Parse the inputs
parse(iP, varargin{:});
pBase = iP.Results.Params;
handles = iP.Results.Handles;
testName = iP.Results.TestName;
dirOut = iP.Results.OutDir;
nRepetitions = iP.Results.NRepetitions;
repetitionMode = validatestring(iP.Results.RepetitionMode, {'MonteCarlo', 'Tuning'});
paramToVary = iP.Results.ParamToVary;
rangeToVary = iP.Results.RangeToVary;
strToVary = iP.Results.StrToVary;
toRandomizeSeed = iP.Results.RandomizeSeed;
firstSeedNumber = iP.Results.FirstSeedNumber;
rngAlgorithmUser = iP.Results.RngAlgorithm;
messageHandler = iP.Results.MessageHandler;
alwaysNew = iP.Results.AlwaysNew;
showFigure = iP.Results.ShowFigure;
copyScript = iP.Results.CopyScript;
plotSampleTraces = iP.Results.PlotSampleTraces;
plotRaster = iP.Results.PlotRaster;
plotJitter = iP.Results.PlotJitter;
plotScatter = iP.Results.PlotScatter;
plotPRC = iP.Results.PlotPhaseResponse;
computeWhiskAmp = iP.Results.ComputeWhiskAmp;
computeWhiskInterval = iP.Results.ComputeWhiskInterval;
saveOutput = iP.Results.SaveOutput;
figTypes = iP.Results.FigTypes;
plotAggregateOnly = iP.Results.PlotAggregateOnly;
useParpool = iP.Results.UseParpool;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preparation
% Define a local message handler for convenience
post_message = @(msg) local_post_message(msg, messageHandler);

% Display message
if isempty(pBase)
    post_message('No parameters provided, loading defaults from virt_moore...');
end

% Create or check the parameters structure
pBase = virt_moore_params('ParamsIn', pBase, 'SeedNumber', firstSeedNumber, ...
                            'RngAlgorithm', rngAlgorithmUser, ...
                            'TestName', testName, 'MessageHandler', messageHandler);

% Update the test name from the parameters structure
testName = pBase.testName;

% Update test name with simulation info if not already present
% If running from script (not GUI), the test name might not have the sweep info.
% We update it here to ensure the output directory has a descriptive name.

% 1. Append nReps if multiple repetitions and not already in name
nRepsPattern = '-nReps-[0-9]+';
if isempty(regexp(testName, nRepsPattern, 'once')) && nRepetitions > 1
    testName = sprintf('%s-nReps-%d', testName, nRepetitions);
    % Update P.testName as well
    pBase.testName = testName;
end

% 2. Append Tuning Parameter info if in Tuning mode and not already in name
if strcmp(repetitionMode, 'Tuning') && ~isempty(paramToVary)
    cleanName = clean_param_name(paramToVary);
    
    % Check if this parameter is already in the test name
    % Look for "-cleanName-" pattern
    paramPattern = ['-', cleanName, '-'];
    
    if ~contains(testName, paramPattern)
        % Construct a suffix string based on the range
        if isnumeric(rangeToVary) && numel(rangeToVary) > 1
            % e.g. "0to100"
            valStr = sprintf('%gto%g', min(rangeToVary), max(rangeToVary));
        elseif isnumeric(rangeToVary)
            valStr = sprintf('%g', rangeToVary);
        else
            valStr = 'varied';
        end
        
        % Replace decimal points with 'p' for cleaner filenames
        valStr = replace(valStr, '.', 'p');
        
        testName = sprintf('%s-%s-%s', testName, cleanName, valStr);
        pBase.testName = testName;
    end
end

% Set up the main output directory for the run
if isempty(dirOut)
    % Get current timestamp
    timeStamp = create_time_stamp;

    % Create directory name
    dirOut = fullfile(pwd, sprintf('%s-%s-%s', repetitionMode, timeStamp, testName));
end

% Create the directory if it doesn't exist
if saveOutput
    check_dir(dirOut);
end 

% Copy this script and its dependencies to the output directory if requested
if copyScript && saveOutput
    archive_dependent_scripts(mfilename, 'OutFolder', dirOut);
end

% Decide on whether to randomize the seed
if isempty(toRandomizeSeed)
    switch repetitionMode
    case 'MonteCarlo'
        % Monte Carlo simulations must have seed randomized
        toRandomizeSeed = true;
    case 'Tuning'
        % By default seed is unchanged when performing parameter sweeps
        toRandomizeSeed = false;
    end
end

% No seed provided by user, use the seed number stroed in parameters or the default
if isempty(firstSeedNumber)
    if isfield(pBase, 'seedNumber')
        firstSeedNumber = pBase.seedNumber;
        post_message(sprintf('Using seed number from parameter file as first seed number: %d', firstSeedNumber));
    else
        firstSeedNumber = defaultFirstSeedNumber;
        post_message(sprintf('Using default first seed number: %d', firstSeedNumber));
    end
end

% Set all seedNumbers to the first seed number first
seedsToUse = firstSeedNumber * ones(nRepetitions, 1);

% Randomize the rest of the seedNumbers if requested and needed
if toRandomizeSeed && nRepetitions > 1
    % Use the first seed number to seed the random number generator
    % Pass the chosen algorithm here as well
    rng(firstSeedNumber, pBase.rngAlgorithm);

    % Generate the rest of the seedNumbers
    seedsToUse(2:end) = randi(2^32 - 1, nRepetitions - 1, 1);
end

% Decide on parameter to vary and range to vary
switch repetitionMode
case 'MonteCarlo'
    % Create warning
    if ~isempty(paramToVary) || ~isempty(rangeToVary)
        warning('ParamToVary and RangeToVary are set to seed numbers in Monte Carlo mode.');
    end

    % Set paramToVary and rangeToVary
    paramToVary = '.seedNumber';
    rangeToVary = seedsToUse;
case 'Tuning'
    % Throw error
    if isempty(paramToVary) || isempty(rangeToVary)
        error('ParamToVary and RangeToVary must be passed in for Tuning mode.');
    end
end

% Clean up parameter name for file names
cleanParamName = clean_param_name(paramToVary);

% Create string to vary if not passed in
switch repetitionMode
case 'MonteCarlo'
    % Create warning
    if ~isempty(strToVary)
        warning('StrToVary is ignored in Monte Carlo mode.');
    end

    % Set strToVary
    strToVary = 'MonteCarlo';
case 'Tuning'
    if isempty(strToVary)
        strToVary = sprintf('%s-varied', cleanParamName);
    end
end

% Create figure title suffix and file and figure name suffix
titleSuffix = sprintf('(%s)', replace(strToVary, '_', ' '));
fileSuffix = sprintf('_%s', strToVary);

% Construct full titles and file names
figTitleJitterSniff = [figTitleBaseJitterSniffPooled, ' ', titleSuffix];
figTitleJitterBasal = [figTitleBaseJitterBasalPooled, ' ', titleSuffix];
figTitleJitterSniffAvg = [figTitleBaseJitterSniffAvg, ' ', titleSuffix];
figTitleJitterBasalAvg = [figTitleBaseJitterBasalAvg, ' ', titleSuffix];
figTitleJitterBasalAmpCorrZScore = [figTitleBaseJitterCorrZScore, ' ', titleSuffix];
figTitleScatterSniff = [figTitleBaseScatterSniffPooled, ' ', titleSuffix];
figTitleScatterBasal = [figTitleBaseScatterBasalPooled, ' ', titleSuffix];
figTitlePRC = [figTitleBasePRCPooled, ' ', titleSuffix];
figTitleWhiskScalar = [figTitleBaseWhiskScalar, ' ', titleSuffix];

figNameJitterSniff = [figNameBaseJitterSniffPooled, fileSuffix];
figNameJitterBasal = [figNameBaseJitterBasalPooled, fileSuffix];
figNameJitterSniffAvg = [figNameBaseJitterSniffAvg, fileSuffix];
figNameJitterBasalAvg = [figNameBaseJitterBasalAvg, fileSuffix];
figNameJitterBasalAmpCorrZScore = [figNameBaseJitterBasalAmpCorrZScore, fileSuffix];
figNameScatterSniff = [figNameBaseScatterSniffPooled, fileSuffix];
figNameScatterBasal = [figNameBaseScatterBasalPooled, fileSuffix];
figNamePRC = [figNameBasePRCPooled, fileSuffix];
figNameWhiskScalarSummary = [figNameBaseWhiskScalarSummary, fileSuffix];

fileNameSniffTable = [fileNameBaseSniffTable, fileSuffix, '.csv'];
fileNameBasalTable = [fileNameBaseBasalTable, fileSuffix, '.csv'];
fileNameWhiskScalarTable = [fileNameBaseWhiskScalarTable, fileSuffix, '.csv'];
fileNameWhiskVectorTable = [fileNameBaseWhiskVectorTable, fileSuffix, '.csv'];

% Convert all color strings to RGB vectors within pBase.Plotting and save as pPlot
pPlot = convert_colors_to_rgb(pBase.Plotting, 'ColorSubStr', 'color');

% Extract analysis parameters
maxWhisksBasalRespToAnalyze = pBase.Analysis.maxWhisksBasalRespToAnalyze;
nCorrToAnalyze = pBase.Analysis.nCorrToAnalyze;
whiskDirForPhase = pBase.Analysis.whiskDirForPhase;

% Handle PlotAggregateOnly logic
if plotAggregateOnly
    % If "Aggregate Only" is true, suppress individual plots and figure showing
    repPlotSampleTraces = false;
    repPlotRaster = false;
    repPlotJitter = false;
    repPlotScatter = false;
    repPlotPRC = false;
    repShowFigure = false;
else
    % Otherwise, use the global settings passed in by the user
    repPlotSampleTraces = plotSampleTraces;
    repPlotRaster = plotRaster;
    repPlotJitter = plotJitter;
    repPlotScatter = plotScatter;
    repPlotPRC = plotPRC;
    repShowFigure = showFigure;
end

%% Initialize Outputs
% Pre-allocate cell arrays for all outputs from each run
pAll = cell(nRepetitions, 1);
analysisAllCell = cell(nRepetitions, 1);
analysisStringAll = cell(nRepetitions, 1);
plotDataAllCell = cell(nRepetitions, 1);
stateLast = struct;

% Bundle individual plot settings into a struct to pass to the helper function
simPlotSettings = struct( ...
    'PlotSampleTraces', repPlotSampleTraces, ...
    'PlotRaster', repPlotRaster, ...
    'PlotJitter', repPlotJitter, ...
    'PlotScatter', repPlotScatter, ...
    'PlotPhaseResponse', repPlotPRC);

%% Run Simulations
post_message(sprintf('Starting %s simulation with %d repetitions...', repetitionMode, nRepetitions));

if useParpool
    % --- PARALLEL EXECUTION (PARFOR) ---
    % When running in parallel, we cannot update a single shared figure handle set.
    % We pass 'Handles' as [] to each worker.
    
    % Create cell array to store 'stateLast' just in the last cell
    stateLastCell = cell(nRepetitions, 1);
    
    fprintf('Running in parallel using parfor. Real-time plot updates disabled.\n');
    
    parfor iRep = 1:nRepetitions
        % Note: Cannot use 'post_message' here safely if it involves GUI handles
        fprintf('--- Repetition %d of %d (Parallel) ---\n', iRep, nRepetitions);
        
        % Run the single repetition logic
        % Note: Handles is [], MessageHandler is []
        [~, pRep, stateRep, analysisRep, analysisStringRep, plotDataRep] = ...
            run_single_repetition(iRep, pBase, seedsToUse, paramToVary, rangeToVary, ...
                cleanParamName, testName, dirOut, [], [], saveOutput, ...
                computeWhiskAmp, computeWhiskInterval, figTypes, simPlotSettings, repShowFigure);
        
        % To save on memory, only store state for the last repetition
        if iRep == nRepetitions
            stateLastCell{iRep} = stateRep;
        else
            stateLastCell{iRep} = [];
        end

        % Store outputs
        pAll{iRep} = pRep;
        analysisAllCell{iRep} = analysisRep;
        analysisStringAll{iRep} = analysisStringRep;
        plotDataAllCell{iRep} = plotDataRep;
    end
    
    % Retrieve stateLast from the last repetition
    if ~isempty(stateLastCell)
        stateLast = stateLastCell{end};
    end
    
else
    % --- SERIAL EXECUTION ---
    % Loop through each repetition
    handlesPrev = handles;
    
    for iRep = 1:nRepetitions
        % Progress message
        post_message(sprintf('--- Repetition %d of %d ---', iRep, nRepetitions));
        
        % Run the single repetition logic
        % Note: Passes handlesPrev and messageHandler
        [handlesRep, pRep, stateRep, analysisRep, analysisStringRep, plotDataRep] = ...
            run_single_repetition(iRep, pBase, seedsToUse, paramToVary, rangeToVary, ...
                cleanParamName, testName, dirOut, handlesPrev, messageHandler, saveOutput, ...
                computeWhiskAmp, computeWhiskInterval, figTypes, simPlotSettings, repShowFigure);

        % Store all outputs from the current repetition
        %   Note: Only the lastest State structure is saved to free up memory
        pAll{iRep} = pRep;
        stateLast = stateRep;
        analysisAllCell{iRep} = analysisRep;
        analysisStringAll{iRep} = analysisStringRep;
        plotDataAllCell{iRep} = plotDataRep;

        % Decide whether to use the same handles for the next round of plots
        if alwaysNew
            handlesPrev = [];
        else
            handlesPrev = handlesRep;
        end

        % If this is the last repetition, save its handles to be returned
        if iRep == nRepetitions
            handles = handlesRep;
        end
    end
end

% Convert cell arrays of structs to structure arrays for output
if ~isempty(pAll) && all(cellfun(@isstruct, pAll))
    pAll = vertcat(pAll{:});
end
if ~isempty(analysisAllCell) && all(cellfun(@isstruct, analysisAllCell))
    analysisAll = vertcat(analysisAllCell{:});
else
    analysisAll = analysisAllCell; % Return as cell if empty or contains non-structs
end
if ~isempty(plotDataAllCell) && all(cellfun(@isstruct, plotDataAllCell))
    plotDataAll = vertcat(plotDataAllCell{:});
else
    plotDataAll = plotDataAllCell;
end

%% Aggregate Data
post_message('All repetitions complete. Aggregating and plotting final results...');

% Gather Sniff Start Windows and Basal Respiration Cycles from all runs
allSniffTables = cellfun(@(x) x.whisk.sniffStartWinTable, analysisAllCell, 'UniformOutput', false);
allBasalTables = cellfun(@(x) x.whisk.basalRespCycleTable, analysisAllCell, 'UniformOutput', false);

% Create a table with identifying variables to add to the aggregated tables
repetitionNumber = (1:nRepetitions)';
seedNumber = seedsToUse;
paramValues = rangeToVary(:);
varsToAdd = table(repetitionNumber, seedNumber);
varsToAdd.(cleanParamName) = paramValues; % Add the varied parameter as a column

% Aggregate Sniff Start Windows and Basal Respiration Cycles into single tables
aggregatedSniffTable = combine_tables(allSniffTables, 'VarsToAdd', varsToAdd);
aggregatedBasalTable = combine_tables(allBasalTables, 'VarsToAdd', varsToAdd);

% Aggregate scalar whisk outputs from all repetitions
aggregatedWhiskScalarTable = aggregate_scalar_outputs(analysisAll, cleanParamName, rangeToVary, seedsToUse);

% Aggregate vector whisk outputs from all repetitions for jitter plots
aggregatedWhiskVectorTable = aggregate_vector_outputs(analysisAll, cleanParamName, rangeToVary, seedsToUse);

% Save in output
aggregatedAnalysis.sniffTable = aggregatedSniffTable;
aggregatedAnalysis.basalTable = aggregatedBasalTable;
aggregatedAnalysis.whiskScalarOutputs = aggregatedWhiskScalarTable;
aggregatedAnalysis.whiskVectorOutputs = aggregatedWhiskVectorTable;

% Save aggregated tables to .csv files if requested
if saveOutput
    post_message('Saving aggregated analysis tables...');

    % Define output file paths
    pathSniffTable = fullfile(dirOut, fileNameSniffTable);
    pathBasalTable = fullfile(dirOut, fileNameBasalTable);
    pathWhiskScalarTable = fullfile(dirOut, fileNameWhiskScalarTable);
    pathWhiskVectorTable = fullfile(dirOut, fileNameWhiskVectorTable);

    % Save the aggregated sniff start window table
    write_table(aggregatedSniffTable, pathSniffTable);
    post_message(sprintf('Aggregated Sniff Start Window table saved to %s!', pathSniffTable));

    % Save the aggregated basal respiration cycle table
    write_table(aggregatedBasalTable, pathBasalTable);
    post_message(sprintf('Aggregated Basal Respiration Cycle table saved to %s!', pathBasalTable));

    % Save the aggregated whisk scalar outputs table
    write_table(aggregatedWhiskScalarTable, pathWhiskScalarTable);
    post_message(sprintf('Aggregated Whisk Scalar table saved to %s!', pathWhiskScalarTable));

    % Save the aggregated whisk vector outputs table
    write_table(aggregatedWhiskVectorTable, pathWhiskVectorTable);
    post_message(sprintf('Aggregated Whisk Vector table saved to %s!', pathWhiskVectorTable));
end

%% Plot aggregated results
% --- FOR GUI: Reset figure state to 'Rendering' before updating ---
if ~isempty(handles)
    allHandles = struct2cell(handles);
    for i = 1:numel(allHandles)
        currentHandleArray = allHandles{i};
        for j = 1:numel(currentHandleArray)
            currentHandle = currentHandleArray(j);
            if isgraphics(currentHandle) && strcmp(get(currentHandle, 'Type'), 'figure')
                set(currentHandle, 'UserData', 'Rendering');
            end
        end
    end
end

% Plot Amplitude Ratio Log Decrement plots for aggregated data
if plotJitter
    % --- Jitter plots of all individual data points (pooled) ---
    % Decide whether to create or update the sniff pooled jitter plot
    if ~alwaysNew && isfield(handles, 'figJitterSniffPooled') && ...
            isgraphics(handles.figJitterSniffPooled)
        post_message('Updating sniff start pooled jitter plot...');
        handlesToUpdate.fig = handles.figJitterSniffPooled;
        handlesToUpdate.axJitter = handles.axJitterSniffPooled;
        handlesToUpdate.hJitter = handles.hJitterSniffPooled;
        handlesToUpdate.hErrorBars = handles.hErrorBarsSniffPooled;
        handlesToUpdate.hMeans = handles.hMeansSniffPooled;
        handlesToUpdate.pTextJitter = handles.pTextJitterSniffPooled;
        handlesToUpdate.hNull = handles.hNullSniffPooled;
        handlesToUpdate.sigMarkerJitter = handles.sigMarkerJitterSniffPooled;
        handlesToUpdate.transformedLabels = handles.transformedLabelsSniffPooled;
    else
        post_message('Plotting sniff start pooled jitter plot...');
        handlesToUpdate = [];
    end
    [~, hOut] = virt_plot_jitter(aggregatedSniffTable, pPlot, ...
                            'DataColumn', 'whiskLogDecrements', ...
                            'DataMode', 'LogDecrement', ...
                            'GroupingColumn', cleanParamName, ...
                            'Handles', handlesToUpdate, ...
                            'FigTitle', figTitleJitterSniff, ...
                            'FigName', figNameJitterSniff, ...
                            'OutDir', dirOut, 'ToSaveOutput', saveOutput, ...
                            'FigTypes', figTypes, 'ShowFigure', showFigure);
    if isfield(hOut, 'axJitter')
        handles.figJitterSniffPooled = hOut.fig;
        handles.axJitterSniffPooled = hOut.axJitter;
        handles.hJitterSniffPooled = hOut.hJitter;
        handles.hErrorBarsSniffPooled = hOut.hErrorBars;
        handles.hMeansSniffPooled = hOut.hMeans;
        handles.hNullSniffPooled = hOut.hNull;
        handles.pTextJitterSniffPooled = hOut.pTextJitter;
        handles.sigMarkerJitterSniffPooled = hOut.sigMarkerJitter;
        handles.transformedLabelsSniffPooled = hOut.transformedLabels;
    end

    % Decide whether to create or update the basal pooled jitter plot
    if ~alwaysNew && isfield(handles, 'figJitterBasalPooled') && ...
            isgraphics(handles.figJitterBasalPooled)
        post_message('Updating basal respiration pooled jitter plot...');
        handlesToUpdate.fig = handles.figJitterBasalPooled;
        handlesToUpdate.axJitter = handles.axJitterBasalPooled;
        handlesToUpdate.hJitter = handles.hJitterBasalPooled;
        handlesToUpdate.hErrorBars = handles.hErrorBarsBasalPooled;
        handlesToUpdate.hMeans = handles.hMeansBasalPooled;
        handlesToUpdate.pTextJitter = handles.pTextJitterBasalPooled;
        handlesToUpdate.hNull = handles.hNullBasalPooled;
        handlesToUpdate.sigMarkerJitter = handles.sigMarkerJitterBasalPooled;
        handlesToUpdate.transformedLabels = handles.transformedLabelsBasalPooled;
    else
        post_message('Plotting basal respiration pooled jitter plot...');
        handlesToUpdate = [];
    end
    [~, hOut] = virt_plot_jitter(aggregatedBasalTable, pPlot, ...
                            'DataColumn', 'whiskLogDecrements', ...
                            'DataMode', 'LogDecrement', ...
                            'GroupingColumn', cleanParamName, ...
                            'Handles', handlesToUpdate, ...
                            'MaxOrders', maxWhisksBasalRespToAnalyze - 1, ...
                            'FigTitle', figTitleJitterBasal, ...
                            'FigName', figNameJitterBasal, ...
                            'OutDir', dirOut, 'ToSaveOutput', saveOutput, ...
                            'FigTypes', figTypes, 'ShowFigure', showFigure);
    if isfield(hOut, 'axJitter')
        handles.figJitterBasalPooled = hOut.fig;
        handles.axJitterBasalPooled = hOut.axJitter;
        handles.hJitterBasalPooled = hOut.hJitter;
        handles.hErrorBarsBasalPooled = hOut.hErrorBars;
        handles.hMeansBasalPooled = hOut.hMeans;
        handles.hNullBasalPooled = hOut.hNull;
        handles.pTextJitterBasalPooled = hOut.pTextJitter;
        handles.sigMarkerJitterBasalPooled = hOut.sigMarkerJitter;
        handles.transformedLabelsBasalPooled = hOut.transformedLabels;
    end
    
    % --- Jitter plots of average statistics across runs ---
    % Decide whether to create or update the sniff average jitter plot
    if ~alwaysNew && isfield(handles, 'figJitterSniffAvg') && ...
            isgraphics(handles.figJitterSniffAvg)
        post_message('Updating sniff start average jitter plot...');
        handlesToUpdate.fig = handles.figJitterSniffAvg;
        handlesToUpdate.axJitter = handles.axJitterSniffAvg;
        handlesToUpdate.hJitter = handles.hJitterSniffAvg;
        handlesToUpdate.hErrorBars = handles.hErrorBarsSniffAvg;
        handlesToUpdate.hMeans = handles.hMeansSniffAvg;
        handlesToUpdate.pTextJitter = handles.pTextJitterSniffAvg;
        handlesToUpdate.hNull = handles.hNullSniffAvg;
        handlesToUpdate.sigMarkerJitter = handles.sigMarkerJitterSniffAvg;
        handlesToUpdate.transformedLabels = handles.transformedLabelsSniffAvg;
    else
        post_message('Plotting sniff start average jitter plot...');
        handlesToUpdate = [];
    end
    [~, hOut] = virt_plot_jitter(aggregatedWhiskVectorTable, pPlot, ...
                    'DataColumn', 'meanLogDecrementsSniff', 'DataMode', 'LogDecrement', ...
                    'GroupingColumn', cleanParamName, ...
                    'Handles', handlesToUpdate, ...
                    'FigTitle', figTitleJitterSniffAvg, ...
                    'FigName', figNameJitterSniffAvg, ...
                    'OutDir', dirOut, 'ToSaveOutput', saveOutput, ...
                    'FigTypes', figTypes, 'ShowFigure', showFigure);
    if isfield(hOut, 'axJitter')
        handles.figJitterSniffAvg = hOut.fig;
        handles.axJitterSniffAvg = hOut.axJitter;
        handles.hJitterSniffAvg = hOut.hJitter;
        handles.hErrorBarsSniffAvg = hOut.hErrorBars;
        handles.hMeansSniffAvg = hOut.hMeans;
        handles.hNullSniffAvg = hOut.hNull;
        handles.pTextJitterSniffAvg = hOut.pTextJitter;
        handles.sigMarkerJitterSniffAvg = hOut.sigMarkerJitter;
        handles.transformedLabelsSniffAvg = hOut.transformedLabels;
    end

    % Decide whether to create or update the basal average jitter plot
    if ~alwaysNew && isfield(handles, 'figJitterBasalAvg') && ...
            isgraphics(handles.figJitterBasalAvg)
        post_message('Updating basal respiration average jitter plot...');
        handlesToUpdate.fig = handles.figJitterBasalAvg;
        handlesToUpdate.axJitter = handles.axJitterBasalAvg;
        handlesToUpdate.hJitter = handles.hJitterBasalAvg;
        handlesToUpdate.hErrorBars = handles.hErrorBarsBasalAvg;
        handlesToUpdate.hMeans = handles.hMeansBasalAvg;
        handlesToUpdate.pTextJitter = handles.pTextJitterBasalAvg;
        handlesToUpdate.hNull = handles.hNullBasalAvg;
        handlesToUpdate.sigMarkerJitter = handles.sigMarkerJitterBasalAvg;
        handlesToUpdate.transformedLabels = handles.transformedLabelsBasalAvg;
    else
        post_message('Plotting basal respiration average jitter plot...');
        handlesToUpdate = [];
    end
     [~, hOut] = virt_plot_jitter(aggregatedWhiskVectorTable, pPlot, ...
                    'DataColumn', 'meanLogDecrementsBasal', 'DataMode', 'LogDecrement', ...
                    'GroupingColumn', cleanParamName, ...
                    'Handles', handlesToUpdate, ...
                    'MaxOrders', maxWhisksBasalRespToAnalyze - 1, ...
                    'FigTitle', figTitleJitterBasalAvg, ...
                    'FigName', figNameJitterBasalAvg, ...
                    'OutDir', dirOut, 'ToSaveOutput', saveOutput, ...
                    'FigTypes', figTypes, 'ShowFigure', showFigure);
    if isfield(hOut, 'axJitter')
        handles.figJitterBasalAvg = hOut.fig;
        handles.axJitterBasalAvg = hOut.axJitter;
        handles.hJitterBasalAvg = hOut.hJitter;
        handles.hErrorBarsBasalAvg = hOut.hErrorBars;
        handles.hMeansBasalAvg = hOut.hMeans;
        handles.hNullBasalAvg = hOut.hNull;
        handles.pTextJitterBasalAvg = hOut.pTextJitter;
        handles.sigMarkerJitterBasalAvg = hOut.sigMarkerJitter;
        handles.transformedLabelsBasalAvg = hOut.transformedLabels;
    end
end

% Plot Amplitude Correlation plots for aggregated data
if plotScatter
    % Decide whether to create or update the sniff scatter plot
    if ~alwaysNew && isfield(handles, 'figScatterSniffMC') && ...
            isgraphics(handles.figScatterSniffMC)
        post_message('Updating sniff start pooled scatter plot...');
        handlesToUpdate.fig = handles.figScatterSniffMC;
        handlesToUpdate.axScatter = handles.axScatterSniffMC;
        handlesToUpdate.hScatters = handles.hScattersSniffMC;
        handlesToUpdate.hBestFits = handles.hBestFitsSniffMC;
        handlesToUpdate.hTextBestFit = handles.hTextBestFitSniffMC;
        handlesToUpdate.hThrOrigs = handles.hThrOrigsSniffMC;
        handlesToUpdate.hTextThrOrig = handles.hTextThrOrigSniffMC;
    else
        post_message('Plotting sniff start pooled scatter plot...');
        handlesToUpdate = [];
    end
    [~, hOut] = virt_plot_amplitude_correlation(aggregatedSniffTable, pPlot, ...
                            'GroupingColumn', cleanParamName, ...
                            'Handles', handlesToUpdate, ...
                            'NCorrelations', nCorrToAnalyze, ...
                            'FigTitle', figTitleScatterSniff, ...
                            'FigName', figNameScatterSniff, ...
                            'OutDir', dirOut, 'ToSaveOutput', saveOutput, ...
                            'FigTypes', figTypes, 'ShowFigure', showFigure);
    if isfield(hOut, 'axScatter')
        handles.figScatterSniffMC = hOut.fig;
        handles.axScatterSniffMC = hOut.axScatter;
        handles.hScattersSniffMC = hOut.hScatters;
        handles.hBestFitsSniffMC = hOut.hBestFits;
        handles.hTextBestFitSniffMC = hOut.hTextBestFit;
        handles.hThrOrigsSniffMC = hOut.hThrOrigs;
        handles.hTextThrOrigSniffMC = hOut.hTextThrOrig;
    end

    % Decide whether to create or update the basal scatter plot
    if ~alwaysNew && isfield(handles, 'figScatterBasalMC') && ...
            isgraphics(handles.figScatterBasalMC)
        post_message('Updating basal respiration pooled scatter plot...');
        handlesToUpdate.fig = handles.figScatterBasalMC;
        handlesToUpdate.axScatter = handles.axScatterBasalMC;
        handlesToUpdate.hScatters = handles.hScattersBasalMC;
        handlesToUpdate.hBestFits = handles.hBestFitsBasalMC;
        handlesToUpdate.hTextBestFit = handles.hTextBestFitBasalMC;
        handlesToUpdate.hThrOrigs = handles.hThrOrigsBasalMC;
        handlesToUpdate.hTextThrOrig = handles.hTextThrOrigBasalMC;
    else
        post_message('Plotting basal respiration pooled scatter plot...');
        handlesToUpdate = [];
    end
    [~, hOut] = virt_plot_amplitude_correlation(aggregatedBasalTable, pPlot, ...
                            'GroupingColumn', cleanParamName, ...
                            'Handles', handlesToUpdate, ...
                            'NCorrelations', nCorrToAnalyze, ...
                            'FigTitle', figTitleScatterBasal, ...
                            'FigName', figNameScatterBasal, ...
                            'OutDir', dirOut, 'ToSaveOutput', saveOutput, ...
                            'FigTypes', figTypes, 'ShowFigure', showFigure);
    if isfield(hOut, 'axScatter')
        handles.figScatterBasalMC = hOut.fig;
        handles.axScatterBasalMC = hOut.axScatter;
        handles.hScattersBasalMC = hOut.hScatters;
        handles.hBestFitsBasalMC = hOut.hBestFits;
        handles.hTextBestFitBasalMC = hOut.hTextBestFit;
        handles.hThrOrigsBasalMC = hOut.hThrOrigs;
        handles.hTextThrOrigBasalMC = hOut.hTextThrOrig;
    end

    % Decide whether to create or update the Z-score jitter plot
    if ~alwaysNew && isfield(handles, 'figJitterBasalZScore') && ...
            isgraphics(handles.figJitterBasalZScore)
        post_message('Updating basal Z-score jitter plot...');
        handlesToUpdate.fig = handles.figJitterBasalZScore;
        handlesToUpdate.axJitter = handles.axJitterBasalZScore;
        handlesToUpdate.hJitter = handles.hJitterBasalZScore;
        handlesToUpdate.hErrorBars = handles.hErrorBarsBasalZScore;
        handlesToUpdate.hMeans = handles.hMeansBasalZScore;
        handlesToUpdate.pTextJitter = handles.pTextJitterBasalZScore;
        handlesToUpdate.hNull = handles.hNullBasalZScore;
        handlesToUpdate.sigMarkerJitter = handles.sigMarkerJitterBasalZScore;
        handlesToUpdate.transformedLabels = handles.transformedLabelsBasalZScore;
    else
        post_message('Plotting basal Z-score jitter plot...');
        handlesToUpdate = [];
    end
    [~, hOut] = virt_plot_jitter(aggregatedWhiskVectorTable, pPlot, ...
                    'DataColumn', 'whiskAmpCorrZScoresBasal', ...
                    'DataMode', 'FisherZScore', ...
                    'MaxOrders', nCorrToAnalyze, ...
                    'GroupingColumn', cleanParamName, ...
                    'Handles', handlesToUpdate, ...
                    'FigTitle', figTitleJitterBasalAmpCorrZScore, ...
                    'FigName', figNameJitterBasalAmpCorrZScore, ...
                    'OutDir', dirOut, 'ToSaveOutput', saveOutput, ...
                    'FigTypes', figTypes, 'ShowFigure', showFigure);
    if isfield(hOut, 'axJitter')
        handles.figJitterBasalZScore = hOut.fig;
        handles.axJitterBasalZScore = hOut.axJitter;
        handles.hJitterBasalZScore = hOut.hJitter;
        handles.hErrorBarsBasalZScore = hOut.hErrorBars;
        handles.hMeansBasalZScore = hOut.hMeans;
        handles.hNullBasalZScore = hOut.hNull;
        handles.pTextJitterBasalZScore = hOut.pTextJitter;
        handles.sigMarkerJitterBasalZScore = hOut.sigMarkerJitter;
        handles.transformedLabelsBasalZScore = hOut.transformedLabels;
    end
end

% Plot Phase Response Curve for aggregated data if requested
if plotPRC
    % Decide whether to create or update the PRC plot
    if ~alwaysNew && isfield(handles, 'figPrcMC') && isgraphics(handles.figPrcMC)
        post_message('Updating phase response curve...');
        handlesToUpdate.fig = handles.figPrcMC;
        handlesToUpdate.axPRC = handles.axPrcMC;
        handlesToUpdate.hPRCScatter = handles.hPRCScatterMC;
        handlesToUpdate.hPRCRegLine = handles.hPRCRegLineMC;
        handlesToUpdate.hPRCRegText = handles.hPRCRegTextMC;
        handlesToUpdate.hBinError = handles.hPRCBinErrorMC;
        handlesToUpdate.hBinSig = handles.hPRCBinSigMC;
    else
        post_message('Plotting phase response curve...');
        handlesToUpdate = [];
    end
    [~, hOut] = virt_plot_phase_response(aggregatedBasalTable, pPlot, ...
                                    'GroupingColumn', cleanParamName, ...
                                    'Handles', handlesToUpdate, ...
                                    'WhiskDir', whiskDirForPhase, ...
                                    'FigTitle', figTitlePRC, 'FigName', figNamePRC, ...
                                    'OutDir', dirOut, 'ToSaveOutput', saveOutput, ...
                                    'FigTypes', figTypes, 'ShowFigure', showFigure);
    if isfield(hOut, 'axPRC')
        handles.figPrcMC = hOut.fig;
        handles.axPrcMC = hOut.axPRC;
        handles.hPRCScatterMC = hOut.hPRCScatter;
        handles.hPRCRegLineMC = hOut.hPRCRegLine;
        handles.hPRCRegTextMC = hOut.hPRCRegText;
        handles.hPRCBinErrorMC = hOut.hBinError;
        handles.hPRCBinSigMC = hOut.hBinSig;
    end
end

% Plot aggregated scalar outputs
% Decide whether to create or update the scalar plot
if ~alwaysNew && isfield(handles, 'scalarPlot') && ...
        isfield(handles.scalarPlot, 'fig') && isgraphics(handles.scalarPlot.fig)
    post_message('Updating aggregated scalar plot...');
    handlesToUpdate = handles.hScalarPlot;
    updateOrCreate = 'updated';
else
    post_message('Plotting aggregated scalar plot...');
    handlesToUpdate = [];
    updateOrCreate = 'created';
end
hOut = plot_scalar_outputs(aggregatedWhiskScalarTable, ...
                            varsToPlot, varLabels, ...
                            cleanParamName, rangeToVary, ...
                            dirOut, saveOutput, figTitleWhiskScalar, ...
                            figNameWhiskScalarSummary, figTypes, ...
                            'Handles', handlesToUpdate, 'AlwaysNew', alwaysNew, ...
                            'ShowFigure', showFigure);
if isfield(hOut, 'ax')
    handles.hScalarPlot = hOut.fig;
    handles.hScalarPlotAx = hOut.ax;
    handles.hScalarPlotSupAx = hOut.supAx;
    handles.hScalarPlotCurves = hOut.hTuning;
end
post_message(['Aggregated scalar plot ', updateOrCreate, '.']);

%% FOR GUI
% --- Signal Completion for waitfor ---
% Set a property on all generated figure handles to signal that all plotting
% is complete. The GUI can use waitfor() on this property.
if ~isempty(handles)
    % Flush the graphics queue and then pause briefly to yield
    % execution to the renderer.
    if showFigure
        drawnow;
        pause(0.1);
    end

    % Set figure handles as ready
    allHandles = struct2cell(handles);
    for i = 1:numel(allHandles)
        currentHandleArray = allHandles{i};
        % Handle cases where a field contains an array of handles (like axes)
        for j = 1:numel(currentHandleArray)
            currentHandle = currentHandleArray(j);
            if isgraphics(currentHandle) && strcmp(get(currentHandle, 'Type'), 'figure')
                % Set handle as ready
                set(currentHandle, 'UserData', 'Ready');
            end
        end
    end
end

post_message('Multiple repetition analysis and plotting complete!');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [handlesRep, pRep, stateRep, analysisRep, analysisStringRep, plotDataRep] = ...
    run_single_repetition(iRep, pBase, seedsToUse, paramToVary, rangeToVary, ...
        cleanParamName, testName, dirOut, handlesPrev, messageHandler, saveOutput, ...
        computeWhiskAmp, computeWhiskInterval, figTypes, plotSettings, showFigure)
%% Runs a single repetition of the simulation
% This helper function encapsulates the logic for setting up parameters,
% directory names, and running the core simulation for a single step.

    % Set the seed for the current repetition
    currentSeed = seedsToUse(iRep);

    % Vary the specified parameter
    paramToVaryValue = rangeToVary(iRep);
    pCurrentRep = setfield_recursive(pBase, paramToVary, paramToVaryValue);

    % Update linked parameters if a period or variability parameter is being tuned
    switch paramToVary
        case '.SquareInput.Pb.period1'
            pCurrentRep.SquareInput.Pb.perVar1 = (pCurrentRep.SquareInput.Pb.perVarPerc1 / 100) * paramToVaryValue;
        case '.SquareInput.Pb.perVar1'
            pCurrentRep.SquareInput.Pb.perVarPerc1 = (paramToVaryValue / pCurrentRep.SquareInput.Pb.period1) * 100;
        case '.SquareInput.Pb.perVarPerc1'
            pCurrentRep.SquareInput.Pb.perVar1 = (paramToVaryValue / 100) * pCurrentRep.SquareInput.Pb.period1;
        case '.SquareInput.Pb.period2'
            pCurrentRep.SquareInput.Pb.perVar2 = (pCurrentRep.SquareInput.Pb.perVarPerc2 / 100) * paramToVaryValue;
        case '.SquareInput.Pb.perVar2'
            pCurrentRep.SquareInput.Pb.perVarPerc2 = (paramToVaryValue / pCurrentRep.SquareInput.Pb.period2) * 100;
        case '.SquareInput.Pb.perVarPerc2'
            pCurrentRep.SquareInput.Pb.perVar2 = (paramToVaryValue / 100) * pCurrentRep.SquareInput.Pb.period2;
    end

    % Create string for the parameter varied
    strParamVaried = sprintf('%s-%g', cleanParamName, paramToVaryValue);

    % Define a unique output directory for the current repetition
    repOutDir = fullfile(dirOut, strParamVaried);

    % Define a unique test name for the current repetition
    repTestName = virt_updated_test_name_for_rep(testName, cleanParamName, strParamVaried);

    % Run simulation
    % Note: pCurrentRep already contains the correct RngAlgorithm from pBase
    [handlesRep, pRep, stateRep, analysisRep, analysisStringRep, plotDataRep] = ...
        virt_moore('Handles', handlesPrev, 'Params', pCurrentRep, ...
                   'TestName', repTestName, 'OutDir', repOutDir, ...
                   'MessageHandler', messageHandler, 'SaveOutput', saveOutput, ...
                   'RandomizeSeed', false, 'SeedNumber', currentSeed, ...
                   'PlotSampleTraces', plotSettings.PlotSampleTraces, ...
                   'PlotRaster', plotSettings.PlotRaster, ...
                   'ComputeWhiskAmp', computeWhiskAmp, ...
                   'ComputeWhiskInterval', computeWhiskInterval, ...
                   'PlotJitter', plotSettings.PlotJitter, ...
                   'PlotScatter', plotSettings.PlotScatter, ...
                   'PlotPhaseResponse', plotSettings.PlotPhaseResponse, ...
                   'FigTypes', figTypes, ...
                   'ShowFigure', showFigure);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cleanParamName = clean_param_name(paramPath)
%% Cleans up a parameter name
% Similar to logic from virt_moore_gui.mlapp but does NOT have leading '-'

% Split the path into parts
parts = strsplit(paramPath, '.');

% Get the param name to place in test name field
%   Since the first part is '' this will start with '-'
cleanParamName = strjoin(parts(1:end), '-');

% Abbreviate common parameter names for a shorter test name
cleanParamName = replace(cleanParamName, '-Analysis', '');
cleanParamName = replace(cleanParamName, '-Plotting', '');
cleanParamName = replace(cleanParamName, '-Connect', '');
cleanParamName = replace(cleanParamName, '-AdEx', '');
cleanParamName = replace(cleanParamName, '-Synaptic', '');
cleanParamName = replace(cleanParamName, '-ISynWeight', '-wt');
cleanParamName = replace(cleanParamName, '-Variability', '-Var');
cleanParamName = replace(cleanParamName, '-ExtCurrent', '-IExt');
cleanParamName = replace(cleanParamName, '-SquareInput-Pb', '-SqInpt');

% Remove leading '-'
cleanParamName = regexprep(cleanParamName, '^-', '');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function s = setfield_recursive(s, fieldPath, val)

% Set nested structure fields using a path string (e.g., 'AdEx.Pb.Cm')
parts = strsplit(fieldPath, '.');
s = setfield(s, parts{2:end}, val); % Start from index 2 to skip initial '.'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function local_post_message(msg, handler)
% Posts a message to the command window and optionally to a GUI handler.

% Always print to command window
fprintf('%s\n', msg);

% Pass to GUI message board if handle exists
if ~isempty(handler) && isa(handler, 'function_handle')
    handler(msg); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function scalarTable = aggregate_scalar_outputs(analysisStructArray, paramToVary, rangeToVary, seedsUsed)
%% Aggregates scalar numeric fields from the .whisk sub-structure across all repetitions.
%
% This function iterates through the analysis results from each simulation
% run, identifies all scalar numeric fields in the 'whisk' analysis
% sub-structure, and compiles them into a single summary table.
%
% Arguments:
%       analysisStructArray - A structure array of analysis results for each run.
%
% Outputs:
%       scalarTable         - A table where each row corresponds to a repetition
%                             and columns correspond to the scalar outputs.

% Get the number of repetitions
nRepetitions = numel(analysisStructArray);

% Return an empty table if there are no analysis results to process
if nRepetitions == 0 || ~isfield(analysisStructArray(1), 'whisk')
    scalarTable = table();
    return;
end

% Get the 'whisk' sub-structure from the first repetition to identify fields
firstWhiskStruct = analysisStructArray(1).whisk;
allFieldNames = fieldnames(firstWhiskStruct);

% Identify which fields are scalar and numeric
isScalarNumeric = cellfun(@(f) isnumeric(firstWhiskStruct.(f)) && ...
                                isscalar(firstWhiskStruct.(f)), ...
                                allFieldNames);
scalarFieldNames = allFieldNames(isScalarNumeric);

% Get the number of scalar fields found
nScalarFields = numel(scalarFieldNames);

% Pre-allocate a cell array to hold all the data for efficiency
data = cell(nRepetitions, nScalarFields);

% Loop through each repetition to extract the scalar data
for iRep = 1:nRepetitions
    for iField = 1:nScalarFields
        fieldName = scalarFieldNames{iField};
        % Extract and store the value, using NaN as a placeholder if a field is missing
        if isfield(analysisStructArray(iRep).whisk, fieldName)
            data{iRep, iField} = analysisStructArray(iRep).whisk.(fieldName);
        else
            data{iRep, iField} = NaN;
        end
    end
end

% Create the final table from the collected data
scalarTable = cell2table(data, 'VariableNames', scalarFieldNames);

% Add seed number if it doesn't already exist
if ~ismember('seedNumber', scalarFieldNames)
    scalarTable = addvars(scalarTable, seedsUsed, 'Before', 1, ...
                            'NewVariableNames', 'seedNumber');
end

% Add parameter varied if it doesn't already exist
if ~strcmp(paramToVary, 'seedNumber') && ~ismember(paramToVary, scalarFieldNames)
    scalarTable = addvars(scalarTable, rangeToVary(:), 'Before', 1, ...
                            'NewVariableNames', paramToVary);
end

% Add repetition number if it doesn't already exist
if ~ismember('repetitionNumber', scalarFieldNames)
    scalarTable = addvars(scalarTable, (1:nRepetitions)', 'Before', 1, ...
                            'NewVariableNames', 'repetitionNumber');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vectorTable = aggregate_vector_outputs(analysisStructArray, paramToVary, rangeToVary, seedsUsed)
%% Aggregates numeric vector fields from the .whisk sub-structure across all repetitions.
%
% This function iterates through the analysis results from each simulation
% run, identifies all numeric vector fields (non-scalar) in the 'whisk'
% analysis sub-structure, and compiles them into a single summary table.
%
% Arguments:
%       analysisStructArray - A structure array of analysis results for each run.
%
% Outputs:
%       vectorTable         - A table where each row corresponds to a repetition
%                             and columns correspond to the vector outputs.

% Get the number of repetitions
nRepetitions = numel(analysisStructArray);

% Return an empty table if there are no analysis results to process
if nRepetitions == 0 || ~isfield(analysisStructArray(1), 'whisk')
    vectorTable = table();
    return;
end

% Get the 'whisk' sub-structure from the first repetition to identify fields
firstWhiskStruct = analysisStructArray(1).whisk;
allFieldNames = fieldnames(firstWhiskStruct);

% Identify which fields are numeric vectors (and not scalars)
isVectorNumeric = cellfun(@(f) isnumeric(firstWhiskStruct.(f)) && ...
                               isvector(firstWhiskStruct.(f)) && ...
                               ~isscalar(firstWhiskStruct.(f)), ...
                               allFieldNames);
vectorFieldNames = allFieldNames(isVectorNumeric);

% Get the number of vector fields found
nVectorFields = numel(vectorFieldNames);

% Pre-allocate a cell array to hold all the data for efficiency
data = cell(nRepetitions, nVectorFields);

% Loop through each repetition to extract the vector data
for iRep = 1:nRepetitions
    for iField = 1:nVectorFields
        fieldName = vectorFieldNames{iField};
        % Extract and store the value, using NaN as a placeholder if a field is missing
        if isfield(analysisStructArray(iRep).whisk, fieldName)
            % Ensure data is a row vector for consistency in the table
            data{iRep, iField} = {analysisStructArray(iRep).whisk.(fieldName)(:)};
        else
            data{iRep, iField} = {NaN};
        end
    end
end

% Create the final table from the collected data
vectorTable = cell2table(data, 'VariableNames', vectorFieldNames);

% Add seed number if it doesn't already exist
if ~ismember('seedNumber', vectorFieldNames)
    vectorTable = addvars(vectorTable, seedsUsed, 'Before', 1, ...
                            'NewVariableNames', 'seedNumber');
end

% Add parameter varied if it doesn't already exist
if ~strcmp(paramToVary, 'seedNumber') && ~ismember(paramToVary, vectorFieldNames)
    vectorTable = addvars(vectorTable, rangeToVary(:), 'Before', 1, ...
                            'NewVariableNames', paramToVary);
end

% Add repetition number if it doesn't already exist
if ~ismember('repetitionNumber', vectorFieldNames)
    vectorTable = addvars(vectorTable, (1:nRepetitions)', 'Before', 1, ...
                            'NewVariableNames', 'repetitionNumber');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function handles = plot_scalar_outputs(dataTable, varsToPlot, varLabels, paramToVary, rangeToVary, outDir, toSave, figTitle, figBaseName, figTypes, varargin)
%% Plots each scalar output against the repetition number or varied parameter.
%
% This function uses plot_table to create a figure where each subplot shows
% how a single scalar metric varies across the different repetitions. It can
% create a new plot or update an existing one.
%
% Arguments:
%       dataTable      - The table of aggregated scalar outputs.
%       ... (other required arguments)
%       varargin       - 'Handles': A structure of graphics handles to update.
%                      - 'AlwaysNew': A logical flag to force creation of a new figure.
%
% Outputs:
%       handles        - A structure of handles to the plot objects.

%% Hard-coded parameters
plotMode = 'parallel'; 
lineStyle = 'none';
markerType = '.';
markerSize = 6;

%% Default values for optional arguments
handlesDefault = [];
alwaysNewDefault = false;
showFigureDefault = false;

%% Input Parser
iP = inputParser;
addParameter(iP, 'Handles', handlesDefault);
addParameter(iP, 'AlwaysNew', alwaysNewDefault, @islogical);
addParameter(iP, 'ShowFigure', showFigureDefault, @islogical);
parse(iP, varargin{:});
handlesIn = iP.Results.Handles;
alwaysNew = iP.Results.AlwaysNew;
showFigure = iP.Results.ShowFigure;

%% Preparation
% Return immediately if the input table is empty
if isempty(dataTable) || width(dataTable) < 3
    handles = struct;
    return;
end

% Configure plot based on paramToVary
switch paramToVary
case 'seedNumber'
    rowLabel = 'Seed Number';
    rowValues = dataTable.repetitionNumber;
    rowTickLabels = arrayfun(@num2str, rangeToVary, 'UniformOutput', false);
    rowTickLocs = rowValues;
otherwise
    rowLabel = paramToVary;
    rowValues = rangeToVary(:);
    rowTickLabels = arrayfun(@num2str, rangeToVary, 'UniformOutput', false);
    rowTickLocs = rowValues;
end

% Get the names of the columns to plot
if isempty(varsToPlot)
    varsToPlot = setdiff(dataTable.Properties.VariableNames, {'repetitionNumber', 'seedNumber', paramToVary});
end

% Get the variable labels for plotting
if isempty(varLabels)
    varLabels = varsToPlot;
end

% Set the full figure path (without extension) for saving if requested
if toSave
    figBaseNamePath = fullfile(outDir, figBaseName);
else
    figBaseNamePath = '';
end

%% Plotting
% Check if we are updating an existing plot or creating a new one
if ~alwaysNew && ~isempty(handlesIn) && isfield(handlesIn, 'fig') && isgraphics(handlesIn.fig)
    % --- UPDATE EXISTING PLOT ---
    handles = handlesIn;
    
    % Check if the number of variables to plot matches the existing axes
    nVars = numel(varsToPlot);
    if isfield(handles, 'ax') && numel(handles.ax) == nVars && ...
       isfield(handles, 'hTuning') && numel(handles.hTuning) == nVars
        
        % Update each subplot
        for iVar = 1:nVars
            axCurrent = handles.ax(iVar);
            varName = varsToPlot{iVar};
            yData = dataTable.(varName);
            
            % Update the curve data
            set(handles.hTuning(iVar).curves, 'XData', rowValues, 'YData', yData);
            
            % Update axis labels and limits
            set(get(axCurrent, 'XLabel'), 'String', rowLabel);
            set(get(axCurrent, 'YLabel'), 'String', varLabels{iVar});
            if ~isempty(rowValues)
                axCurrent.XLim = [min(rowValues), max(rowValues)];
            end
            axCurrent.YLimMode = 'auto'; % Let MATLAB recalculate Y limits
            axCurrent.XTick = rowTickLocs;
            axCurrent.XTickLabel = rowTickLabels;
        end
        
        % Update overall figure title
        if isfield(handles, 'supAx') && isgraphics(handles.supAx)
            title(handles.supAx, figTitle);
        end
        
        % Save the figure to file if requested
        if toSave
            save_all_figtypes(handles.fig, figBaseNamePath, figTypes);
        end
    else
        % Fallback: if handles are invalid or don't match, replot
        if isfield(handles, 'fig') && isgraphics(handles.fig)
            delete(handles.fig);
        end
        % Call plot_table to generate a new figure
        handles = plot_table(dataTable, 'PlotMode', plotMode, ...
            'VarsToPlot', varsToPlot, 'VarLabels', varLabels, ...
            'RowValues', rowValues, 'RowLabel', rowLabel, ...
            'RowTickLabels', rowTickLabels, 'RowTickLocs', rowTickLocs, ...
            'ReadoutLabel', 'suppress', 'LineStyle', lineStyle, ...
            'Marker', markerType, 'MarkerSize', markerSize, ...
            'FigTitle', figTitle, 'FigName', figBaseNamePath, ...
            'FigTypes', figTypes, 'AlwaysNew', true, 'ShowFigure', showFigure);
    end
else
    % --- CREATE NEW PLOT ---
    % Call plot_table to generate the figure with specified options
    handles = plot_table(dataTable, ...
                'PlotMode', plotMode, ...
                'VarsToPlot', varsToPlot, ...
                'VarLabels', varLabels, ...
                'RowValues', rowValues, ...
                'RowLabel', rowLabel, ...
                'RowTickLabels', rowTickLabels, ...
                'RowTickLocs', rowTickLocs, ...
                'ReadoutLabel', 'suppress', ...
                'LineStyle', lineStyle, ...
                'Marker', markerType, ...
                'MarkerSize', markerSize, ...
                'FigTitle', figTitle, ...
                'FigName', figBaseNamePath, ...
                'FigTypes', figTypes, ...
                'AlwaysNew', true, 'ShowFigure', showFigure);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%