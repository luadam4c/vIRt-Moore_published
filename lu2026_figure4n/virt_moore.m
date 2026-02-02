function [handles, P, State, analysis, analysisString, plotData, Conn] = virt_moore (varargin)
%% Simulate vIRt model once
% Usage: [handles, P, State, analysis, analysisString, plotData, Conn] = virt_moore (varargin)
% Explanation:
%
% --- Hierarchical AdEx Network Simulation ---
% This script simulates a network of four interconnected neuron pools
% (pre-Botzinger Complex, vIRt-retraction, vIRt-protraction, facial motor) 
% using Adaptive Exponential Integrate-and-Fire (AdEx) neurons.
%
% The network structure and interactions are inspired by brainstem circuits
% involved in rhythmic motor patterns like whisking and its coordination
% with breathing.
%
% All parameters are defined at the beginning of the script, allowing for
% easy modification and exploration.
%
% Example(s):
%       handles = virt_moore;
%       handles = virt_moore('Params', P);
%       handles = virt_moore('Params', P, 'Handles', handles);
%       handles = virt_moore('SaveOutput', false, 'PlotRaster', false);
%
% Outputs:
%       handles     - handles to all created graphics objects.
%                   specified as a structure with fields:
%                       .fig1:  Figure handle for sample traces.
%                       .ax1:   Axes handles for sample trace subplots.
%                       .hPlot: Line handles for voltage/current traces.
%                       .hShade1: Patch handles for input pulse shades.
%                       .hEffectorOverlay: Line handles for the effector overlay.
%                       .fig2:  Figure handle for the raster plot.
%                       .ax2:   Axes handles for raster/effector subplots.
%                       .hRaster: Handles for raster plot markers.
%                       .textPools: Text handles for pool name labels.
%                       .hEffector: Line handle for the effector position trace.
%                       .hShade2: Patch handles for shades on the effector plot.
%                       (additional fields like .fig3, .ax3, etc. are added
%                        for analysis plots if generated).
%       P           - updated parameters structure
%                   specified as a structure
%       State       - simulation results, containing the final state of all 
%                   neurons and full time-series data.
%                   specified as a structure with fields:
%                       .(poolName): Sub-structure for each neuron pool 
%                                    (e.g., State.Pb, State.Rr) containing:
%                           .V:         Membrane potential (mV) of all cells 
%                                       at the final time step.
%                           .w:         Adaptation current (pA) of all cells 
%                                       at the final time step.
%                           .ISynIn:    Synaptic current (pA) of all cells 
%                                       at the final time step.
%                           .stepsBeforeChange: Refractory step counter for 
%                                       all cells.
%                           .spikeTimes: Cell array of spike times (ms) for 
%                                       each neuron in the pool.
%                       .Sample: Sub-structure with sample voltage traces for 
%                                the first neuron of each pool.
%                           .(poolName).Cell1: Full voltage trace (mV).
%                       .X:         Final value of the effector's intermediate 
%                                   variable.
%                       .Y:         Final value of the effector's position.
%                       .effectorPosition: Full time-series of the effector 
%                                   position (a.u.).
%                       .tVecMs:    Time vector for the simulation (ms).
%                       .ExtCurrents: Structure with the full external current 
%                                   time-series for each neuron in each pool.
%       analysis    - results from whisking and phase-response analysis
%                   specified as a structure with fields:
%                       .whisk:     Contains whisk cycle statistics like peak 
%                                   amplitudes, intervals, correlations, etc.
%                       .PRC:       Contains phase response curve data, 
%                                   including phases and IEI ratios.
%       analysisString
%                   - A summary of the analysis results.
%                   specified as a string scalar
%       plotData    - data for plots
%                   specified as a structure
%
% Arguments:
%       varargin    - 'Params': parameters structure from virt_moore_params.m
%                   must be a structure
%                   default == Set by virt_moore_params.m
%                   - 'ParamOld': previous parameters structure
%                   must be a structure
%                   default == []
%                   - 'StateOld': previous state structure
%                   must be a structure
%                   default == []
%                   - 'Handles': handles to graphics objects
%                   specified as a structure with fields:
%                       fig1, fig2, ax1, ax2, etc.
%                   - 'TestName': name of test
%                   must be a string scalar or a character vector
%                   default == 'Jeff081825'
%                   - 'OutDir': output directory
%                   must be a string scalar or a character vector
%                   default == fullfile(pwd, sprintf('Test-%s-%s', timeStamp, testName))
%                   - 'MessageHandler': function handle for displaying messages
%                   must be a function handle
%                   default == [] (print to command window only)
%                   - 'RandomizeSeed': Whether to use a randomize seed
%                   must be a logical scalar
%                   default == false
%                   - 'SeedNumber': A specific seed for the random number
%                   generator. This is only used if 'RandomizeSeed' is false.
%                   If a seed is already present in the 'Params' structure,
%                   a warning will be displayed if they differ.
%                   must be a numeric scalar or empty `[]`
%                   default == [] (set in virt_moore_params.m)
%                   - 'RngAlgorithm': The random number generator algorithm to use.
%                       Options: 'twister', 'simdTwister', 'combRecursive',
%                                'multFibonacci', 'philox', 'threefry',
%                                'v4', 'v5uniform', 'v5normal'
%                   must be a string scalar or a character vector
%                   default == [] (set in virt_moore_params.m)
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
%                   - 'LoadParametersOnly': Whether to only load default
%                   parameters and return without running the simulation.
%                   must be a logical scalar
%                   default == false
%                   - 'FigTypes': Figure type(s) for saving.
%                   must be a string scalar, character vector, or a cell array of character vectors
%                   default == {'png'}
%
% Requires (Please add \Shared\Code\Adams_Functions to path):
%       \Shared\Code\Adams_Functions\archive_dependent_scripts.m
%       \Shared\Code\Adams_Functions\cell2num.m
%       \Shared\Code\Adams_Functions\check_dir.m
%       \Shared\Code\Adams_Functions\compute_trial_numbers.m
%       \Shared\Code\Adams_Functions\convert_colors_to_rgb.m
%       \Shared\Code\Adams_Functions\create_subplots.m
%       \Shared\Code\Adams_Functions\create_time_stamp.m
%       \Shared\Code\Adams_Functions\create_pulse_train_times.m
%       \Shared\Code\Adams_Functions\find_closest.m
%       \Shared\Code\Adams_Functions\isequalstructs.m
%       \Shared\Code\Adams_Functions\plot_raster.m
%       \Shared\Code\Adams_Functions\plot_vertical_shade.m
%       \Shared\Code\Adams_Functions\save_all_figtypes.m
%       \Shared\Code\Adams_Functions\set_figure_properties.m
%       \Shared\Code\vIRt-Moore\virt_analyze_whisk.m
%       \Shared\Code\vIRt-Moore\virt_moore_params.m
%       \Shared\Code\vIRt-Moore\virt_plot_whisk_analysis.m
%
% Used by:
%       \Shared\Code\vIRt-Moore\virt_moore_gui.mlapp
%       \Shared\Code\vIRt-Moore\virt_moore_multiple_reps.m
%       \Shared\Code\vIRt-Moore\virt_run_monte_carlo_simulations.m
%       \Shared\Code\vIRt-Moore\jm_run_virt_sim_replicates.m
%       \Shared\Code\vIRt-Moore\jm_run_virt_sim_vary_ext_current.m

% File History:
% 2025-07-31 Commented and streamlined code, also changed square waves
%               to reflect variable periods for each cycle
% 2025-08-01 Added code for plotting sample traces
% 2025-08-05 Fined tuned parameters
% 2025-08-14 Now uses create_pulse_train_trace.m
% 2025-08-14 Convert into a function that takes a parameters structure as 
%               an argument
% 2025-08-21 Updated to work with GUI and updated default parameters
% 2025-08-21 Fixed bug on not updating Ipulse times
% 2025-08-22 Added toSaveOutput
% 2025-08-22 Converted several flags to optional parameter-value pair arguments
% 2025-08-22 Updated analysis code
% 2025-08-22 Updated analysis to detect and save all peak amplitudes
% 2025-08-22 Updated valley detection to be independent of pulse cycles
% 2025-08-22 Updated interval statistics code
% 2025-08-24 Now plots the effector trace overlayed on sample voltage traces
% 2025-08-26 Fixed plot save order
% 2025-09-03 Updated external current to include across-neuron and
%               across-time variability (low-pass filtered white noise)
% 2025-09-04 Now uses compute_fundamental_frequency.m
% 2025-09-04 Now uses freqfilter.m
% 2025-09-04 Now uses parse_oscillation.m
% 2025-09-05 Added fundFreqRange
% 2025-09-16 Update whisk peak detection to match virt_analyze_sniff_whisk.m
% 2025-09-16 Now runs simulation 10000 ms by default and made plot width larger
% 2025-09-16 Update whisk analysis to include amplitude correlations 
%               and logarithmic decrements per Gemini 
% 2025-09-18 Now plots regression lines
% 2025-09-24 Refactored analysis and plotting into separate functions.
% 2025-09-24 Added phase response analysis and plotting.
% 2025-09-25 Completed code pull out to virt_analyze_whisk.m 
%               and virt_plot_whisk_analysis.m
% 2025-09-25 Added phase detection plotting with horizontal bars.
% 2025-09-30 Allowed simulation to not be rerun if simulation parameters
%               are unchanged
% 2025-09-30 Allowed PB input to have basal respiration to sniffing transitions
% 2025-10-01 Updated default parameters to use Jeff-20250928-SquareInput-period2-160
% 2025-10-02 Now analyzes sniff start windows
% 2025-10-02 Now allows plotting parameters to be changed on the GUI
% 2025-10-03 Now saves the seed in the parameters structure
% 2025-10-03 Added 'LoadParametersOnly' as an optional argument
% 2025-10-05 Added relativeTransitionTime and relativeAnalysisStart parameters
%               and backward compatibility logic.
% 2025-10-05 Added parameters for pre-analysis shade.
% 2025-10-06 Now checks for installation of Statistics and Machine Learning Toolbox, 
%           Signal Processing Toolbox and optionally Parallel Processing Toolbox
% 2025-10-08 Now shifts sample window when transition time is near extremums
% 2025-10-15 Added perVarPerc1 and perVarPerc2 as linked parameters by Gemini
% 2025-10-15 Refactored parameter definition into virt_moore_params.m by Gemini
% 2025-10-16 Moved seed determination logic to virt_moore_params.m by Gemini
% 2025-10-16 Added code to save analysis summary by Gemini
% 2025-10-16 Added 'FigTypes' as an optional argument
% 2025-10-17 Made 'AlwaysNew' an optional argument, default == false, pass to virt_plot_whisk_analysis, and pass in from the GUI
% 2025-10-17 Made 'ShowFigure' an optional argument, default == true, and pass in from the GUI 
% 2026-01-17 Added connectivity matrix as last output
% 2026-01-17 Added 'RngAlgorithm' optional argument

%% Compatibility Check
% Check MATLAB version for function flag support
v = ver('matlab');
versionStr = regexp(v.Version, '^\d+\.\d+', 'match', 'once');
matlabVersion = str2double(versionStr);
supportsRGBAColor = matlabVersion >= 9.7;   % R2019b introduced RGBA Color support for plots

%% Toolbox Dependency Check
% Define required and optional toolboxes
requiredToolboxes = {'Statistics and Machine Learning Toolbox', 'Signal Processing Toolbox', 'Image Processing Toolbox'};
optionalToolboxes = {'Parallel Computing Toolbox'};
missingRequired = {};

% Check for required toolboxes
for i = 1:length(requiredToolboxes)
    toolboxName = requiredToolboxes{i};
    if ~any(any(contains(struct2cell(ver), toolboxName)))
        missingRequired{end+1} = toolboxName;
    end
end

% Throw an error if any required toolboxes are missing
if ~isempty(missingRequired)
    error('virt_moore:MissingToolbox', ...
          'The following required toolbox(es) are not installed:\n%s\nPlease install them to continue.', ...
          strjoin(missingRequired, '\n'));
end

% Check for optional toolboxes and issue a warning if missing
for i = 1:length(optionalToolboxes)
    toolboxName = optionalToolboxes{i};
    if ~any(any(contains(struct2cell(ver), toolboxName)))
        warning('virt_moore:OptionalToolboxMissing', ...
                ['Optional toolbox "%s" not found. ', ...
                 'The script will run, but certain features or performance improvements may be unavailable.'], ...
                toolboxName);
    end
end


%% Check required function
scriptRequired = 'find_closest.m';
folderRequired = fullfile('Shared', 'Code', 'Adams_Functions');
if exist(scriptRequired, 'file') ~= 2
    error('Please add %s to path for required script %s!!\n', folderRequired, scriptRequired);
end

%% Hard-Coded Information For This Script Only
% File names
timeStamp = create_time_stamp;
nameModel = 'virt_moore';
nameSampleFig = 'Sample';
nameRasterFig = 'Raster';
nameParamsFile = 'Params';
nameAnalysisFile = 'Analysis';
nameAnalysisSummaryFile = 'Analysis_Summary';
validRngAlgorithms = {'twister', 'simdTwister', 'combRecursive', ...
                      'multFibonacci', 'philox', 'threefry', ...
                      'v4', 'v5uniform', 'v5normal'};

%% Default values for optional arguments
paramsDefault = [];             % default parameters will be set later
paramOldDefault = struct.empty; % no previous parameters by default
stateOldDefault = struct.empty; % no previous state by default
handlesDefault = [];            % no handles by default
testNameDefault = '';           % set later
dirOutDefault = '';             % default directory name will be set later
messageHandlerDefault = [];
timeLimitsSampleDefault = [];   % default will be calculated based on simDur
toRandomizeSeedDefault = false;
seedNumberDefault = [];
rngAlgorithmDefault = [];       % set later
toComputeWhiskAmpDefault = true;
toComputeWhiskIntervalDefault = true;
toSaveOutputDefault = true;
plotSampleTracesFlagDefault = true;
plotRasterFlagDefault = true;
plotJitterFlagDefault = true;
plotScatterFlagDefault = true;
plotPRCFlagDefault = true;
alwaysNewDefault = false;
showFigureDefault = true;
toCopyScriptDefault = [];       % decide based on whether P is provided
loadParamsOnlyDefault = false;
figTypesDefault = {'png'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Functions

% Function to generate a connectivity matrix
%   Values are either meanISyn * (1 +/- ISynVariability),
%       with probability probConn, or 0 (with probability 1 - probConn)
function connMatrix = generate_conn_matrix(nTarget, nSource, probConn, meanISyn, ISynVariability)
    connMatrix = (rand(nTarget, nSource) < probConn) .* ...
                 (meanISyn + meanISyn * ISynVariability * ...
                                (2*rand(nTarget, nSource) - 1));
end

% Function to post messages
function post_message(msg)
    % Always print to command window
    fprintf('%s\n', msg);
    
    % Also send to GUI if handler exists
    if ~isempty(messageHandler)
        messageHandler(msg);
    end
end

% Function to replace raster data
function hRaster = replace_raster(hRaster, spikeTimesCell)
    % Compute trial numbers
    [trialNos, ~, nTrialsTotal] = compute_trial_numbers(spikeTimesCell);

    % Compute y midpoints
    yMids = cellfun(@(x) nTrialsTotal - x + 1, trialNos, ...
                    'UniformOutput', false);

    % Replace raster data
    hRaster = cellfun(@replace_raster_helper, hRaster, spikeTimesCell, ...
                        yMids, 'UniformOutput', false);

end

% Function to replace raster data
function hPlot = replace_raster_helper(hPlot, spikeTimes, yMids)
    if ~isempty(hPlot) && isvalid(hPlot)
        % Prepare data for a single plot call
        xDots = spikeTimes(:);
        
        % Create a matrix of y-values matching the data matrix size
        yMatrix = repmat(yMids', size(spikeTimes, 1), 1);
        yDots = yMatrix(:);
        
        % Remove any NaNs used for padding
        nanMask = isnan(xDots);
        xDots(nanMask) = [];
        yDots(nanMask) = [];
        
        hPlot.XData = xDots;
        hPlot.YData = yDots;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Params', paramsDefault);
addParameter(iP, 'ParamOld', paramOldDefault);
addParameter(iP, 'StateOld', stateOldDefault);
addParameter(iP, 'Handles', handlesDefault);
addParameter(iP, 'TestName', testNameDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'OutDir', dirOutDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'MessageHandler', messageHandlerDefault);
addParameter(iP, 'TimeLimitsSample', timeLimitsSampleDefault, ...
    @(x) isempty(x) || (isnumeric(x) && isvector(x) && numel(x) == 2 && x(1) < x(2))); % New Parameter
addParameter(iP, 'RandomizeSeed', toRandomizeSeedDefault, ...
    @(x) validateattributes(x, {'logical'}, {'scalar'}));
addParameter(iP, 'SeedNumber', seedNumberDefault, ...
    @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
addParameter(iP, 'RngAlgorithm', rngAlgorithmDefault, ...
    @(x) isempty(x) || any(validatestring(x, validRngAlgorithms)));
addParameter(iP, 'ComputeWhiskAmp', toComputeWhiskAmpDefault, ...
    @(x) validateattributes(x, {'logical'}, {'scalar'}));
addParameter(iP, 'ComputeWhiskInterval', toComputeWhiskIntervalDefault, ...
    @(x) validateattributes(x, {'logical'}, {'scalar'}));
addParameter(iP, 'SaveOutput', toSaveOutputDefault, ...
    @(x) validateattributes(x, {'logical'}, {'scalar'}));
addParameter(iP, 'PlotSampleTraces', plotSampleTracesFlagDefault, ...
    @(x) validateattributes(x, {'logical'}, {'scalar'}));
addParameter(iP, 'PlotRaster', plotRasterFlagDefault, ...
    @(x) validateattributes(x, {'logical'}, {'scalar'}));
addParameter(iP, 'PlotJitter', plotJitterFlagDefault, ...
    @(x) validateattributes(x, {'logical'}, {'scalar'}));
addParameter(iP, 'PlotScatter', plotScatterFlagDefault, ...
    @(x) validateattributes(x, {'logical'}, {'scalar'}));
addParameter(iP, 'PlotPhaseResponse', plotPRCFlagDefault, ...
    @(x) validateattributes(x, {'logical'}, {'scalar'}));
addParameter(iP, 'AlwaysNew', alwaysNewDefault, ...
    @(x) validateattributes(x, {'logical'}, {'scalar'}));
addParameter(iP, 'ShowFigure', showFigureDefault, ...
    @(x) validateattributes(x, {'logical'}, {'scalar'}));
addParameter(iP, 'CopyScript', toCopyScriptDefault, ...
    @(x) isempty(x) || (islogical(x) && isscalar(x)));
addParameter(iP, 'LoadParametersOnly', loadParamsOnlyDefault, ...
    @(x) validateattributes(x, {'logical'}, {'scalar'}));
addParameter(iP, 'FigTypes', figTypesDefault, ...
    @(x) ischar(x) || isstring(x) || iscell(x));

% Read from the Input Parser
parse(iP, varargin{:});
paramsIn = iP.Results.Params;
ParamOld = iP.Results.ParamOld;
StateOld = iP.Results.StateOld;
handles = iP.Results.Handles;
testName = iP.Results.TestName;
dirOut = iP.Results.OutDir;
messageHandler = iP.Results.MessageHandler;
tLimSample = iP.Results.TimeLimitsSample;
toRandomizeSeed = iP.Results.RandomizeSeed;
seedNumberUser = iP.Results.SeedNumber;
rngAlgorithmUser = iP.Results.RngAlgorithm;
toComputeWhiskAmp = iP.Results.ComputeWhiskAmp;
toComputeWhiskInterval = iP.Results.ComputeWhiskInterval;
toSaveOutput = iP.Results.SaveOutput;
plotSampleTracesFlag = iP.Results.PlotSampleTraces;
plotRasterFlag = iP.Results.PlotRaster;
plotJitterFlag = iP.Results.PlotJitter;
plotScatterFlag = iP.Results.PlotScatter;
plotPRCFlag = iP.Results.PlotPhaseResponse;
alwaysNew = iP.Results.AlwaysNew;
showFigure = iP.Results.ShowFigure;
toCopyScript = iP.Results.CopyScript;
loadParamsOnly = iP.Results.LoadParametersOnly;
figTypes = iP.Results.FigTypes;

%% Decide whether to re-run simulation
toRerunSimulation = true;
if ~isempty(StateOld) && isfield(StateOld, 'X') && ...
        ~isempty(ParamOld) && isfield(ParamOld, 'dt')
    % Use isequalstructs to compare simulation parameters, ignoring analysis fields
    if isequalstructs(paramsIn, ParamOld, 'IgnoreFields', {'Analysis', 'Plotting'})
        toRerunSimulation = false;
        post_message('Simulation parameters unchanged. Re-running analysis only...');
    else
        post_message('Simulation parameters changed. Re-running full simulation...');
    end
end

%% Decide on default behaviors
% By default, copy script only if parameters are set by script
if isempty(toCopyScript)
    toCopyScript = isempty(paramsIn);
end

%% I. MODEL PARAMETER DEFINITION SECTION
% =========================================================================

% Create or check the parameters structure. This now handles seed determination.
P = virt_moore_params('ParamsIn', paramsIn, 'RandomizeSeed', toRandomizeSeed, ...
                      'SeedNumber', seedNumberUser, 'RngAlgorithm', rngAlgorithmUser, ...
                      'TestName', testName, 'MessageHandler', messageHandler);

% Update local testName variable in case it was created by the function
testName = P.testName;

% Handle LoadParametersOnly flag
if loadParamsOnly
    % Set all other outputs to empty and return only the parameters
    handles = [];
    State = [];
    analysis = [];
    analysisString = '';
    plotData = [];
    return;
end

% Get pool names and number of pools
poolNames = fieldnames(P.N);
nPools = length(poolNames);

% Decide on output directory name
if isempty(dirOut)
    dirOut = fullfile(pwd, sprintf('Test-%s-%s', timeStamp, testName));
end

% Create output directory if needed
if toSaveOutput
    check_dir(dirOut);
end

% Create a copy of all dependent files
if toCopyScript
    archive_dependent_scripts(mfilename, 'OutFolder', dirOut);
end

% Set the RNG with the determined seed and algorithm
% virt_moore_params already did this if RandomizeSeed was true, 
% but we do it again here to be safe and consistent for all cases.
rng(P.seedNumber, P.rngAlgorithm);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if toRerunSimulation

%% II. INITIALIZATION SECTION
% =========================================================================

% Create time vector
dt = P.dt;
tMax = P.simDur;
tVecMs = (0:dt:tMax)';
nSteps = length(tVecMs);
nPb = P.N.Pb;
nRr = P.N.Rr;
nRp = P.N.Rp;
nFm = P.N.Fm;

% Extract constants for simplicity
maxExpVal = P.maxExpVal;
tauSyn = P.Synaptic.tauSyn;
tauAlpha = P.Effector.tauAlpha;
probPbToRr = P.Connect.probPbToRr;
probRrToRp = P.Connect.probRrToRp;
probRpToRr = P.Connect.probRpToRr;
probRrToFm = P.Connect.probRrToFm;
probRrIntra = P.Connect.probRrIntra;
probRpIntra = P.Connect.probRpIntra;
ISynPbToRr = P.ISynWeight.PbToRr;
ISynRrToRp = P.ISynWeight.RrToRp;
ISynRpToRr = P.ISynWeight.RpToRr;
ISynRrToFm = P.ISynWeight.RrToFm;
ISynRrIntra = P.ISynWeight.RrIntra;
ISynRpIntra = P.ISynWeight.RpIntra;
ISynVar = P.ISynWeight.Variability;
XperSpike = P.Effector.XperSpike;

% Initialize state array for all cells
State = struct();
for iPool = 1:nPools
    pool = poolNames{iPool};
    nCells = P.N.(pool);

    % Membrane potential is initialized to leak reversal potential (EL)
    State.(pool).V = ones(nCells, 1) * P.AdEx.(pool).EL;

    % No adaptive current initially
    State.(pool).w = zeros(nCells, 1);

    % No synaptic input current initially
    State.(pool).ISynIn = zeros(nCells, 1);

    % No cells are in the absolute refractory period
    State.(pool).stepsBeforeChange = zeros(nCells, 1);

    % Initialize array to store spike times
    State.(pool).spikeTimes = cell(nCells, 1);
    for iCell = 1:nCells
        State.(pool).spikeTimes{iCell} = [];
    end
end

% Initialize arrays for sample traces
for iPool = 1:nPools
    pool = poolNames{iPool};
    State.Sample.(pool).Cell1 = zeros(nSteps, 1);
end

% Initialize Effector State
State.X = P.Effector.initialIntermediateVariable;
State.Y = P.Effector.initialEffectorPosition;
State.effectorPosition = zeros(nSteps, 1);
State.effectorPosition(1) = State.Y;

% Generate Connectivity Weight Matrix Between Populations
post_message('Pre-computing Synaptic Connectivity ...');
Conn = struct();
Conn.PbToRr = generate_conn_matrix(nRr, nPb, probPbToRr, ...
                                    ISynPbToRr, ISynVar);
Conn.RrToRp = generate_conn_matrix(nRp, nRr, probRrToRp, ...
                                    ISynRrToRp, ISynVar);
Conn.RpToRr = generate_conn_matrix(nRr, nRp, probRpToRr, ...
                                    ISynRpToRr, ISynVar);
Conn.RrToFm = generate_conn_matrix(nFm, nRr, probRrToFm, ...
                                    ISynRrToFm, ISynVar);

% Generate Connectivity Synaptic Current Strengths Within Populations
if probRrIntra > 0
    % Generate all synaptic current strengths (pA)
    Conn.RrIntra = generate_conn_matrix(nRr, nRr, probRrIntra, ...
                                        ISynRrIntra, ISynVar);

    % Set the diagonal to be 0 (no neurons connects to itself)
    Conn.RrIntra(1:nRr+1:end) = 0;
else
    Conn.RrIntra = zeros(nRr, nRr);
end
if probRpIntra > 0
    % Generate all synaptic current strengths (pA)
    Conn.RpIntra = generate_conn_matrix(nRp, nRp, probRpIntra, ...
                                        ISynRpIntra, ISynVar);

    % Set the diagonal to be 0 (no neurons connects to itself)
    Conn.RpIntra(1:nRp+1:end) = 0;
else
    Conn.RpIntra = zeros(nRp, nRp);
end
post_message('Synaptic Connectivity generation complete.');

% Precompute external current vectors for each neuron
post_message('Pre-computing external currents with noise...');
ExtCurrents = struct();
for iPool = 1:nPools
    pool = poolNames{iPool};
    nCells = P.N.(pool);
    
    % Get parameters for this pool
    params = P.ExtCurrent.(pool);
    meanI = params.mean;
    stdPop = params.stdPop;
    stdTime = params.stdTime;
    tauTime = params.tauTime;

    % 1. Generate mean current for each neuron (across-population variability)
    meanCurrentsPop = meanI + stdPop * randn(1, nCells);

    % 2. Generate time-varying noise component (low-pass filtered white noise)
    noiseVecs = zeros(nSteps, nCells);
    if stdTime > 0
        if tauTime > 0 && isfinite(tauTime)
            % Use the exact update formula (see rationale under Notes and Dextexhe et al 2001)
            %   for an Ornstein-Uhlenbeck process dX = −(1/tau) Xt*dt + sigma * sqrt(2/tau) * dW
            %   where tau is the time constant and 
            %       sigma is the long-term standard deviation of the white noise
            expDecay = exp(-dt / tauTime);
            noiseStdDevThisStep = stdTime * sqrt(1 - exp(-2 * dt / tauTime));
            for iStep = 1:nSteps - 1
                noiseVecs(iStep + 1, :) = noiseVecs(iStep, :) * expDecay + ...
                                           noiseStdDevThisStep * randn(1, nCells);
            end
        else
            % Uncorrelated (white) noise if tauTime is zero or invalid
            noiseVecs = stdTime * randn(nSteps, nCells);
        end
    end

    % 3. Combine to get the final external current matrix for the pool
    ExtCurrents.(pool) = repmat(meanCurrentsPop, nSteps, 1) + noiseVecs;
end
post_message('External current generation complete.');

% Recompute PB input times if needed
post_message('Pre-computing square wave currents ...');
if isfield(P.SquareInput.Pb, 'timesStart') && ~P.SquareInput.Pb.toRecompute
    % Extract square wave times
    timesStart = P.SquareInput.Pb.timesStart;
    timesEnd = P.SquareInput.Pb.timesEnd;
else
    % Check if periods are the same or different
    if P.SquareInput.Pb.period1 == P.SquareInput.Pb.period2 && ...
       P.SquareInput.Pb.perVar1 == P.SquareInput.Pb.perVar2
        % If periods are the same, generate for the full duration
        [timesStart, timesEnd] = ...
            create_pulse_train_times(P.SquareInput.Pb.pulseWidth, ...
                                    P.SquareInput.Pb.period1, tMax, ...
                                    P.SquareInput.Pb.perVar1, 'TimeFirst', 0);
    else
        % If periods are different, generate in two steps

        % Get the transition time
        tTransition = tMax * P.SquareInput.Pb.relativeTransitionTime;

        % Generate start and end times for the first epoch, 
        %   but generate one extra start time for the next epoch
        [timesStart1, timesEnd1] = ...
            create_pulse_train_times(P.SquareInput.Pb.pulseWidth, ...
                                    P.SquareInput.Pb.period1, tTransition, ...
                                    P.SquareInput.Pb.perVar1, 'TimeFirst', 0, ...
                                    'GenerateOneBeyond', true);

        % Second epoch
        % If at least one pulse is created, use the last pulse start
        %   as the previous pulse start. Otherwise start at the transition
        %   time
        if ~isempty(timesStart1)
            [timesStart2, timesEnd2] = ...
                create_pulse_train_times(P.SquareInput.Pb.pulseWidth, ...
                                        P.SquareInput.Pb.period2, tMax, ...
                                        P.SquareInput.Pb.perVar2, ...
                                        'TimePreviousStart', timesStart1(end), ...
                                        'GenerateOneBeyond', false);
        else
            [timesStart2, timesEnd2] = ...
                create_pulse_train_times(P.SquareInput.Pb.pulseWidth, ...
                                        P.SquareInput.Pb.period2, tMax, ...
                                        P.SquareInput.Pb.perVar2, ...
                                        'TimeFirst', tTransition, ...
                                        'GenerateOneBeyond', false);
        end
        
        % Combine results from both epochs
        timesStart = [timesStart1; timesStart2];
        timesEnd = [timesEnd1; timesEnd2];

        % Remove any pulses that are out of limits
        outOfBounds = timesStart < 0 | timesEnd > tMax;
        timesStart = timesStart(~outOfBounds);
        timesEnd = timesEnd(~outOfBounds);
    end
    
    % Save randomized square wave times
    P.SquareInput.Pb.timesStart = timesStart;
    P.SquareInput.Pb.timesEnd = timesEnd;
end

% Generate square wave current vector
low = P.SquareInput.Pb.baseline;
high = P.SquareInput.Pb.baseline + P.SquareInput.Pb.amplitude;
squareWaveIVec = ones(nSteps, 1) .* low;
if ~isempty(timesStart) && ~isempty(timesEnd)
    % Get square wave indices
    iTimesStart = find_closest(tVecMs, timesStart);
    iTimesEnd = find_closest(tVecMs, timesEnd);

    % Add square waves
    for iPulse = 1:length(iTimesStart)
        squareWaveIVec(iTimesStart(iPulse):iTimesEnd(iPulse)) = high;
    end
end
post_message('Square wave current generation complete.');

% Add square wave input for the Pb pool if enabled
%   subtract the original mean amplitude to not double count the amplitude
%   but still maintain inter-cell variability
if P.SquareInput.Pb.enable
    ExtCurrents.Pb = ExtCurrents.Pb + repmat(squareWaveIVec, 1, nPb) - P.ExtCurrent.Pb.mean;
end

% Display Message
post_message('Initialization complete. Starting simulation...');

%% III. MAIN SIMULATION LOOP
% =========================================================================
% Run each time step
for iStep = 1:nSteps-1
    % Get Current Time
    tNow = tVecMs(iStep);
   
    % Initialize Boolean structure containing whether each cell is spiking
    spikesCurrentStep = struct();
    for iPool = 1:nPools
        pool = poolNames{iPool};
        nCells = P.N.(pool);
        spikesCurrentStep.(pool) = false(nCells, 1);
    end
    
    % Integrate states of cells from one pool at a time
    for iPool = 1:nPools
        % Get current pool name, number of cells
        pool = poolNames{iPool}; 
        nCells = P.N.(pool);

        % Extract parameters
        params = P.AdEx.(pool);
        Cm = params.Cm;           % Membrane capacitance (pF)
        gL = params.gL;           % Leak conductance (nS)
        EL = params.EL;           % Reversal potential (mV)
        VT = params.VT;           % Voltage threshold (mV)
        DeltaT = params.DeltaT;   % Slope factor (mV)
        tauW = params.tauW;       % Time constant for adaptation current
        a = params.a;             % Subthreshold adaptive conductance (nS)
        b = params.b;             % Spike triggered adaptive current (pA)
        Vreset = params.Vreset;   % Reset potential after a spike (mV)
        Vspike = params.Vspike;   % Spike detection threshold (mV)
        TRef = params.TRef;       % Absolute refractory period (ms)

        % Get previous states of each cell
        VLast = State.(pool).V;   % Membrane potential (mV)
        wLast = State.(pool).w;   % Adaptation current (pA)
        ISynLast = State.(pool).ISynIn;  % Synaptic input current (pA)
        stepsBeforeChange = State.(pool).stepsBeforeChange; 
                                  % Number of steps before absolute
                                  % refractory period ends
        spikeTimes = State.(pool).spikeTimes;
                                  % Recorded spike times for each cell
        
        % Get pre-computed external current for all cells in this pool
        IExtNowAllCells = ExtCurrents.(pool)(iStep, :)';
               
        % Get whether each cell is in an absolute refractory period
        isNotRefractory = State.(pool).stepsBeforeChange == 0;

        % Initiate changes
        dVdt = zeros(nCells, 1);
        dwdt = zeros(nCells, 1);

        % Extract last states of cells not in the refractory period
        VisNR = VLast(isNotRefractory);
        wIsNR = wLast(isNotRefractory);
        IisNR = ISynLast(isNotRefractory);
        IExtNowIsNR = IExtNowAllCells(isNotRefractory);

        % For each cell not in the refractory period:
        %   (1) Compute exponent value for exponential term
        expVal = (VisNR - VT) / DeltaT;

        %   (2) Compute active current/exponential term:
        %   This is the active current from HH-based sodium channels
        %   modelled to be exponentially related to how close 
        %   the current membrane potential is to the voltage threshold
        Iactive = gL * DeltaT .* exp(min(expVal, maxExpVal));

        %   (3) Compute passive current:
        %   This is the passive leak current from all passive channels
        Ipassive = -gL * (VisNR - EL);

        %   (4) Compute changes in membrane potential
        dVdt(isNotRefractory) = (Ipassive + Iactive - wIsNR + ...
                                    IisNR + IExtNowIsNR) / Cm;

        %   (5) Compute changes in adaptive (subthreshold) currents:
        %   This regresses to 0 exponentially with time constant tauW,
        %       but also changes based on membrane potential
        %   'a' is the sensitivity of the adaptive currents to 
        %       changes in the membrane voltage
        dwdt(isNotRefractory) = (a * (VisNR - EL) - wIsNR) / tauW;
       
        % Update the states for each cell
        VNow = VLast + dt * dVdt;
        wNow = wLast + dt * dwdt;

        % Save states for plotting
        State.Sample.(pool).Cell1(iStep, 1) = VNow(1);

        % Find indices of cells that fired an action potential
        %   (membrane voltage crosses Vspike)
        newlySpiked = find(VNow >= Vspike & isNotRefractory);

        % Update states of cells that have spiked
        VNext = VNow;
        wNext = wNow;
        if ~isempty(newlySpiked)
            % Reset membrane potential of next step to reset potential
            VNext(newlySpiked) = Vreset;

            % Add on spike-triggered adaptive current to next step
            wNext(newlySpiked) = wNext(newlySpiked) + b;

            % Update boolean structure
            spikesCurrentStep.(pool)(newlySpiked) = true;

            % Update absolute refractory period:
            %   This is the number of steps before voltage next integrated
            %   Before then, V = Vreset
            stepsBeforeChange(newlySpiked) = round(TRef / dt);

            % Record spike times
            for iSpike = 1:length(newlySpiked)
                spikeTimes{newlySpiked(iSpike)}(end + 1, 1) = tNow;
            end
        end

        % Decrement the non-zero absolute refractory periods
        stepsBeforeChange(stepsBeforeChange > 0) = ...
            stepsBeforeChange(stepsBeforeChange > 0) - 1;

        % Update synaptic input current based on decay:
        %   This decays exponentially with time constant tauSyn
        ISynNext = ISynLast * exp(-dt / tauSyn);

        % Update the states for each cell in the pool
        State.(pool).V = VNext;
        State.(pool).w = wNext;
        State.(pool).ISynIn = ISynNext;
        State.(pool).stepsBeforeChange = stepsBeforeChange; 
        State.(pool).spikeTimes = spikeTimes;
    end

    % Update synaptic currents based on new spikes
    if any(spikesCurrentStep.Pb)
        State.Rr.ISynIn = State.Rr.ISynIn + Conn.PbToRr * spikesCurrentStep.Pb; 
    end
    if any(spikesCurrentStep.Rp)
        State.Rr.ISynIn = State.Rr.ISynIn + Conn.RpToRr * spikesCurrentStep.Rp;
    end
    if probRrIntra > 0 && any(spikesCurrentStep.Rr)
        State.Rr.ISynIn = State.Rr.ISynIn + Conn.RrIntra * spikesCurrentStep.Rr;
    end
    if any(spikesCurrentStep.Rr)
        State.Rp.ISynIn = State.Rp.ISynIn + Conn.RrToRp * spikesCurrentStep.Rr;
    end
    if probRpIntra > 0 && any(spikesCurrentStep.Rp)
        State.Rp.ISynIn = State.Rp.ISynIn + Conn.RpIntra * spikesCurrentStep.Rp;
    end
    if any(spikesCurrentStep.Rr)
        State.Fm.ISynIn = State.Fm.ISynIn + Conn.RrToFm * spikesCurrentStep.Rr; 
    end
    
    % --- Update Effector Position (Alpha Function) ---
    if isfield(P, 'Effector')
        % Extract X variable
        X = State.X;
        Y = State.Y;
        nFmSpikesThisStep = sum(spikesCurrentStep.Fm);
        
        % Increment intermediate variable X by spike contributions
        % Each spike adds 'XperSpike' to X
        X = X + XperSpike * nFmSpikesThisStep;
        
        % ODE for X (intermdiate variable): dX/dt = -X / tauAlpha
        % Solved by Euler method
        dXdt = -X / tauAlpha;
        X = X + dXdt * dt;
        
        % ODE for Y (effector_position): dY/dt = (X - Y) / tauAlpha
        % Solved by Euler method
        dYdt = (X - Y) / tauAlpha;
        Y = Y + dYdt * dt;
        
        % Store trace
        State.X = X;
        State.Y = Y;
        if iStep < nSteps 
            % iStep goes up to nSteps-1, so iStep + 1 goes up to nSteps
            State.effectorPosition(iStep + 1) = Y;
        end
    end
    % --- End Effector Update ---

    % Display progress from time to time
    if mod(iStep, round(nSteps/20)) == 0 || iStep == 1
        post_message(sprintf('Simulation progress: %.0f%% (Time: %.1f ms)', ...
                (iStep/(nSteps-1))*100, tNow));
    end
end

% Display Message
post_message('Simulation finished!');

% Store variables needed for analysis/plotting that are not saved in P
State.tVecMs = tVecMs;
State.ExtCurrents = ExtCurrents;

else

% If simulation is not re-run, the state structure stays the same
State = StateOld;

end

% Save parameters if requested
if toSaveOutput
    save(fullfile(dirOut, sprintf('%s_%s_%s.mat', nameModel, nameParamsFile, timeStamp)), 'P');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% IV. PLOTTING SECTION
% =========================================================================

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

% Extract variables needed for analysis/plotting
% This ensures they exist even if the simulation was not re-run
tMax = P.simDur;
tVecMs = State.tVecMs;
ExtCurrents = State.ExtCurrents;

% Convert all color strings to RGB vectors within P.Plotting and save as pPlot
pPlot = convert_colors_to_rgb(P.Plotting, 'ColorSubStr', 'color');

% Extract plotting parameters for local use
ILimSample = pPlot.ILimSample;
vLimSample = pPlot.vLimSample;
eLimSample = pPlot.eLimSample;
colorPb = pPlot.colorPb;
colorRr = pPlot.colorRr;
colorRp = pPlot.colorRp;
colorFm = pPlot.colorFm;
colorEffectorOverlay = pPlot.colorEffectorOverlay;
transparencyOverlay = pPlot.transparencyOverlay;
markerSizeRaster = pPlot.markerSizeRaster;
lineWidthSampleTraces = pPlot.lineWidthSampleTraces;
lineWidthOverlay = pPlot.lineWidthOverlay;
sampleDuration = pPlot.sampleDuration;

% Extract transition time
relativeTransitionTime = P.SquareInput.Pb.relativeTransitionTime;

% Decide on default time limits for the sample plot if not provided
if isempty(tLimSample)
    % Get simulation duration in seconds
    simDurSec = P.simDur / 1000;

    % Maximize the amount of data examined
    if simDurSec <= sampleDuration
        % If simulation is shorter than sample duration, 
        %   just use simulation duration
        tLimSample = [0, simDurSec];
    else
        % First try to plot around the transition time
        tLimSample = simDurSec * relativeTransitionTime + ...
                        sampleDuration * [-1/2, 1/2];

        % Shift sample window to maximize the amount of data examined
        if tLimSample(1) < 0
            tLimSample = [0, sampleDuration];
        elseif tLimSample(end) > simDurSec
            tLimSample = [simDurSec - sampleDuration, simDurSec];
        end
    end
end

% Save in pPlot
pPlot.tLimSample = tLimSample;

% Reconstruct the poolColors matrix from individual color settings
poolColors = [colorPb; colorRr; colorRp; colorFm];

% Compute the total number of cells
nCellsTotal = sum(structfun(@(x) x, P.N));

% Create time vector in seconds
tVecSec = tVecMs / 1000;

% Extract Square Wave Current times and convert to seconds
timesStart = P.SquareInput.Pb.timesStart;
timesEnd = P.SquareInput.Pb.timesEnd;
timesInputPb = [timesStart, timesEnd]';
timesInputPbSec = timesInputPb / 1000;   

% Extract Current Trace
IExtPbCell1 = ExtCurrents.Pb(:, 1);

% Extract Voltage Traces
vPbCell1 = State.Sample.Pb.Cell1;
vRrCell1 = State.Sample.Rr.Cell1;
vRpCell1 = State.Sample.Rp.Cell1;
vFmCell1 = State.Sample.Fm.Cell1;

% Extract whisk position vector
effectorVec = State.effectorPosition;

% Extract spiketimes for each pool as a cell array
spikeTimesCell = cell(nPools, 1);
for iPool = 1:nPools
    % Get name of this pool
    pool = poolNames{iPool};

    % Extract spike times for this pool stored as a cell array
    spikeTimesThisPool = State.(pool).spikeTimes;

    % Convert to numeric array, patch extra entries with NaNs
    spikeTimesCell{iPool} = cell2num(spikeTimesThisPool, 'CombineMethod', 'padNaN');
end

% Define color for effector trace, with fallback for older MATLAB
if supportsRGBAColor
    % Use true transparency
    colorOverlay = [colorEffectorOverlay, transparencyOverlay]; 
else
    % Blend with white instead
    colorOverlay = colorEffectorOverlay * transparencyOverlay + ...
                    [1, 1, 1] * (1 - transparencyOverlay); 
end

% Compute whisk position vector for overlaying on sample voltage traces
effectorOverlayI = ILimSample(1) + (effectorVec - eLimSample(1)) * ...
                                    diff(ILimSample) / diff(eLimSample);
effectorOverlayV = vLimSample(1) + (effectorVec - eLimSample(1)) * ...
                                    diff(vLimSample) / diff(eLimSample);

% Package data for external plotting
plotData.tVecSec = tVecSec;
plotData.timesInputPbSec = timesInputPbSec;
plotData.IExtPbCell1 = IExtPbCell1;
plotData.vPbCell1 = vPbCell1;
plotData.vRrCell1 = vRrCell1;
plotData.vRpCell1 = vRpCell1;
plotData.vFmCell1 = vFmCell1;
plotData.effectorVec = effectorVec;
plotData.effectorOverlayI = effectorOverlayI;
plotData.effectorOverlayV = effectorOverlayV;
plotData.spikeTimesCell = spikeTimesCell;

% FIGURE 1: Sample Traces
if plotSampleTracesFlag
    % Display message
    post_message('Plotting sample traces ...');

    % Update or create figure 1
    % Note: There is an unidentified bug when updating cause the 
    %   external figures to be flashing and in the meantime freezes the GUI
    if ~alwaysNew && ~isempty(handles) && isfield(handles, 'fig1') && ...
            isgraphics(handles.fig1) && isfield(handles, 'hPlot')
        % Update plots by swapping the data sources
        for iPlot = 1:5
            handles.hPlot(iPlot).XData = tVecSec';
        end
        handles.hPlot(1).YData = IExtPbCell1';
        handles.hPlot(2).YData = vPbCell1';
        handles.hPlot(3).YData = vRrCell1';
        handles.hPlot(4).YData = vRpCell1';
        handles.hPlot(5).YData = vFmCell1';

        for iPlot = 1:5
            handles.hEffectorOverlay(iPlot).XData = tVecSec';
            if iPlot == 1
                handles.hEffectorOverlay(iPlot).YData = effectorOverlayI';
            else                
                handles.hEffectorOverlay(iPlot).YData = effectorOverlayV';
            end
        end

        % Update shades by deleting old ones and plotting new ones
        cellfun(@delete, handles.hShade1);
        for iShade = 1:numel(handles.hShade1)
            if iShade == 1
                yLimits = ILimSample;
            else
                yLimits = vLimSample;
            end
            handles.hShade1{iShade} = ...
                plot_vertical_shade(timesInputPbSec, yLimits(1), yLimits(2), ...
                                    'AxesHandle', handles.ax1(iShade));
        end

        % Clear previous analysis if any
        delete(findobj(handles.ax1(4), 'Tag', 'analysis_points'));
    else
        % Remove any previous figure 1
        if ~isempty(handles) && isfield(handles, 'fig1') && ...
                isgraphics(handles.fig1)
            delete(handles.fig1);
        end

        % Create subplots
        [fig1, ax1] = create_subplots(5, 1, 'FigExpansion', [0.4, 2], ...
                                    'AlwaysNew', true, 'ShowFigure', showFigure);

        % Initiate handles to plots
        hShade1 = cell(5, 1);
        hPlot = gobjects(5, 1);
        hEffectorOverlay = gobjects(5, 1);

        % Plot external input current to PbCell1
        hold(ax1(1), 'on');
        hShade1{1} = plot_vertical_shade(timesInputPbSec, ILimSample(1), ILimSample(2), ...
                                            'AxesHandle', ax1(1));
        hPlot(1) = plot(ax1(1), tVecSec, IExtPbCell1, 'Color', poolColors(1, :), ...
                            'LineWidth', lineWidthSampleTraces);
        hEffectorOverlay(1) = plot(ax1(1), tVecSec, effectorOverlayI, ...
                                    'Color', colorOverlay, 'LineWidth', lineWidthOverlay);
        hold(ax1(1), 'off');
        xlim(ax1(1), tLimSample); ylim(ax1(1), ILimSample);
        ylabel(ax1(1), 'I_{extToPBotC} (pA)');
        
        % Plot voltage of PbCell1
        hold(ax1(2), 'on');
        hShade1{2} = plot_vertical_shade(timesInputPbSec, vLimSample(1), vLimSample(2), ...
                                            'AxesHandle', ax1(2));
        hPlot(2) = plot(ax1(2), tVecSec, vPbCell1, 'Color', poolColors(1, :), ...
                            'LineWidth', lineWidthSampleTraces);
        hEffectorOverlay(2) = plot(ax1(2), tVecSec, effectorOverlayV, ...
                                    'Color', colorOverlay, 'LineWidth', lineWidthOverlay);
        hold(ax1(2), 'off');
        xlim(ax1(2), tLimSample); ylim(ax1(2), vLimSample);
        ylabel(ax1(2), 'vPBotC (mV)');

        % Plot voltage of RrCell1
        hold(ax1(3), 'on');
        hShade1{3} = plot_vertical_shade(timesInputPbSec, vLimSample(1), vLimSample(2), ...
                                            'AxesHandle', ax1(3));
        hPlot(3) = plot(ax1(3), tVecSec, vRrCell1, 'Color', poolColors(2, :), ...
                            'LineWidth', lineWidthSampleTraces);
        hEffectorOverlay(3) = plot(ax1(3), tVecSec, effectorOverlayV, ...
                                    'Color', colorOverlay, 'LineWidth', lineWidthOverlay);
        hold(ax1(3), 'off');
        xlim(ax1(3), tLimSample); ylim(ax1(3), vLimSample);
        ylabel(ax1(3), 'vIRt_{ret} (mV)');

        % Plot voltage of RpCell1
        hold(ax1(4), 'on');
        hShade1{4} = plot_vertical_shade(timesInputPbSec, vLimSample(1), vLimSample(2), ...
                                            'AxesHandle', ax1(4));
        hPlot(4) = plot(ax1(4), tVecSec, vRpCell1, 'Color', poolColors(3, :), ...
                            'LineWidth', lineWidthSampleTraces);
        hEffectorOverlay(4) = plot(ax1(4), tVecSec, effectorOverlayV, ...
                                    'Color', colorOverlay, 'LineWidth', lineWidthOverlay);
        hold(ax1(4), 'off');
        xlim(ax1(4), tLimSample); ylim(ax1(4), vLimSample);
        ylabel(ax1(4), 'vIRt_{pro} (mV)');

        % Plot voltage of FmCell1
        hold(ax1(5), 'on');
        hShade1{5} = plot_vertical_shade(timesInputPbSec, vLimSample(1), vLimSample(2), ...
                                            'AxesHandle', ax1(5));
        hPlot(5) = plot(ax1(5), tVecSec, vFmCell1, 'Color', poolColors(4, :), ...
                            'LineWidth', lineWidthSampleTraces);
        hEffectorOverlay(5) = plot(ax1(5), tVecSec, effectorOverlayV, ...
                                    'Color', colorOverlay, 'LineWidth', lineWidthOverlay);
        hold(ax1(5), 'off');
        xlim(ax1(5), tLimSample); ylim(ax1(5), vLimSample);
        xlabel(ax1(5), 'Time (s)');
        ylabel(ax1(5), 'vFMN (mV)');
    
        % Link X axes of all subplots
        linkaxes(ax1, 'x'); 

        % Output handles to graphics objects
        handles.fig1 = fig1;
        handles.ax1 = ax1;
        handles.hPlot = hPlot;
        handles.hShade1 = hShade1;
        handles.hEffectorOverlay = hEffectorOverlay;
    end
    
    % Update figure
    if showFigure
        drawnow
        pause(0.1);
    end
    post_message('Example Traces Plot updated!');

end

% FIGURE 2: Raster Plot with Effector Position
% Create figure and axes if needed
if plotRasterFlag  
    % Display message
    post_message('Plotting Raster Plot and Effector Trace ...');

    % Update or create figure 2
    % Note: There is an unidentified bug when updating cause the 
    %   external figures to be flashing and in the meantime freezes the GUI
    % Note: If a pool has no spikes, the raster handle will not exist
    isValidAndExisting = @(h) isempty(h) || isgraphics(h);
    if ~alwaysNew && isfield(handles, 'fig2') && isgraphics(handles.fig2) && ...
       isfield(handles, 'hRaster') && all(cellfun(isValidAndExisting, handles.hRaster))

        % Update x axis limits
        handles.ax2(1).XLim = [0, tMax];
        handles.ax2(2).XLim = [0, tMax];

        % Update pool names location
        for iPool = 1:nPools
            handles.textPools(iPool).Position(1) = tMax * 1.02;
        end

        % Replace raster data
        handles.hRaster = replace_raster(handles.hRaster, spikeTimesCell);

        % Replace effector data
        handles.hEffector.XData = tVecMs';
        handles.hEffector.YData = effectorVec';

        % Update square wave time windows
        delete(handles.hShade2)
        if ~isempty(timesInputPb)
            handles.hShade2 = plot_vertical_shade(timesInputPb, ...
                                            'AxesHandle', handles.ax2(2));
        else
            handles.hShade2 = gobjects;
        end

        % Clear previous analysis if any
        delete(findobj(handles.ax2(2), 'Tag', 'analysis_points'));
        delete(findobj(handles.ax2(2), 'Tag', 'analysis_basal_shades'));
        delete(findobj(handles.ax2(2), 'Tag', 'analysis_sniff_shades'));
    else
        % Delete any previous figure 2
        if ~isempty(handles) && isfield(handles, 'fig2') && ...
                isgraphics(handles.fig2)
            delete(handles.fig2);
        end

        % Clear figure
        fig2 = set_figure_properties('AlwaysNew', true, 'FigExpansion', [2, 1], ...
                                    'ShowFigure', showFigure);
    
        % --- Top Subplot: Spike Raster Plot ---
        % Raster plot takes top 2/3
        ax2(1) = subplot(3, 1, [1 2]); 
        hRaster = plot_raster(spikeTimesCell, 'YLabel', 'Neuron Index', ...
                    'FigTitle', 'Network Activity & Effector Position', ...
                    'ColorMap', poolColors, 'XTickLocs', 'suppress', ...
                    'XLabel', 'suppress', 'PlotMode', 'Dot', ...
                    'MarkerSize', markerSizeRaster, 'AxesHandle', ax2(1));
        hold(ax2(1), 'on');
        rasterOffset = 0;
        textPools = gobjects(nPools, 1);
        for iPool = 1:nPools
            pool = poolNames{iPool};
            nCells = P.N.(pool);
            thisPoolColor = poolColors(iPool, :); 
        
            % Show text of each pool name
            textPools(iPool) = ...
                text(ax2(1), tMax * 1.02, nCellsTotal - (rasterOffset + nCells/2), ...
                    strrep(pool, '_', '-'),...
                'Color', thisPoolColor, 'FontWeight', 'bold');
        
            % Plot a dotted line in between each pool
            % if iPool < nPools
            %     lineYPos = rasterOffset + nCells + 0.5;
            %     plot(ax2(1), get(ax2(1), 'XLim'), [lineYPos lineYPos], 'k:');
            % end
        
            % Increase raster plot offset
            rasterOffset = rasterOffset + nCells;
        end
        hold(ax2(1), 'off');
        grid(ax2(1), 'on');
        
        % --- Bottom Subplot: Effector Position ---
        % Effector position takes bottom 1/3
        ax2(2) = subplot(3, 1, 3);
        hEffector = plot(ax2(2), tVecMs, effectorVec, 'm', 'LineWidth', 1.5);
        xlabel(ax2(2), 'Time (ms)');
        ylabel(ax2(2), 'Effector Position (a.u.)');
        xlim(ax2(2), [0, tMax]);
        grid(ax2(2), 'on');

        % Plot square wave time windows
        if ~isempty(timesInputPb)
            hShade2 = plot_vertical_shade(timesInputPb, 'AxesHandle', ax2(2));
        else
            hShade2 = gobjects;
        end

        % Link X axes of all subplots
        linkaxes(ax2, 'x'); 

        % Output handles to graphics objects
        handles.fig2 = fig2;
        handles.ax2 = ax2;
        handles.hRaster = hRaster;
        handles.textPools = textPools;
        handles.hEffector = hEffector;
        handles.hShade2 = hShade2;
    end
    
    % Update figure
    if showFigure
        drawnow
        pause(0.1);
    end
    post_message('Raster Plot updated!');    
end

%% V. FURTHER ANALYSIS
% =========================================================================

% Set default analysis string and analysis structure as empty
analysisString = '';
analysis = struct.empty;

% Perform whisk analysis
if toComputeWhiskAmp || toComputeWhiskInterval
    % Analyze whisking effector trace
    [analysis, analysisString] = virt_analyze_whisk(P, State, tVecMs, ...
                                    'ToComputeWhiskAmp', toComputeWhiskAmp, ...
                                    'ToComputeWhiskInterval', toComputeWhiskInterval, ...
                                    'ToPostMessage', true, 'MessageHandler', messageHandler);


    % Plot analysis results
    handles = virt_plot_whisk_analysis(analysis, pPlot, ...
                                'Handles', handles, ...
                                'PlotSpikeDetectionFlag', plotSampleTracesFlag && toComputeWhiskInterval, ...
                                'PlotWhiskDetectionFlag', plotRasterFlag && toComputeWhiskAmp, ...
                                'PlotJitterFlag', plotJitterFlag, ...
                                'PlotScatterFlag', plotScatterFlag, ...
                                'PlotPRCFlag', plotPRCFlag, ...
                                'OutDir', dirOut, 'FigTypes', figTypes, ...
                                'TimeStamp', timeStamp, ...
                                'MessageHandler', messageHandler, ...
                                'ToSaveOutput', toSaveOutput, ...
                                'NameModel', nameModel, ...
                                'AlwaysNew', alwaysNew, 'ShowFigure', showFigure);
end

% Save figures after analysis to include analysis points
if plotSampleTracesFlag && toSaveOutput
    post_message('Saving Example Traces Plot ...');
    figPath1 = fullfile(dirOut, sprintf('%s_%s_%s', nameModel, nameSampleFig, timeStamp));
    save_all_figtypes(handles.fig1, figPath1, figTypes);
    post_message('Example Traces Plot saved!');
end
if plotRasterFlag && toSaveOutput
    post_message('Saving Raster Plot ...');
    figPath2 = fullfile(dirOut, sprintf('%s_%s_%s', nameModel, nameRasterFig, timeStamp));
    save_all_figtypes(handles.fig2, figPath2, figTypes);
    post_message('Raster Plot saved!');
end

% TODO: Option to create and save table output as well by replicating 
%       scalars and concatenating contents of cell arrays

% Save analysis results if analysis performed
if ~isempty(analysisString) && toSaveOutput
    % Save analysis structure
    save(fullfile(dirOut, sprintf('%s_%s_%s.mat', nameModel, nameAnalysisFile, timeStamp)), 'analysis');
    
    % Save the analysis string to a text file
    summaryFilePath = fullfile(dirOut, sprintf('%s_%s_%s.txt', nameModel, nameAnalysisSummaryFile, timeStamp));
    try
        fileID = fopen(summaryFilePath, 'w');
        fprintf(fileID, '%s\n', analysisString);
        fclose(fileID);
        post_message(sprintf('Analysis summary saved to %s', summaryFilePath));
    catch ME
        post_message(sprintf('Error saving analysis summary: %s', ME.message));
    end
end

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

post_message('Single simulation analysis and plotting complete!');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%