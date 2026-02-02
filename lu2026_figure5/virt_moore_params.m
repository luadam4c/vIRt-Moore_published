function P = virt_moore_params (varargin)
%% Create or check/update a parameters structure for the vIRt model
% Usage: P = virt_moore_params (varargin)
% Explanation:
%       This function either creates a default parameter structure (P) for
%       the virt_moore.m simulation or checks an existing structure for
%       missing fields to ensure backward compatibility. This centralizes
%       parameter definition and management.
%
% Outputs:
%       P           - A complete and checked parameters structure.
%                   specified as a structure
%
% Arguments:
%       varargin    - 'ParamsIn': An existing parameters structure to check. If empty or
%                       not provided, a default structure is created.
%                   must be a structure with fields:
%                       testName    - Name for the test run
%                       seedNumber  - Seed for the random number generator
%                       rngAlgorithm - Algorithm for the random number generator
%                       simDur      - Total simulation time (ms)
%                       dt          - Time step for Euler integration (ms)
%                       maxExpVal   - Maximum exponent value to prevent Inf/NaN
%                       N           - Structure with neuron counts for each pool:
%                                       Pb, Rr, Rp, Fm
%                       AdEx        - Structure with AdEx neuron parameters for each pool:
%                                     (Fields: Pb, Rr, Rp, Fm)
%                                       Cm:     Membrane capacitance (pF)
%                                       gL:     Leak conductance (nS)
%                                       EL:     Leak reversal potential (mV)
%                                       VT:     Effective spike threshold (mV)
%                                       DeltaT: Slope factor (mV)
%                                       a:      Subthreshold adaptation conductance (nS)
%                                       tauW:   Adaptation time constant (ms)
%                                       b:      Spike-triggered adaptation increment (pA)
%                                       Vreset: Reset potential after spike (mV)
%                                       Vspike: Spike detection threshold (mV)
%                                       TRef:   Absolute refractory period (ms)
%                       Synaptic    - Structure with synaptic parameters:
%                                       tauSyn: Synaptic current time constant (ms)
%                       Connect     - Structure with connection probabilities:
%                                       probPbToRr, probRrToRp, probRpToRr,
%                                       probRrToFm, probRrIntra, probRpIntra
%                       ISynWeight  - Structure with synaptic current weights (pA):
%                                       PbToRr, RrToRp, RpToRr, RrToFm,
%                                       RrIntra, RpIntra, Variability
%                       ExtCurrent  - Structure with external current parameters for each pool:
%                                     (Fields: Pb, Rr, Rp, Fm)
%                                       mean:    Mean external current (pA)
%                                       stdPop:  Std dev across population (pA)
%                                       stdTime: Std dev for time-varying noise (pA)
%                                       tauTime: Time constant for noise filter (ms)
%                       SquareInput - Structure for square wave input to Pb pool:
%                                       enable:      (true/false)
%                                       toRecompute: (true/false) recompute pulse times
%                                       baseline:    Baseline current (pA)
%                                       amplitude:   Pulse amplitude (pA)
%                                       period1:     Period 1 (ms)
%                                       perVar1:     Period variability 1 (ms)
%                                       perVarPerc1: Period variability 1 (%)
%                                       period2:     Period 2 (ms)
%                                       perVar2:     Period variability 2 (ms)
%                                       perVarPerc2: Period variability 2 (%)
%                                       pulseWidth:  Pulse duration (ms)
%                                       relativeTransitionTime: Relative time for period change
%                       Effector    - Structure for effector dynamics:
%                                       initialIntermediateVariable: Initial value for X
%                                       initialEffectorPosition:     Initial effector position
%                                       tauAlpha:    Time constant for alpha function (ms)
%                                       XperSpike:   Contribution of each spike (a.u.)
%                       Analysis    - Structure with whisk analysis parameters:
%                                       relativeAnalysisStart: Relative time to start analysis
%                                       amplitudeDefinition: 'peak-to-equilibrium', 'peak-to-avgvalley', 'peak-to-prevalley', 'peak-to-postvalley'
%                                       fundFreqRange:   [min, max] Hz for fundamental frequency
%                                       fCutoffWhisk:    [low, high] Hz for bandpass filter
%                                       fCutoffRelToFund: Cutoff relative to fundamental
%                                       filterOrderWhisk: Order of Butterworth filter
%                                       promThresholdPercWhisk: Min peak prominence (%)
%                                       minPeakPromWhisk: Min peak prominence (a.u.)
%                                       maxWhiskDurationMs: Max inter-valley interval (ms)
%                                       minPeakDistanceMsWhisk: Min distance between peaks (ms)
%                                       sniffFreqThreshold: Sniffing frequency threshold (Hz)
%                                       basalFreqThreshold: Basal respiration threshold (Hz)
%                                       nWhisksSniffStartToAnalyze: # whisks to analyze
%                                       minWhisksBasalRespToAnalyze: Min # whisks in basal cycle
%                                       maxWhisksBasalRespToAnalyze: Max # whisks in basal cycle
%                                       nCorrToAnalyze: Max # of amplitude correlations
%                                       whiskDirForPhase: 'retraction' or 'protraction'
%                                       breathOnsetLatencyMs: Latency for breath onset (ms)
%                       Plotting    - Structure with plotting parameters:
%                                       sampleDuration, ILimSample, vLimSample, eLimSample
%                                       vMarkerChosenSpikes, colorPb, colorRr, colorRp, colorFm
%                                       markerTypePeak, colorFirstPeakUsed, colorFirstPeakNotUsed
%                                       colorOtherPeak, colorValley, lineStylePeakAmplitude
%                                       colorPeakAmplitude, markerSizeFirstPeak, markerSizeOtherPeak
%                                       markerSizeValley, colorEffectorOverlay, colorBasalCycleShade
%                                       colorSniffStartWin, colorPreAnalysisShade, transparencyOverlay
%                                       transparencyBasalCycle, transparencySniffStartWin,
%                                       transparencyPreAnalysisShade, colorT0, colorT1, colorReset
%                                       barFaceAlpha, markerSizeRaster, markerSizeChosenSpikes
%                                       lineWidthSampleTraces, lineWidthOverlay, lineWidthChosenSpikes
%                                       markerSizeJitter, jitterWidth, colorScatter, markerTypeScatter
%                                       markerSizeScatter, markerLineWidthScatter, colorBestFit
%                                       lineStyleBestFit, lineWidthBestFit, colorThrOrig
%                                       lineStyleThrOrig, lineWidthThrOrig
%                   default == Set in this file
%                   - 'RandomizeSeed': Whether to use a randomize seed
%                   must be a logical scalar
%                   default == false
%                   - 'SeedNumber': A specific seed for the random number
%                       generator, used when creating a new default structure.
%                   must be a numeric scalar
%                   default == [] (User Input > P.seedNumber > 1)
%                   - 'RngAlgorithm': The random number generator algorithm to use.
%                       Options: 'twister', 'simdTwister', 'combRecursive',
%                                'multFibonacci', 'philox', 'threefry',
%                                'v4', 'v5uniform', 'v5normal'
%                   must be a string scalar or a character vector
%                   default == [] (User Input > P.rngAlgorithm > 'twister')
%                   - 'ParamFilePath': The full path to a loaded .mat parameter file.
%                       Used to derive the test name if 'TestName' is not provided.
%                   must be a string scalar or a character vector
%                   default == ''
%                   - 'TestName': A specific test name, used when creating a
%                       new default structure.
%                   must be a string scalar or a character vector
%                   default == ''
%                   - 'MessageHandler': Function handle for displaying messages to a GUI.
%                   must be a function handle
%                   default == []
%                   
% Requires:
%
% Used by:
%       \Shared\Code\vIRt-Moore\jm_automate_virt_sim_vary_ext_current.m
%       \Shared\Code\vIRt-Moore\jm_automate_virt_sim_vary_freq_perVar.m
%       \Shared\Code\vIRt-Moore\jm_automate_virt_sim_vary_freq_reps.m
%       \Shared\Code\vIRt-Moore\virt_moore.m
%       \Shared\Code\vIRt-Moore\virt_moore_multiple_reps.m
%       \Shared\Code\vIRt-Moore\virt_moore_gui.mlapp

% File History:
% 2025-10-15 - Created by Gemini by pulling code from virt_moore.m.
% 2025-10-16 - Moved seed determination logic from virt_moore.m by Gemini.
% 2025-10-16 - Moved test name determination logic from virt_moore_gui.mlapp
% 2025-11-23 - Changed default amplitudeDefinition to 'peak-to-prevalley'
% 2025-11-26 - Changed default amplitudeDefinition back to 'peak-to-avgvalley'
% 2026-01-05 - Changed default parameters to match Jeff-20251121
% 2026-01-14 - Changed default parameters to match Jeff-20260113
%               but simDur-40000-relativeAnalysisStart-0.375-relativeTransitionTime-0.75
% 2026-01-17 - Changed default toRecompute to true for safety
% 2026-01-17 - Save the random number generator algorithm used 
%                   in the parameters structure with default 'twister'

%% Hard-coded parameters
% Default randomization parameters
defaultSeedNumber = 1;              % Default seed is 1 if there is no user input
defaultRngAlgorithm = 'twister';    % Default random number generator algorithm is 'twister'
validRngAlgorithms = {'twister', 'simdTwister', 'combRecursive', ...
                      'multFibonacci', 'philox', 'threefry', ...
                      'v4', 'v5uniform', 'v5normal'};

%% Default analysis parameters
%   Note: Keep this consistent with virt_analyze_sniff_whisk.m
amplitudeDefinition = 'peak-to-avgvalley';        
fundFreqRange = [0.5, 20];          % range of possible fundamental frequencies to be detected
fCutoffWhisk = [3, 25];             % bandpass filter cutoff for whisk trace (Moore et al 2013 used [3, 25] Hz)
fCutoffRelToFund = [];              % don't use this
filterOrderWhisk = 3;               % Butterworth filter order for whisk trace (Moore et al 2013 used 3)
promThresholdPercWhisk = 5;         % Percentage of data range for minimum peak prominence for whisk peaks
minPeakPromWhisk = 5;               % Minimum whisk angle change (degrees) to detect as a peak (Moore et al 2013 used 5 degrees)
maxWhiskDurationMs = 250;           % Maximum whisk inter-valley interval (ms) (Moore et al 2013 used whisk duration < 250 ms)
minPeakDistanceMsWhisk = 30;        % Minimum peak distance (ms) for whisk peaks
sniffFreqThreshold = 4;             % Frequency threshold for sniffing in Hz
basalFreqThreshold = 3;             % Frequency threshold for basal respiration in Hz
nWhisksSniffStartToAnalyze = 5;     % Number of whisks at the start of a sniff period to be analyzed
minWhisksBasalRespToAnalyze = 3;    % Minimum number of whisks at the start of a basal respiration cycle to be analyzed
maxWhisksBasalRespToAnalyze = 7;    % Maximum number of whisks at the start of a basal respiration cycle to be analyzed
nCorrToAnalyze = 4;                 % Maximum number of whisk amplitude correlations to analyze
whiskDirForPhase = 'protraction';   % whisk direction for phase calculations
breathOnsetLatencyMs = 0;           % No latency for simulated data

%% Default plotting parameters
sampleDuration = 2.8;            % Sample duration (sec)
ILimSample = [-50, 550];         % Pb to Rr synaptic current (pA)
vLimSample = [-90, 0];           % Membrane potential (mV)
eLimSample = [0, 250];           % Effector Position Limits
vMarkerChosenSpikes = -10;       % voltage for plotting chosen spikes (mV)

% Colors for Raster Plot Neuron Pools
colorPb = 'red';                 % Neuron pool color for Pb
colorRr = 'blue';                % Neuron pool color for Rr
colorRp = 'darkgreen';           % Neuron pool color for Rp
colorFm = 'black';               % Neuron pool color for Fm

% Colors & Markers for Peak/Valley Detection
markerTypePeak = 'o';
colorFirstPeakUsed = 'blue';
colorFirstPeakNotUsed = 'gray';
colorOtherPeak = 'cyan';
colorValley = 'green';
lineStylePeakAmplitude = '--';
colorPeakAmplitude = 'red';
markerSizeFirstPeak = 5;
markerSizeOtherPeak = 3;
markerSizeValley = 2;

% Colors & Transparency for Shades and Overlays
colorEffectorOverlay = 'orange';
colorBasalCycleShade = 'lightorange';
colorSniffStartWin = 'aquamarine';
colorPreAnalysisShade = 'lightgray';
transparencyOverlay = 0.25;
transparencyBasalCycle = 0.8;
transparencySniffStartWin = 0.8;
transparencyPreAnalysisShade = 0.8;

% Colors & Transparency for Phase Bars
colorT0 = 'darkblue';    % Color for T0 horizontal bar (Dark Blue)
colorT1 = 'cyan';        % Color for T1 horizontal bar (Cyan)
colorReset = 'red';      % Color for reset time bar (Red)
barFaceAlpha = 0.8;      % Transparency for the filled bars

% Marker & Line Sizes
markerSizeRaster = 1;
markerSizeChosenSpikes = 3;
lineWidthSampleTraces = 1;
lineWidthOverlay = 0.5;
lineWidthChosenSpikes = 1;

% Parameters for Analysis Plots (Jitter, Scatter)
markerSizeJitter = 8;
jitterWidth = 0.3;
colorScatter = 'blue';
markerTypeScatter = 'o';
markerSizeScatter = 8;
markerLineWidthScatter = 1;
colorBestFit = 'default';       % Use default: if significant 'r', otherwise 'k'
lineStyleBestFit = '--';
lineWidthBestFit = 1;
colorThrOrig = 'gray';
lineStyleThrOrig = '--';
lineWidthThrOrig = 0.5;

%% Default values for optional arguments
pInDefault = [];
toRandomizeSeedDefault = false;
seedNumberDefault = [];
rngAlgorithmDefault = [];
paramFilePathDefault = '';
testNameDefault = '';
messageHandlerDefault = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Functions

% Function to post messages to command window and/or GUI
function post_message(msg)
    % Print to command window
    fprintf('%s\n', msg);

    % Also send to GUI if handler exists
    if ~isempty(messageHandler)
        messageHandler(msg);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ParamsIn', pInDefault);
addParameter(iP, 'RandomizeSeed', toRandomizeSeedDefault, ...
    @(x) validateattributes(x, {'logical'}, {'scalar'}));
addParameter(iP, 'SeedNumber', seedNumberDefault, @isnumeric);
addParameter(iP, 'RngAlgorithm', rngAlgorithmDefault, ...
    @(x) isempty(x) || any(validatestring(x, validRngAlgorithms)));
addParameter(iP, 'ParamFilePath', paramFilePathDefault, @ischar);
addParameter(iP, 'TestName', testNameDefault, @ischar);
addParameter(iP, 'MessageHandler', messageHandlerDefault);

% Read from the Input Parser
parse(iP, varargin{:});
paramsIn = iP.Results.ParamsIn;
toRandomizeSeed = iP.Results.RandomizeSeed;
seedNumberUser = iP.Results.SeedNumber;
rngAlgorithmUser = iP.Results.RngAlgorithm;
paramFilePath = iP.Results.ParamFilePath;
testName = iP.Results.TestName;
messageHandler = iP.Results.MessageHandler;

%% Determine Random Number Generator Algorithm
% Priority: 
% 1. User Input ('RngAlgorithm')
% 2. Existing Parameters Structure (paramsIn.rngAlgorithm)
% 3. Default ('twister')

% Check if an algorithm exists in the passed parameters structure
rngAlgorithmInParams = '';
if ~isempty(paramsIn) && isfield(paramsIn, 'rngAlgorithm')
    rngAlgorithmInParams = paramsIn.rngAlgorithm;
end

if ~isempty(rngAlgorithmUser)
    % User provided an algorithm
    rngAlgorithmToUse = rngAlgorithmUser;
    
    if ~isempty(rngAlgorithmInParams) && ~strcmp(rngAlgorithmUser, rngAlgorithmInParams)
        msg = sprintf(['Warning: An RNG algorithm is already present in the parameters structure (%s), ' ...
                       'but a different algorithm was passed as an argument (%s). ' ...
                       'Using the argument value.'], rngAlgorithmInParams, rngAlgorithmUser);
        post_message(msg);
    else
        post_message(sprintf('Using user-provided RNG algorithm: %s', rngAlgorithmToUse));
    end
elseif ~isempty(rngAlgorithmInParams)
    % User did not provide one, but it exists in parameters
    rngAlgorithmToUse = rngAlgorithmInParams;
    post_message(sprintf('Using RNG algorithm from existing parameters structure: %s', rngAlgorithmToUse));
else
    % Neither provided, use default
    rngAlgorithmToUse = defaultRngAlgorithm;
    post_message(sprintf('Using default RNG algorithm: %s', rngAlgorithmToUse));
end

%% Set seed and algorithm
if toRandomizeSeed
    % If randomizing, always generate a new seed and ignore any user input
    % rng('shuffle', generator) sets the generator and seeds it based on current time
    seed = rng('shuffle', rngAlgorithmToUse);

    % Get the seed number from the seed struct
    seedNumberToUse = seed.Seed;

    % Post a message
    post_message(sprintf('Using a new randomized seed: %d with algorithm: %s', ...
        seedNumberToUse, rngAlgorithmToUse));
else
    % If not randomizing, follow the precedence rules
    if ~isempty(seedNumberUser)
        % User passed in a specific seed number
        % Display message if inconsistent
        if ~isempty(paramsIn) && isfield(paramsIn, 'seedNumber') && ...
                ~isempty(paramsIn.seedNumber) && paramsIn.seedNumber ~= seedNumberUser
            msg = sprintf(['Warning: A seed number is already present in the parameters structure (P.seedNumber = %d), ' ...
                           'but a different seed number was passed as an argument (SeedNumber = %d). ' ...
                           'Using the argument value.'], paramsIn.seedNumber, seedNumberUser);
            post_message(msg);
        end

        % Use the user-provided seed number
        seedNumberToUse = seedNumberUser;
        post_message(sprintf('Using user-provided seed number: %d', seedNumberToUse));
    elseif ~isempty(paramsIn) && isfield(paramsIn, 'seedNumber') && ~isempty(paramsIn.seedNumber)
        % User did not pass a seed, but one exists in the parameters structure
        seedNumberToUse = paramsIn.seedNumber;
        post_message(sprintf('Using seed number from existing parameters structure: %d', seedNumberToUse));
    else
        % No seed provided by user or in parameters, so use the default
        seedNumberToUse = defaultSeedNumber;
        post_message(sprintf('Using default seed number: %d', seedNumberToUse));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% I. MODEL PARAMETER DEFINITION / ASSIGNMENT
% =========================================================================
% Modify parameters in this section to explore different network behaviors.
% All currents (w, I_syn, I_ext, b) are in pA.
% All conductances (gL, a) are in nS.
% All capacitances (Cm) are in pF.
% All voltages (V, EL, VT, Vreset, Vspike, DeltaT) are in mV.
% All times (simDur, dt, tauW, TRef, tauSyn) are in ms.
% Membrane Time Constant = Cm / gL
% =========================================================================

if isempty(paramsIn)
    % Initialize
    P = struct();

    % --- 0. Set and save test name and seed ----
    if isempty(testName)
        P.testName = 'Lu-2026-Fig5';
    else
        P.testName = testName;
    end

    % Save the used seed into the parameters structure for record-keeping
    P.seedNumber = seedNumberToUse;
    
    % Save the used algorithm into the parameters structure
    P.rngAlgorithm = rngAlgorithmToUse;

    % --- 1. Simulation Control ---
    P.simDur = 40000;       % ms, Total simulation time
    P.dt = 0.1;             % ms, Time step for Euler integration
    P.maxExpVal = 50;       % maximum exponent value to prevent Inf/NaN
                            % max(expVal) around 50-70 is usually safe.
    
    % --- 2. Neuron Pool Sizes ---
    P.N.Pb = 200;      % Number of neurons in pre-Botzinger Complex (Pb)
    P.N.Rr = 800;      % Number of neurons in vIRt-retraction pool (Rr)
    P.N.Rp = 400;      % Number of neurons in vIRt-protraction pool (Rp)
    P.N.Fm = 100;      % Number of neurons in facial motor pool (Fm)
    
    % --- 3. AdEx Neuron Parameters (for each pool) ---
    
    % Default paramaters
    %  Notes:
    %  (1) Slope factor is usually around 1 mV for cortical pyramidal neurons
    %           However, Naud et al. 2008 used 2 for most neurons
    %  (2) Taken from Chapter 6.2 of Neuronal Dynamics Online Version:
    %       Type        Fig     τm (ms) a (nS)  τw (ms) b (pA)  Vm (mV)
    %       Tonic       6.3A    20      0.0     30.0    60      -55
    %       Adapting    6.3B    20      0.0     100     5.0     -55
    %       Init. burst 6.4A    5.0     0.5     100     7.0     -51
    %       Bursting    6.4C    5.0     -0.5    100     7.0     -46
    %       Irregular   6.5A    9.9     -0.5    100     7.0     -46
    %       Transient   6.9A    10      1.0     100     10      -60
    %       Delayed             5.0     -1.0    100     10      -60
    %       Note: τm = R * Cm = Cm / gL
    %  (3) Taken from Naud et al. 2008:
    %       Type        C (pF)  gL (nS) EL (mV) VT (mV) ΔT (mV) a (nS)  τw (ms) b (pA)  Vr (mV) C(β)    I (pA)
    %       Fig. 8, cNA 59      2.9     −62     −42     3.0     1.8     16      61      −54     8.4     184
    %       Fig. 8, cAD 83      1.7     −59     −56     5.5     2.0     41      55      −54     10.4    116
    %       Fig. 8, RS  104     4.3     −65     −52     0.8     −0.8    88      65      −53     10.4    98
    %       Fig. 4a     200     10      −70     −50     2       2       30      0       −58     −       500
    %       Fig. 4b     200     12      −70     −50     2       2       300     60      −58     −       500
    %       Fig. 4c     130     18      −58     −50     2       4       150     120     −50     −       400
    %       Fig. 4d     200     10      −58     −50     2       2       120     100     −46     −       210
    %       Fig. 4e     200     12      −70     −50     2       −10     300     0       −58     −       300
    %       Fig. 4f     200     12      −70     −50     2       −6      300     0       −58     −       110
    %       Fig. 4g     100     10      −65     −50     2       −10     90      30      −47     −       350
    %       Fig. 4h, 5  100     12      −60     −50     2       −11     130     30      −48     −       160
    %  (4) Examples of spike-triggered adaptive current: Ca-activated K
    %  (5) Golomb used (for both vIRT and FMN cells):
    %           Cm = 1 uF/cm2, implying a surface area of 10^-4 cm2
    %           gL = (0.12 mS/cm2 +/- 0.06 mS/cm2) * 10^-4 cm2 * 10^6 nS/mS = (12 +/- 6) nS
    %           EL = -70 mV
    %           gNa = 100 mS/cm2 * 10^-4 cm2 * 10^6 nS/mS = 10000 nS
    %           gKdr = 20 mS/cm2 * 10^-4 cm2 * 10^6 nS/mS = 2000 nS
    %           gNaP = 0.04 mS/cm2 * 10^-4 cm2 * 10^6 nS/mS = 4 nS
    %           vIRT: gAdapt = 0.7 mS/cm2 * 10^-4 cm2 * 10^6 nS/mS = 70 nS
    %           FMN: gAdapt = 0.3 mS/cm2 * 10^-4 cm2 * 10^6 nS/mS = 30 nS
    %           gh = 0.05 mS/cm2 * 10^-4 cm2 * 10^6 nS/mS = 5 nS
    
    % Nonadaptive neuron default parameters
    default.NonAdap.Cm     = 100;  % pF, Membrane capacitance
    default.NonAdap.gL     = 12;   % nS, Leak conductance
    default.NonAdap.EL     = -70;  % mV, Leak reversal potential
    default.NonAdap.VT     = -50;  % mV, Effective spike initiation threshold
    default.NonAdap.DeltaT = 2;    % mV, Slope factor or Sharpness of action 
                                   %         potential initiation 
    default.NonAdap.a      = 0;    % nS, subthreshold adaptive conductance
    default.NonAdap.tauW   = 300;  % ms, adaptive current time constant
    default.NonAdap.b      = 0;    % pA, spike triggered adaptive current
    default.NonAdap.Vreset = -58;  % mV, Reset potential after a spike
    default.NonAdap.Vspike = 0;    % mV, Spike detection threshold
    default.NonAdap.TRef   = 0;    % ms, Absolute refractory period
    
    % Adaptive neuron default parameters
    default.Adap.Cm      = 100; % pF, Membrane capacitance
    default.Adap.gL      = 12;  % nS, Leak conductance
    default.Adap.EL     = -70;  % mV, Leak reversal potential
    default.Adap.VT     = -50;  % mV, Effective spike initiation threshold
    default.Adap.DeltaT = 2;    % mV, Slope factor or Sharpness of action 
                                %         potential initiation 
    default.Adap.a      = -11;  % nS, subthreshold adaptive conductance
    default.Adap.tauW   = 130;  % ms, adaptive current time constant
    default.Adap.b      = 30;    % pA, spike triggered adaptive current
    default.Adap.Vreset = -48;  % mV, Reset potential after a spike
    default.Adap.Vspike = 0;    % mV, Spike detection threshold
    default.Adap.TRef   = 0;    % ms, Absolute refractory period
    
    % pre-Botzinger Complex, make non-adaptive
    P.AdEx.Pb.Cm      = default.NonAdap.Cm;
    P.AdEx.Pb.gL      = default.NonAdap.gL;
    P.AdEx.Pb.EL      = default.NonAdap.EL; 
    P.AdEx.Pb.VT      = default.NonAdap.VT;
    P.AdEx.Pb.DeltaT  = default.NonAdap.DeltaT; 
    P.AdEx.Pb.a       = default.NonAdap.a;
    P.AdEx.Pb.tauW    = default.NonAdap.tauW;
    P.AdEx.Pb.b       = default.NonAdap.b;
    P.AdEx.Pb.Vreset  = default.NonAdap.Vreset;
    P.AdEx.Pb.Vspike  = default.NonAdap.Vspike;
    P.AdEx.Pb.TRef    = default.NonAdap.TRef;
    
    % vIRt-retraction, make adaptive
    P.AdEx.Rr.Cm      = default.Adap.Cm; 
    P.AdEx.Rr.gL      = default.Adap.gL;
    P.AdEx.Rr.EL      = default.Adap.EL; 
    P.AdEx.Rr.VT      = default.Adap.VT;
    P.AdEx.Rr.DeltaT  = default.Adap.DeltaT; 
    P.AdEx.Rr.a       = default.Adap.a;
    P.AdEx.Rr.tauW    = default.Adap.tauW;
    P.AdEx.Rr.b       = default.Adap.b;
    P.AdEx.Rr.Vreset  = default.Adap.Vreset;
    P.AdEx.Rr.Vspike  = default.Adap.Vspike;
    P.AdEx.Rr.TRef    = default.Adap.TRef;
    
    % vIRt-protraction, make adaptive
    P.AdEx.Rp.Cm      = default.Adap.Cm; 
    P.AdEx.Rp.gL      = default.Adap.gL;
    P.AdEx.Rp.EL      = default.Adap.EL; 
    P.AdEx.Rp.VT      = default.Adap.VT;
    P.AdEx.Rp.DeltaT  = default.Adap.DeltaT;
    P.AdEx.Rp.a       = default.Adap.a;
    P.AdEx.Rp.tauW    = default.Adap.tauW;
    P.AdEx.Rp.b       = default.Adap.b;
    P.AdEx.Rp.Vreset  = default.Adap.Vreset;
    P.AdEx.Rp.Vspike  = default.Adap.Vspike;
    P.AdEx.Rp.TRef    = default.Adap.TRef;
    
    % Facial motor neurons, make non-adaptive
    P.AdEx.Fm.Cm      = default.NonAdap.Cm;
    P.AdEx.Fm.gL      = default.NonAdap.gL;
    P.AdEx.Fm.EL      = default.NonAdap.EL; 
    P.AdEx.Fm.VT      = default.NonAdap.VT;
    P.AdEx.Fm.DeltaT  = default.NonAdap.DeltaT; 
    P.AdEx.Fm.a       = default.NonAdap.a;
    P.AdEx.Fm.tauW    = default.NonAdap.tauW;
    P.AdEx.Fm.b       = default.NonAdap.b;
    P.AdEx.Fm.Vreset  = default.NonAdap.Vreset;
    P.AdEx.Fm.Vspike  = default.NonAdap.Vspike;
    P.AdEx.Fm.TRef    = default.NonAdap.TRef;
    
    % --- 4. Synaptic Parameters ---
    % Note: currently the synaptic current is NOT dependent on membrane
    %           potential
    %   Golomb used:
    %           E_GABAA = -80 mV
    %           tau_GABAA = 10 ms
    %           g_intra = 0.48 mS/cm2
    %           g_PbToRr = 0.5 mS/cm2
    %           g_Adapt = 7 mS/cm2
    %           I_ext_R = 20 uA/cm2
    %           I_ext_F = 3.1 uA/cm2
    %           probConnect = 0.25
    %   Typical cell surface area: 10^-4 cm2
    %           
    P.Synaptic.tauSyn       = 10;     % GABAA time constant (ms)
    P.Connect.probPbToRr    = 0.4;    % Golomb: 1
    P.Connect.probRrToRp    = 0.1;    % Golomb: 0.25
    P.Connect.probRpToRr    = 0.1;    % Golomb: 0.25
    P.Connect.probRrToFm    = 0.1;    % Golomb: 0.25
    P.Connect.probRrIntra   = 0.1;    % Golomb: 0.25, 0 in Fig4
    P.Connect.probRpIntra   = 0.1;    % Golomb: 0.25, 0 in Fig4
    
    % Synaptic currents in between populations (pA)
    % ISyn = meanISyn * (1 +/- ISynVariability)
    P.ISynWeight.PbToRr     = -0.04;  % Golomb: ~1350 pA = 0.5 mS/cm2 * ~27 mV * 10^-4 cm2 * 10^6 pA/uA
    P.ISynWeight.RrToRp     = -3;     % Golomb: ~2160 pA = 0.8 mS/cm2 * ~27 mV * 10^-4 cm2 * 10^6 pA/uA, Fig 4A: 0.24 mS/cm2 = 650 pA
    P.ISynWeight.RpToRr     = -3;     % Golomb: ~2160 pA = 0.8 mS/cm2 * ~27 mV * 10^-4 cm2 * 10^6 pA/uA, Fig 4A: 0.24 mS/cm2 = 650 pA 
    P.ISynWeight.RrToFm     = -3;     % Golomb: ~120 pA = 0.12 mS/cm2 * ~20 mV * 10^-4 cm2 * 10^6 pA/uA
    P.ISynWeight.RrIntra    = -2;     % Golomb: ~1300 pA = 0.48 mS/cm2 * ~27 mV * 10^-4 cm2 * 10^6 pA/uA, Fig 4A: 0
    P.ISynWeight.RpIntra    = -2;     % Golomb: ~1300 pA = 0.48 mS/cm2 * ~27 mV * 10^-4 cm2 * 10^6 pA/uA, Fig 4A: 0
    P.ISynWeight.Variability = 0;
    
    % --- 5. External Currents ---
    % --- 5a. External Currents Amplitudes (in pA) ---
    % 'mean': Mean amplitude of the external current.
    % 'stdPop': Standard deviation of the mean current across the neuron population.
    % 'stdTime': Standard deviation for the time-varying white noise component.
    % 'tauTime': Time constant (ms) for low-pass filtering the time-varying noise.
    %   Mainen & Sejnowski 1995 used convolution of white noise with t*exp(-t/tau) with mu = 150 pA, sigma = 100 pA, tau = 3 ms
    %   Here, we will use the formula for an Ornstein-Uhlenbeck process
    P.ExtCurrent.Pb.mean    = 500;
    P.ExtCurrent.Pb.stdPop  = 0;
    P.ExtCurrent.Pb.stdTime = 0;
    P.ExtCurrent.Pb.tauTime = 0;        % default as GABAA time constant (ms)

    P.ExtCurrent.Rr.mean    = 120;      % Golomb: 20 uA/cm2 * 10^-4 cm2 = 2000 pA
    P.ExtCurrent.Rr.stdPop  = 0;
    P.ExtCurrent.Rr.stdTime = 0;
    P.ExtCurrent.Rr.tauTime = 0;

    P.ExtCurrent.Rp.mean    = 120;      % Golomb: 20 uA/cm2 * 10^-4 cm2 = 2000 pA
    P.ExtCurrent.Rp.stdPop  = 0;
    P.ExtCurrent.Rp.stdTime = 0;
    P.ExtCurrent.Rp.tauTime = 0;

    P.ExtCurrent.Fm.mean    = 280;      % Golomb: 3.1 uA/cm2 * 10^-4 cm2 = 310 pA
    P.ExtCurrent.Fm.stdPop  = 0;
    P.ExtCurrent.Fm.stdTime = 0;
    P.ExtCurrent.Fm.tauTime = 0;
    
    % --- 5b. Square Wave External Current for Pb pool ---
    % Note: Golomb uses a variable frequency/period for each cycle
    %               dur = (700 ms) +/- (150 ms / 2)
    %       and square wave is fixed length 70 ms
    P.SquareInput.Pb.enable = true;         % If false, will just be tonic
    P.SquareInput.Pb.toRecompute = true;    % If true, will recompute time starts and time ends
    P.SquareInput.Pb.baseline  = 0;
    P.SquareInput.Pb.amplitude = P.ExtCurrent.Pb.mean;
    P.SquareInput.Pb.relativeTransitionTime = 0.75;  % Relative time for period change
    P.SquareInput.Pb.period1 = 800;         % Period (ms) for first half of simulation
    P.SquareInput.Pb.perVar1 = 0;           % Period variability (ms) for first half
    if P.SquareInput.Pb.period1 > 0
        P.SquareInput.Pb.perVarPerc1 = (P.SquareInput.Pb.perVar1 / P.SquareInput.Pb.period1) * 100;
    else
        P.SquareInput.Pb.perVarPerc1 = 0;
    end
                                            % Period variability 1 (%)
    P.SquareInput.Pb.period2 = 166;         % Period (ms) for second half of simulation
    P.SquareInput.Pb.perVar2 = 0;           % Period variability (ms) for second half
    if P.SquareInput.Pb.period2 > 0
        P.SquareInput.Pb.perVarPerc2 = (P.SquareInput.Pb.perVar2 / P.SquareInput.Pb.period2) * 100;
    else
        P.SquareInput.Pb.perVarPerc2 = 0;
    end
                                            % Period variability 2 (%)
    P.SquareInput.Pb.pulseWidth = 50;       % Pulse Width (ms), Golomb used 70 ms
    
    % --- 6. Effector Parameters ---
    %     - These parameters control the dynamics of a simulated effector whose
    %       position is driven by the summed activity of facial motor neurons 
    %       using an alpha function response for each spike.
    %     - P.Effector.tauAlpha (ms):
    %           The time constant for the alpha function. This single time constant
    %           governs both the rise and decay phases of the alpha function shape
    %           ($t \cdot \exp(-t/\tau_{\alpha})$). A larger $\tau_{\alpha}$ means the
    %           effector movement is slower to peak and has a longer overall period.
    %           The peak of an individual alpha response occurs at $t = \tau_{\alpha}$.
    %     - P.Effector.XperSpike (arbitrary units):
    %           Scales the impact of each summed facial motor neuron spike. Specifically,
    %           each spike causes an instantaneous increase of this amount in an
    %           intermediate variable $X$. This $X$ then decays exponentially and drives
    %           the actual effector position $E$, resulting in $E$ following an
    %           alpha function shape. Adjust this to get a desired range of
    %           movement for the effector display. The peak height of the alpha
    %           response to a single spike will be $P.Effector.spike\_contribution / e$
    %           (where $e \approx 2.718$ is Euler's number), occurring at $t=\tau_{\alpha}$
    %           after the spike.
    %     - The effector's position $E(t)$ and an intermediate variable $X(t)$ are
    %       updated based on facial motor spikes $S(t)$ (number of spikes in a time step):
    %       1. $X \leftarrow X + P.Effector.spike\_contribution \cdot S(t)$ (instantaneous update due to spikes)
    %       2. $dX/dt = -X / \tau_{\alpha}$
    %       3. $dE/dt = (X - E) / \tau_{\alpha}$
    %
    P.Effector.initialIntermediateVariable = 0; % Initial intermediate variable X at t = 0
    P.Effector.initialEffectorPosition = 0;     % Initial position of Effector at t = 0
    P.Effector.tauAlpha = 10;   % ms, Time constant for alpha function shape of effector response
    P.Effector.XperSpike = 2;   % Arbitrary units, scales the contribution of each motor spike
                                % to the intermediate variable X, which then shapes the alpha response.

    % --- 7. Analysis Parameters ---
    P.Analysis.relativeAnalysisStart = 0.375;
    P.Analysis.amplitudeDefinition = amplitudeDefinition;
    P.Analysis.fundFreqRange = fundFreqRange;
    P.Analysis.fCutoffWhisk = fCutoffWhisk;
    P.Analysis.fCutoffRelToFund = fCutoffRelToFund;
    P.Analysis.filterOrderWhisk = filterOrderWhisk;
    P.Analysis.promThresholdPercWhisk = promThresholdPercWhisk;
    P.Analysis.minPeakPromWhisk = minPeakPromWhisk;
    P.Analysis.maxWhiskDurationMs = maxWhiskDurationMs;
    P.Analysis.minPeakDistanceMsWhisk = minPeakDistanceMsWhisk;
    P.Analysis.sniffFreqThreshold = sniffFreqThreshold;
    P.Analysis.basalFreqThreshold = basalFreqThreshold;
    P.Analysis.nWhisksSniffStartToAnalyze = nWhisksSniffStartToAnalyze;
    P.Analysis.minWhisksBasalRespToAnalyze = minWhisksBasalRespToAnalyze;
    P.Analysis.maxWhisksBasalRespToAnalyze = maxWhisksBasalRespToAnalyze;
    P.Analysis.nCorrToAnalyze = nCorrToAnalyze;
    P.Analysis.whiskDirForPhase = whiskDirForPhase;
    P.Analysis.breathOnsetLatencyMs = breathOnsetLatencyMs;

    % --- 8. Plotting Parameters ---
    P.Plotting.sampleDuration = sampleDuration;
    P.Plotting.ILimSample = ILimSample;
    P.Plotting.vLimSample = vLimSample;
    P.Plotting.eLimSample = eLimSample;
    P.Plotting.vMarkerChosenSpikes = vMarkerChosenSpikes;
    P.Plotting.colorPb = colorPb;
    P.Plotting.colorRr = colorRr;
    P.Plotting.colorRp = colorRp;
    P.Plotting.colorFm = colorFm;
    P.Plotting.markerTypePeak = markerTypePeak;
    P.Plotting.colorFirstPeakUsed = colorFirstPeakUsed;
    P.Plotting.colorFirstPeakNotUsed = colorFirstPeakNotUsed;
    P.Plotting.colorOtherPeak = colorOtherPeak;
    P.Plotting.colorValley = colorValley;
    P.Plotting.lineStylePeakAmplitude = lineStylePeakAmplitude;
    P.Plotting.colorPeakAmplitude = colorPeakAmplitude;
    P.Plotting.markerSizeFirstPeak = markerSizeFirstPeak;
    P.Plotting.markerSizeOtherPeak = markerSizeOtherPeak;
    P.Plotting.markerSizeValley = markerSizeValley;
    P.Plotting.colorEffectorOverlay = colorEffectorOverlay;
    P.Plotting.colorBasalCycleShade = colorBasalCycleShade;
    P.Plotting.colorSniffStartWin = colorSniffStartWin;
    P.Plotting.colorPreAnalysisShade = colorPreAnalysisShade;
    P.Plotting.transparencyOverlay = transparencyOverlay;
    P.Plotting.transparencyBasalCycle = transparencyBasalCycle;
    P.Plotting.transparencySniffStartWin = transparencySniffStartWin;
    P.Plotting.transparencyPreAnalysisShade = transparencyPreAnalysisShade;
    P.Plotting.colorT0 = colorT0;
    P.Plotting.colorT1 = colorT1;
    P.Plotting.colorReset = colorReset;
    P.Plotting.barFaceAlpha = barFaceAlpha;
    P.Plotting.markerSizeRaster = markerSizeRaster;
    P.Plotting.markerSizeChosenSpikes = markerSizeChosenSpikes;
    P.Plotting.lineWidthSampleTraces = lineWidthSampleTraces;
    P.Plotting.lineWidthOverlay = lineWidthOverlay;
    P.Plotting.lineWidthChosenSpikes = lineWidthChosenSpikes;
    P.Plotting.markerSizeJitter = markerSizeJitter;
    P.Plotting.jitterWidth = jitterWidth;
    P.Plotting.colorScatter = colorScatter;
    P.Plotting.markerTypeScatter = markerTypeScatter;
    P.Plotting.markerSizeScatter = markerSizeScatter;
    P.Plotting.markerLineWidthScatter = markerLineWidthScatter;
    P.Plotting.colorBestFit = colorBestFit;
    P.Plotting.lineStyleBestFit = lineStyleBestFit;
    P.Plotting.lineWidthBestFit = lineWidthBestFit;
    P.Plotting.colorThrOrig = colorThrOrig;
    P.Plotting.lineStyleThrOrig = lineStyleThrOrig;
    P.Plotting.lineWidthThrOrig = lineWidthThrOrig;
else
    P = paramsIn;
end

%% II. BACKWARD COMPATIBILITY CHECKS
% --- Handle backwards compatibility for old param files ---

% Add relativeTransitionTime if it doesn't exist (default for old files is 0.5)
if isfield(P, 'SquareInput') && isfield(P.SquareInput, 'Pb') && ...
   ~isfield(P.SquareInput.Pb, 'relativeTransitionTime')
    P.SquareInput.Pb.relativeTransitionTime = 0.5;
end

% Add relativeAnalysisStart if it doesn't exist (default for old files is 0)
if isfield(P, 'Analysis') && ~isfield(P.Analysis, 'relativeAnalysisStart')
    P.Analysis.relativeAnalysisStart = 0;
end

% Update test name field in parameters structure
if ~isfield(P, 'testName') && ~isempty(testName)
    % Create field with passed in test name
    P.testName = testName;

    % Show message
    post_message(sprintf('Test name %s save in params based on user request!', testName));
elseif ~isfield(P, 'testName') && isempty(testName)
    if ~isempty(paramFilePath)
        % Get the path of the directory and the filename without its extension
        [paramDirPath, paramFileBase, ~] = fileparts(paramFilePath);

        % Get the name of the directory containing the param file
        [~, paramDirName, ~] = fileparts(paramDirPath);

        % Check for the 'Test-<timestamp>-<testname>' pattern
        extractedTestName = ''; % Initialize as empty
        if startsWith(paramDirName, 'Test-')
            parts = split(paramDirName, '-');
            if numel(parts) >= 3
                % The test name is everything after the second hyphen.
                % Re-join in case the test name itself has hyphens.
                extractedTestName = strjoin(parts(3:end), '-');
            end
        end

        % If extraction was successful, use the result
        % As a fallback, use the filename without its extension
        if ~isempty(extractedTestName)
            testName = extractedTestName;
        else
            testName = paramFileBase;
        end

        % Show message
        post_message(sprintf('Test name extracted as %s!', testName));
    else
        % Make sure a test name or a parameter file path is passed in
        post_message('Please pass in a test name or parameter file path!');
        warning('No test name passed in or saved! P.testName set as empty!');
    end

    % Update with extracted test name
    P.testName = testName;
elseif isfield(P, 'testName') && ~isempty(testName)
    if ~strcmp(P.testName, testName)
        % Show message
        post_message(sprintf('Updating test name in params from "%s" to "%s" based on user request!', ...
                     P.testName, testName));

        % Update the test name in the parameters structure based on user request
        P.testName = testName;
    else
        % No update needed if they are the same
    end
elseif isfield(P, 'testName') && isempty(testName)
    % Update local testName variable (not needed unless this variable is used later)
    % testName = P.testName;
end

% Abbreviate 'SquareInput' to 'SqInpt' in testName for brevity
if isfield(P, 'testName') && ischar(P.testName) && ~isempty(P.testName)
    originalTestName = P.testName;
    P.testName = strrep(P.testName, 'SquareInput', 'SqInpt');
    if ~strcmp(originalTestName, P.testName)
        post_message(sprintf('Abbreviated test name from "%s" to "%s"', ...
                     originalTestName, P.testName));
    end
end

% Save or update the used seed into the parameters structure for record-keeping
P.seedNumber = seedNumberToUse;

% Save or update the used algorithm into the parameters structure
P.rngAlgorithm = rngAlgorithmToUse;

% Check Square Input parameters
if isfield(P, 'SquareInput') && isfield(P.SquareInput, 'Pb')
    % Check for and create the toRecompute flag
    if ~isfield(P.SquareInput.Pb, 'toRecompute')
        P.SquareInput.Pb.toRecompute = false;
    end

    % Convert P.SquareInput periods to dual periods for older files
    if isfield(P.SquareInput.Pb, 'period') && isfield(P.SquareInput.Pb, 'perVar')
        if ~isfield(P.SquareInput.Pb, 'period1')
            P.SquareInput.Pb.period1 = P.SquareInput.Pb.period;
            P.SquareInput.Pb.period2 = P.SquareInput.Pb.period;
            P.SquareInput.Pb.perVar1 = P.SquareInput.Pb.perVar;
            P.SquareInput.Pb.perVar2 = P.SquareInput.Pb.perVar;
        else
            post_message(sprintf('Old period %g will be ignored in favor of new period1 %g!', ...
                            P.SquareInput.Pb.period, P.SquareInput.Pb.period1));
        end

        % Remove old period and perVar
        P.SquareInput.Pb = rmfield(P.SquareInput.Pb, {'period', 'perVar'});
    end

    % Add perVarPerc fields if they don't exist
    if ~isfield(P.SquareInput.Pb, 'perVarPerc1')
        if P.SquareInput.Pb.period1 > 0
            P.SquareInput.Pb.perVarPerc1 = (P.SquareInput.Pb.perVar1 / P.SquareInput.Pb.period1) * 100;
        else
            P.SquareInput.Pb.perVarPerc1 = 0;
        end
    end
    if ~isfield(P.SquareInput.Pb, 'perVarPerc2')
        if P.SquareInput.Pb.period2 > 0
            P.SquareInput.Pb.perVarPerc2 = (P.SquareInput.Pb.perVar2 / P.SquareInput.Pb.period2) * 100;
        else
            P.SquareInput.Pb.perVarPerc2 = 0;
        end
    end

    % Enforce consistency between perVar and perVarPerc (perVarPerc has priority)
    tolerance = 1e-9; % Tolerance for floating point comparison
    
    % Check and update perVar1
    expectedPerVar1 = (P.SquareInput.Pb.perVarPerc1 / 100) * P.SquareInput.Pb.period1;
    if abs(P.SquareInput.Pb.perVar1 - expectedPerVar1) > tolerance
        oldPerVar1 = P.SquareInput.Pb.perVar1;
        P.SquareInput.Pb.perVar1 = expectedPerVar1;
        post_message(sprintf(['Warning: Inconsistent perVar1 (%.2f) and perVarPerc1 (%.2f). ' ...
            'Updating perVar1 to %.2f based on percentage.'], oldPerVar1, P.SquareInput.Pb.perVarPerc1, expectedPerVar1));
    end

    % Check and update perVar2
    expectedPerVar2 = (P.SquareInput.Pb.perVarPerc2 / 100) * P.SquareInput.Pb.period2;
    if abs(P.SquareInput.Pb.perVar2 - expectedPerVar2) > tolerance
        oldPerVar2 = P.SquareInput.Pb.perVar2;
        P.SquareInput.Pb.perVar2 = expectedPerVar2;
        post_message(sprintf(['Warning: Inconsistent perVar2 (%.2f) and perVarPerc2 (%.2f). ' ...
            'Updating perVar2 to %.2f based on percentage.'], oldPerVar2, P.SquareInput.Pb.perVarPerc2, expectedPerVar2));
    end
end

% If P.ExtCurrent.(pool) is not a struct, convert it to the new format
if isfield(P, 'ExtCurrent')
    poolNamesForCurrents = fieldnames(P.ExtCurrent);
    for iPool = 1:length(poolNamesForCurrents)
        pool = poolNamesForCurrents{iPool};
        if ~isstruct(P.ExtCurrent.(pool))
            % Old format detected (scalar value), convert to new struct format
            msg = sprintf('Old P.ExtCurrent format detected for %s pool. Converting to new format with zero variance.', pool);
            post_message(msg);
            oldMean = P.ExtCurrent.(pool);
            P.ExtCurrent.(pool) = struct('mean', oldMean, 'stdPop', 0, 'stdTime', 0, 'tauTime', 0);
        end
    end
end

% Add analysis and plotting parameters if they don't exist
if ~isfield(P, 'Analysis')
    P.Analysis.relativeAnalysisStart = 0.25;
    P.Analysis.amplitudeDefinition = amplitudeDefinition;
    P.Analysis.fundFreqRange = fundFreqRange;
    P.Analysis.fCutoffWhisk = fCutoffWhisk;
    P.Analysis.fCutoffRelToFund = fCutoffRelToFund;
    P.Analysis.filterOrderWhisk = filterOrderWhisk;
    P.Analysis.promThresholdPercWhisk = promThresholdPercWhisk;
    P.Analysis.minPeakPromWhisk = minPeakPromWhisk;
    P.Analysis.maxWhiskDurationMs = maxWhiskDurationMs;
    P.Analysis.minPeakDistanceMsWhisk = minPeakDistanceMsWhisk;
    P.Analysis.sniffFreqThreshold = sniffFreqThreshold;
    P.Analysis.basalFreqThreshold = basalFreqThreshold;
    P.Analysis.nWhisksSniffStartToAnalyze = nWhisksSniffStartToAnalyze;
    P.Analysis.minWhisksBasalRespToAnalyze = minWhisksBasalRespToAnalyze;
    P.Analysis.maxWhisksBasalRespToAnalyze = maxWhisksBasalRespToAnalyze;
    P.Analysis.nCorrToAnalyze = nCorrToAnalyze;
    P.Analysis.whiskDirForPhase = whiskDirForPhase;
    P.Analysis.breathOnsetLatencyMs = breathOnsetLatencyMs;
end
if ~isfield(P, 'Plotting')
    P.Plotting.sampleDuration = sampleDuration;
    P.Plotting.ILimSample = ILimSample;
    P.Plotting.vLimSample = vLimSample;
    P.Plotting.eLimSample = eLimSample;
    P.Plotting.vMarkerChosenSpikes = vMarkerChosenSpikes;
    P.Plotting.colorPb = colorPb;
    P.Plotting.colorRr = colorRr;
    P.Plotting.colorRp = colorRp;
    P.Plotting.colorFm = colorFm;
    P.Plotting.markerTypePeak = markerTypePeak;
    P.Plotting.colorFirstPeakUsed = colorFirstPeakUsed;
    P.Plotting.colorFirstPeakNotUsed = colorFirstPeakNotUsed;
    P.Plotting.colorOtherPeak = colorOtherPeak;
    P.Plotting.colorValley = colorValley;
    P.Plotting.lineStylePeakAmplitude = lineStylePeakAmplitude;
    P.Plotting.colorPeakAmplitude = colorPeakAmplitude;
    P.Plotting.markerSizeFirstPeak = markerSizeFirstPeak;
    P.Plotting.markerSizeOtherPeak = markerSizeOtherPeak;
    P.Plotting.markerSizeValley = markerSizeValley;
    P.Plotting.colorEffectorOverlay = colorEffectorOverlay;
    P.Plotting.colorBasalCycleShade = colorBasalCycleShade;
    P.Plotting.colorSniffStartWin = colorSniffStartWin;
    P.Plotting.colorPreAnalysisShade = colorPreAnalysisShade;
    P.Plotting.transparencyOverlay = transparencyOverlay;
    P.Plotting.transparencyBasalCycle = transparencyBasalCycle;
    P.Plotting.transparencySniffStartWin = transparencySniffStartWin;
    P.Plotting.transparencyPreAnalysisShade = transparencyPreAnalysisShade;
    P.Plotting.colorT0 = colorT0;
    P.Plotting.colorT1 = colorT1;
    P.Plotting.colorReset = colorReset;
    P.Plotting.barFaceAlpha = barFaceAlpha;
    P.Plotting.markerSizeRaster = markerSizeRaster;
    P.Plotting.markerSizeChosenSpikes = markerSizeChosenSpikes;
    P.Plotting.lineWidthSampleTraces = lineWidthSampleTraces;
    P.Plotting.lineWidthOverlay = lineWidthOverlay;
    P.Plotting.lineWidthChosenSpikes = lineWidthChosenSpikes;
    P.Plotting.markerSizeJitter = markerSizeJitter;
    P.Plotting.jitterWidth = jitterWidth;
    P.Plotting.colorScatter = colorScatter;
    P.Plotting.markerTypeScatter = markerTypeScatter;
    P.Plotting.markerSizeScatter = markerSizeScatter;
    P.Plotting.markerLineWidthScatter = markerLineWidthScatter;
    P.Plotting.colorBestFit = colorBestFit;
    P.Plotting.lineStyleBestFit = lineStyleBestFit;
    P.Plotting.lineWidthBestFit = lineWidthBestFit;
    P.Plotting.colorThrOrig = colorThrOrig;
    P.Plotting.lineStyleThrOrig = lineStyleThrOrig;
    P.Plotting.lineWidthThrOrig = lineWidthThrOrig;
end

% Add pre-analysis shade parameters if they don't exist
if isfield(P, 'Plotting') && ~isfield(P.Plotting, 'colorPreAnalysisShade')
    P.Plotting.colorPreAnalysisShade = colorPreAnalysisShade;
end
if isfield(P, 'Plotting') && ~isfield(P.Plotting, 'transparencyPreAnalysisShade')
    P.Plotting.transparencyPreAnalysisShade = transparencyPreAnalysisShade;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%