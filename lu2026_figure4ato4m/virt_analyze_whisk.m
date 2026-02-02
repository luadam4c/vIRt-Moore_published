function [analysis, analysisString] = virt_analyze_whisk (P, State, tVecMs, varargin)
%% Analyze whisking motion of a vIRt simulation
% Usage: [analysis, analysisString] = virt_analyze_whisk (P, State, tVecMs, varargin)
% Explanation:
%       This function analyzes the whisking output of a vIRt simulation.
%
% Example(s):
%       [analysis, analysisString] = virt_analyze_whisk(P, State, tVecMs);
%
% Outputs:
%       analysis    - A structured array containing all analysis results.
%                   specified as a structured array
%       analysisString  - A formatted string of analysis results for display.
%                   specified as a string scalar
%
% Arguments:
%       P           - The parameters structure from virt_moore.m.
%                   specified as a structure
%       State       - The State structure from virt_moore.m.
%                   specified as a structure
%       tVecMs      - The time vector in milliseconds.
%                   specified as a numeric vector
%       varargin    - 'ToComputeWhiskAmp': Whether to compute amplitude statistics.
%                   must be a logical scalar
%                   default == true
%                   - 'ToComputeWhiskInterval': Whether to compute interval statistics.
%                   must be a logical scalar
%                   default == true
%                   - 'ToPostMessage': Whether to display messages.
%                   must be a logical scalar
%                   default == true
%                   - 'MessageHandler': Function handle for displaying messages to a GUI.
%                   must be a function handle
%                   default == []
%
% Requires:
%       \Shared\Code\Adams_Functions\compute_combined_trace.m
%       \Shared\Code\Adams_Functions\extract_elements.m
%       \Shared\Code\Adams_Functions\extract_subvectors.m
%       \Shared\Code\Adams_Functions\force_matrix.m
%       \Shared\Code\Adams_Functions\parse_oscillation.m
%       \Shared\Code\Adams_Functions\virt_detect_whisk_analysis_windows.m
%
% Used by:
%       \Shared\Code\vIRt-Moore\virt_moore.m
%

% File History:
% 2025-09-24 Pulled Code from virt_moore.m by Gemini
% 2025-10-01 Added sniff start window detection
% 2025-10-05 Refactored basal cycle analysis to save a basalRespCycleTable
% 2025-10-05 Now restricts vectors to after tAnalysisStartMs before parse_oscillation.
% 2025-10-06 Refactored to use virt_detect_whisk_analysis_windows.m
% 2025-10-09 Now stores and passes seedNumber
% 2025-10-15 Now performs Fisher z-transformation on correlations

%% Compatibility Check
% Check MATLAB version for function flag support
v = ver('matlab');
versionStr = regexp(v.Version, '^\d+\.\d+', 'match', 'once');
matlabVersion = str2double(versionStr);
supportsOmitNaN = matlabVersion >= 8.5;     % R2015a introduced 'omitnan' for std() and mean()

%% Hard-coded parameters

%% Default values for optional arguments
toComputeWhiskAmpDefault = true;
toComputeWhiskIntervalDefault = true;
toPostMessageDefault = true;
messageHandlerDefault = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Functions

% Function to post messages to command window and/or GUI
function post_message(msg)
    if toPostMessage
        % Print to command window
        fprintf('%s\n', msg);

        % Also send to GUI if handler exists
        if ~isempty(messageHandler)
            messageHandler(msg);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ToComputeWhiskAmp', toComputeWhiskAmpDefault, ...
    @(x) validateattributes(x, {'logical'}, {'scalar'}));
addParameter(iP, 'ToComputeWhiskInterval', toComputeWhiskIntervalDefault, ...
    @(x) validateattributes(x, {'logical'}, {'scalar'}));
addParameter(iP, 'ToPostMessage', toPostMessageDefault, ...
    @(x) validateattributes(x, {'logical'}, {'scalar'}));
addParameter(iP, 'MessageHandler', messageHandlerDefault);

% Read from the Input Parser
parse(iP, varargin{:});
toComputeWhiskAmp = iP.Results.ToComputeWhiskAmp;
toComputeWhiskInterval = iP.Results.ToComputeWhiskInterval;
toPostMessage = iP.Results.ToPostMessage;
messageHandler = iP.Results.MessageHandler;

%% Preparation
% Initialize output
analysis = struct;
analysisCellStr = {};

% Extract relevant parameters from the main P structure for whisk analysis
dt = P.dt;
seedNumber = P.seedNumber;
analysisParams = P.Analysis;
relativeAnalysisStart = analysisParams.relativeAnalysisStart;
nCorrToAnalyze = analysisParams.nCorrToAnalyze;
amplitudeDefinition = analysisParams.amplitudeDefinition;
fundFreqRange = analysisParams.fundFreqRange;
fCutoffWhisk = analysisParams.fCutoffWhisk;
filterOrderWhisk = analysisParams.filterOrderWhisk;
promThresholdPercWhisk = analysisParams.promThresholdPercWhisk;
minPeakPromWhisk = analysisParams.minPeakPromWhisk;
maxWhiskDurationMs = analysisParams.maxWhiskDurationMs;
minPeakDistanceMsWhisk = analysisParams.minPeakDistanceMsWhisk;
whiskDirForPhase = analysisParams.whiskDirForPhase;

% Extract whisk position vector and PB input pulse start times
effectorVec = State.effectorPosition;
pulseStartTimesMs = P.SquareInput.Pb.timesStart;
nRp = P.N.Rp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Do the job
% Analyze the whisking movement and detect whisk peaks
try
    % Starting message
    post_message('Starting Whisk Movement Peak Detection ...');

    % --- Whisk Peak Detection ---
    analysisCellStr{end+1} = '--- Whisk Movement Analysis ---';

    % Get the analysis start time in ms and save it in sec
    tAnalysisStartMs = P.simDur * relativeAnalysisStart;

    % Find the index from which to start the analysis
    if tAnalysisStartMs > 0
        startIndex = find(tVecMs >= tAnalysisStartMs, 1);
        post_message(sprintf('Restricting analysis to t >= %.1f ms.', tAnalysisStartMs));
    else
        startIndex = 1;
    end

    % Create subvectors for analysis if necessary
    if startIndex > 1
        subVecs = extract_subvectors({effectorVec, tVecMs}, ...
                                    'EndPoints', [startIndex, Inf], ...
                                    'ForceInRange', true);
        effectorVecForAnalysis = subVecs{1};
        tVecMsForAnalysis = subVecs{2};
    else
        effectorVecForAnalysis = effectorVec;
        tVecMsForAnalysis = tVecMs;
    end

    % Use parse_oscillation to find peaks, valleys, and fundamental frequency
    [peakTable, valleyTable, otherResults] = ...
        parse_oscillation(effectorVecForAnalysis, ...
                        'TimeVec', tVecMsForAnalysis, 'SamplingIntervalMs', dt, ...
                        'AmpMode', amplitudeDefinition, ...
                        'FundFreqRange', fundFreqRange, ...
                        'FilterCutoffs', fCutoffWhisk, ...
                        'FilterOrder', filterOrderWhisk, ...
                        'PromThresholdPerc', promThresholdPercWhisk, ...
                        'MinPeakProminence', minPeakPromWhisk, ...
                        'MinValleyProminence', minPeakPromWhisk, ...
                        'MinPeakDistanceMs', minPeakDistanceMsWhisk, ...
                        'MaxDurationForAmplitudeMs', maxWhiskDurationMs);

    % Extract the fundamental frequency of the oscillation (Hz)
    freqFundamental = otherResults.freqFundamental;

    % Print and save results to output
    lineStr1 = sprintf('Detected Fundamental Frequency: %.2f Hz', freqFundamental);
    post_message(lineStr1);
    analysisCellStr{end+1} = lineStr1;

    % Extract peak and valley times and values from tables
    allPeakTimes = peakTable.peakTime;
    allPeakValues = peakTable.peakValue;
    allValleyTimes = valleyTable.valleyTime;
    allValleyValues = valleyTable.valleyValue;
    
    % Group whisk data by simulation pulse cycles
    nPulses = length(pulseStartTimesMs);
    pulseNumbers = (1:nPulses)';
    pulseCycleStartTimes = nan(nPulses, 1);
    pulseCycleEndTimes = nan(nPulses, 1);
    nPeaksByPulse = zeros(nPulses, 1);
    peakTablesByPulse = cell(nPulses, 1);
    firstPeakTimesByPulse = nan(nPulses, 1);
    firstPreValleyTimesByPulse = nan(nPulses, 1);
    firstPostValleyTimesByPulse = nan(nPulses, 1);
    firstPeakAmplitudesByPulse = nan(nPulses, 1);
    peakTimesByPulse = cell(nPulses, 1);
    preValleyTimesByPulse = cell(nPulses, 1);
    postValleyTimesByPulse = cell(nPulses, 1);
    peakAmplitudesByPulse = cell(nPulses, 1);
    for iPulse = 1:nPulses
        % Define the time window for the current pulse cycle
        startTime = pulseStartTimesMs(iPulse);
        if iPulse < nPulses
            endTime = pulseStartTimesMs(iPulse + 1);
        else
            endTime = tVecMs(end);
        end

        % Find which peaks from the table fall into this cycle
        peaksInPulseCycle = peakTable(allPeakTimes >= startTime & ...
                                        allPeakTimes < endTime, :);
            
        % Count number of peaks in each cycle
        nPeaksThisPulse = height(peaksInPulseCycle);

        % Store the results for the first peak for this cycle
        if nPeaksThisPulse > 0
            firstPeakTimesByPulse(iPulse) = peaksInPulseCycle.peakTime(1);
            firstPreValleyTimesByPulse(iPulse) = peaksInPulseCycle.preValleyTime(1);
            firstPostValleyTimesByPulse(iPulse) = peaksInPulseCycle.postValleyTime(1);
            firstPeakAmplitudesByPulse(iPulse) = peaksInPulseCycle.amplitude(1);
        end

        pulseCycleStartTimes(iPulse) = startTime;
        pulseCycleEndTimes(iPulse) = endTime;
        nPeaksByPulse(iPulse) = nPeaksThisPulse;
        peakTablesByPulse{iPulse} = peaksInPulseCycle;
        peakTimesByPulse{iPulse} = peaksInPulseCycle.peakTime;
        peakAmplitudesByPulse{iPulse} = peaksInPulseCycle.amplitude;
        preValleyTimesByPulse{iPulse} = peaksInPulseCycle.preValleyTime;
        postValleyTimesByPulse{iPulse} = peaksInPulseCycle.postValleyTime;
    end
    
    % Detect Sniff and Basal Analysis Windows
    [sniffStartWinTable, basalRespCycleTable] = ...
        virt_detect_whisk_analysis_windows(peakTable, valleyTable, analysisParams, ...
            'SniffDefinitionMode', 'pulseCycle', ...
            'BasalDefinitionMode', 'pulseCycle', ...
            'PulseCycleStartTimes', pulseCycleStartTimes, ...
            'PulseCycleEndTimes', pulseCycleEndTimes, ...
            'PulseNumbers', pulseNumbers, ...
            'SeedNumber', seedNumber, ...
            'TVecMs', tVecMs, ...
            'NPeaksByPulse', nPeaksByPulse, ...
            'PeakAmplitudesByPulse', peakAmplitudesByPulse, ...
            'PeakTablesByPulse', peakTablesByPulse, ...
            'ToPostMessage', toPostMessage, ...
            'MessageHandler', messageHandler);

    % Extract results from returned tables
    nSniffStarts = height(sniffStartWinTable);
    nBasals = height(basalRespCycleTable);

    % --- Save analysis parameters to the main analysis structure ---
    analysis.PbInput.pulseWidth = P.SquareInput.Pb.pulseWidth;
    analysis.PbInput.pulseStartTimes = pulseStartTimesMs;

    analysis.whisk.seedNumber = seedNumber;
    analysis.whisk.relativeAnalysisStart = relativeAnalysisStart;
    analysis.whisk.nCorrToAnalyze = nCorrToAnalyze;
    analysis.whisk.amplitudeDefinition = amplitudeDefinition;
    analysis.whisk.fundFreqRange = fundFreqRange;
    analysis.whisk.fCutoffWhisk = fCutoffWhisk;
    analysis.whisk.filterOrderWhisk = filterOrderWhisk;
    analysis.whisk.promThresholdPercWhisk = promThresholdPercWhisk;
    analysis.whisk.minPeakPromWhisk = minPeakPromWhisk;
    analysis.whisk.maxWhiskDurationMs = maxWhiskDurationMs;
    analysis.whisk.minPeakDistanceMsWhisk = minPeakDistanceMsWhisk;
    analysis.whisk.whiskDirForPhase = whiskDirForPhase;

    % --- Save all results to the main analysis structure ---
    % Scalars
    analysis.whisk.tAnalysisStartSec = tAnalysisStartMs / 1000;
    analysis.whisk.fundamentalFrequency = freqFundamental;
    analysis.whisk.nBasals = nBasals;
    analysis.whisk.nSniffStarts = nSniffStarts;

    % Vectors with length nPulses
    analysis.whisk.pulseNumbers = pulseNumbers;
    analysis.whisk.pulseCycleStartTimes = pulseCycleStartTimes;
    analysis.whisk.pulseCycleEndTimes = pulseCycleEndTimes;
    analysis.whisk.firstPeakTimesByPulse = firstPeakTimesByPulse;
    analysis.whisk.firstPreValleyTimesByPulse = firstPreValleyTimesByPulse;
    analysis.whisk.firstPostValleyTimesByPulse = firstPostValleyTimesByPulse;
    analysis.whisk.firstPeakAmplitudesByPulse = firstPeakAmplitudesByPulse;
    analysis.whisk.nPeaksByPulse = nPeaksByPulse;

    % Cell arrays with length nPulses
    analysis.whisk.peakTimesByPulse = peakTimesByPulse;
    analysis.whisk.preValleyTimesByPulse = preValleyTimesByPulse;
    analysis.whisk.postValleyTimesByPulse = postValleyTimesByPulse;
    analysis.whisk.peakAmplitudesByPulse = peakAmplitudesByPulse;
    analysis.whisk.peakTablesByPulse = peakTablesByPulse;

    % Tables
    analysis.whisk.sniffStartWinTable = sniffStartWinTable;
    analysis.whisk.basalRespCycleTable = basalRespCycleTable;

    % Other Vectors
    analysis.whisk.allPeakTimes = allPeakTimes;
    analysis.whisk.allPeakValues = allPeakValues;
    analysis.whisk.allValleyTimes = allValleyTimes;
    analysis.whisk.allValleyValues = allValleyValues;
    analysis.whisk.pulseNumbersSniffStart = sniffStartWinTable.pulseNumber;
    analysis.whisk.pulseNumbersBasal = basalRespCycleTable.pulseNumber;

    % Ending message
    post_message('Whisk Movement Peak Detection Complete!');
catch ME
    errorStr = sprintf('An error occurred during whisk peak detection: %s', ME.message);
    post_message(errorStr); 
    post_message(ME.getReport);
    analysisCellStr{end+1} = errorStr;
end

% --- Amplitude Analysis (Logarithmic Decrement & Correlation) ---
if toComputeWhiskAmp
    try
        % Starting message
        post_message('Starting Amplitude Analysis...');

        % --- First Peak Amplitude Analysis ---
        analysisCellStr{end+1} = ''; % Add a blank line
        analysisCellStr{end+1} = '--- First Peak Amplitude Analysis ---';
        
        % Print all first peak amplitudes
        post_message(sprintf('Pulse#\tStartTime(ms)\tPeak1Amp\tNPeaks\n'));
        analysisCellStr{end+1} = sprintf('Pulse#\tStartTime(ms)\tPeak1Amp\tNPeaks');
        for iPulse = 1:nPulses
            lineStrTable = sprintf('%d\t%.0f\t\t%.2f\t\t%d', ...
                pulseNumbers(iPulse), pulseCycleStartTimes(iPulse), ...
                firstPeakAmplitudesByPulse(iPulse), nPeaksByPulse(iPulse));
            post_message(lineStrTable);
            analysisCellStr{end+1} = lineStrTable;
        end

        % Get the mean breath-induced whisk peak amplitude (NaNs are ignored)
        if supportsOmitNaN
            meanFirstPeakAmpToAnalyze = mean(firstPeakAmplitudesByPulse, 'omitnan');
        else
            meanFirstPeakAmpToAnalyze = nanmean(firstPeakAmplitudesByPulse);
        end

        % Save to analysis structure output
        analysis.whisk.meanFirstPeakAmplitudeToAnalyze = meanFirstPeakAmpToAnalyze;

        % Print and save results to output
        lineStr2 = sprintf('Mean Amplitude of First Peak (analyzed section): %.2f', meanFirstPeakAmpToAnalyze);
        post_message(sprintf('%s\n', lineStr2));
        analysisCellStr{end+1} = ''; % Add a blank line
        analysisCellStr{end+1} = lineStr2;

        % --- Basal Respiration Logarithmic Decrement Analysis ---
        post_message('Starting Whisk Logarithmic Decrement Analysis for Basal Respiration ...');
        analysisCellStr{end+1} = ''; % Add a blank line
        analysisCellStr{end+1} = '--- Whisk Logarithmic Decrement Analysis (Basal Respiration) ---';

        % Extract log decrements from the table if it exists and is not empty
        if isfield(analysis.whisk, 'basalRespCycleTable') && ~isempty(analysis.whisk.basalRespCycleTable)
            % Extract log decrements from the table
            logDecrementsByBasal = basalRespCycleTable.whiskLogDecrements;

            % Compute the maximum number of log decrements
            maxDecrementsBasal = max(cellfun(@numel, logDecrementsByBasal));

            % Aggregate and compute statistics
            if maxDecrementsBasal > 0
                % Create a matrix to hold all decrements (rows=cycles, cols=decrement order)
                logDecrementsMatrixBasal = transpose(force_matrix(logDecrementsByBasal, ...
                                                'CombineMethod', 'leftAdjustPad'));

                % Compute stats for each decrement order, ignoring NaNs
                meanLogDecrementsBasal = compute_combined_trace(logDecrementsMatrixBasal', 'mean');
                stdLogDecrementsBasal = compute_combined_trace(logDecrementsMatrixBasal', 'std');
                stderrLogDecrementsBasal = compute_combined_trace(logDecrementsMatrixBasal', 'stderr');
                err95LogDecrementsBasal = compute_combined_trace(logDecrementsMatrixBasal', 'err95');
                lower95LogDecrementsBasal = meanLogDecrementsBasal - err95LogDecrementsBasal;
                upper95LogDecrementsBasal = meanLogDecrementsBasal + err95LogDecrementsBasal;
                
                % Calculate geometric mean of amplitude ratios from mean log decrements
                avgWhiskAmpRatiosBasal = exp(meanLogDecrementsBasal);
                
                % Save results to analysis structure
                analysis.whisk.maxDecrementsBasal = maxDecrementsBasal;
                analysis.whisk.logDecrementsByBasal = logDecrementsByBasal;
                analysis.whisk.logDecrementsMatrixBasal = logDecrementsMatrixBasal;
                analysis.whisk.meanLogDecrementsBasal = meanLogDecrementsBasal;
                analysis.whisk.stdLogDecrementsBasal = stdLogDecrementsBasal;
                analysis.whisk.stderrLogDecrementsBasal = stderrLogDecrementsBasal;
                analysis.whisk.err95LogDecrementsBasal = err95LogDecrementsBasal;
                analysis.whisk.lower95LogDecrementsBasal = lower95LogDecrementsBasal;
                analysis.whisk.upper95LogDecrementsBasal = upper95LogDecrementsBasal;
                analysis.whisk.avgWhiskAmpRatiosBasal = avgWhiskAmpRatiosBasal;
                
                % Format and print results
                analysisCellStr{end+1} = sprintf('Order#\tmean[ln(A_n+1/A_n)]\tSEM\t\tGeoMean[Ratio]');
                for iDec = 1:maxDecrementsBasal
                    lineStr = sprintf('%d\t%.3f\t\t\t%.3f\t\t%.3f', iDec, ...
                        meanLogDecrementsBasal(iDec), stderrLogDecrementsBasal(iDec), avgWhiskAmpRatiosBasal(iDec));
                    analysisCellStr{end+1} = lineStr;
                end
            else
                % Ending message
                post_message('Not enough sequential whisks in basal respiration cycles to compute decrements!');
                analysisCellStr{end+1} = 'Not enough sequential whisks in basal respiration cycles to compute decrements.';
            end
        else
            post_message('No basal respiration cycles found for analysis!');
            analysisCellStr{end+1} = 'No basal respiration cycles found for analysis.';
        end

        % Ending message
        post_message('Basal Respiration Whisk Logarithmic Decrement Analysis complete!');
    
        % --- Basal Respiration Amplitude Correlation Analysis ---
        post_message('Starting Basal Respiration Successive Whisk Amplitude Correlation Analysis ...');
        analysisCellStr{end+1} = ''; % Add a blank line
        analysisCellStr{end+1} = '--- Successive Whisk Amplitude Correlation (Basal Respiration) ---';
    
        % Extract amplitudes from the table if it exists and is not empty
        if isfield(analysis.whisk, 'basalRespCycleTable') && ~isempty(analysis.whisk.basalRespCycleTable)
            % Extract peak amplitudes from the new table
            peakAmplitudesByBasal = basalRespCycleTable.whiskPeakAmplitudes;

            % Count the maximum number of peak orders across all basal respiration cycles
            maxPeaksForAnalysis = max(cellfun(@numel, peakAmplitudesByBasal));

            % Determine number of correlations to compute
            nCorrelations = min(nCorrToAnalyze, maxPeaksForAnalysis - 1);

            % Analyze amplitude correlations
            if maxPeaksForAnalysis > 1
                % Create a matrix of whisk amplitudes from valid cycles (rows=cycles, cols=whisk peak order)
                whiskAmplitudesMatrix = transpose(force_matrix(peakAmplitudesByBasal, ...
                                                    'CombineMethod', 'leftAdjustPad'));
                
                % Compute correlations
                corrCoeffs = nan(nCorrelations, 1);
                pValues = nan(nCorrelations, 1);
                zScores = nan(nCorrelations, 1);
                if nCorrelations > 0
                    analysisCellStr{end+1} = sprintf('Amps\tCorrCoeff\tpValue\tFisherZScore');
                    for iCorr = 1:nCorrelations
                        % Extract amplitude data for whisk N and whisk N+1
                        ampCurrent = whiskAmplitudesMatrix(:, iCorr);
                        ampNext = whiskAmplitudesMatrix(:, iCorr + 1);
                        
                        % Remove pairs with NaN values to ensure proper correlation
                        toKeep = ~isnan(ampCurrent) & ~isnan(ampNext);
                        ampCurrentValid = ampCurrent(toKeep);
                        ampNextValid = ampNext(toKeep);

                        % Compute correlation with p values and z-score
                        if numel(ampCurrentValid) > 2 % Need at least 3 points to correlate
                            [corrMatrix, pValuesThis] = corrcoef(ampCurrentValid, ampNextValid);
                            corrCoeffs(iCorr) = corrMatrix(1, 2);
                            pValues(iCorr) = pValuesThis(1, 2);

                            % Apply the Fisher Z-transformation for normalized statistics later
                            % Clamp r to prevent Inf from atanh
                            r = corrMatrix(1, 2);
                            rClamped = max(min(r, 1 - 1e-9), -1 + 1e-9);
                            zScores(iCorr) = atanh(rClamped);
                        end
            
                        % Format and print results
                        lineStr = sprintf('A%d-A%d\t%.3f\t\t%.4f\t%.3f', iCorr, iCorr+1, ...
                            corrCoeffs(iCorr), pValues(iCorr), zScores(iCorr));
                        analysisCellStr{end+1} = lineStr;
                    end
                else
                    post_message('Not enough sequential whisks in basal respiration cycles to compute correlations!');
                    analysisCellStr{end+1} = 'Not enough sequential whisks in basal respiration cycles to compute correlations.';
                end

                % Save results to analysis structure
                analysis.whisk.peakAmplitudesByBasal = peakAmplitudesByBasal;
                analysis.whisk.nCorrToAnalyzeBasal = nCorrToAnalyze;
                analysis.whisk.whiskAmpMatrixBasal = whiskAmplitudesMatrix;
                analysis.whisk.whiskAmpCorrsBasal = corrCoeffs;
                analysis.whisk.whiskAmpCorrPValuesBasal = pValues;
                analysis.whisk.whiskAmpCorrZScoresBasal = zScores;
            else
                post_message('Not enough sequential whisks in basal respiration cycles to compute correlations!');
                analysisCellStr{end+1} = 'Not enough sequential whisks in basal respiration cycles to compute correlations.';
            end
        else
            post_message('No basal respiration cycles found for analysis!');
            analysisCellStr{end+1} = 'No basal respiration cycles found for analysis.';
        end

        % Ending message
        post_message('Basal Respiration Successive Whisk Amplitude Correlation Analysis complete!');

        % --- Sniff Start Logarithmic Decrement Analysis ---
        post_message('Starting Whisk Logarithmic Decrement Analysis for Sniff Start Windows...');
        analysisCellStr{end+1} = ''; % Add a blank line
        analysisCellStr{end+1} = '--- Whisk Logarithmic Decrement Analysis (Sniff Start) ---';

        % Extract amplitudes and log decrements from the table if it exists and is not empty
        if isfield(analysis.whisk, 'sniffStartWinTable') && ~isempty(analysis.whisk.sniffStartWinTable)
            % Extract amplitudes and log decrements
            peakAmplitudesBySniff = analysis.whisk.sniffStartWinTable.whiskPeakAmplitudes;
            logDecrementsBySniff = analysis.whisk.sniffStartWinTable.whiskLogDecrements;

            % Compute the maximum number of log decrements
            maxDecrementsSniff = max(cellfun(@numel, logDecrementsBySniff));

            if maxDecrementsSniff > 0
                % Create a matrix to hold all decrements (rows=windows, cols=decrement order)
                logDecrementsMatrixSniff = transpose(force_matrix(logDecrementsBySniff, ...
                                                        'CombineMethod', 'leftAdjustPad'));

                % Compute stats for each decrement order, ignoring NaNs
                meanLogDecrementsSniff = compute_combined_trace(logDecrementsMatrixSniff', 'mean');
                stdLogDecrementsSniff = compute_combined_trace(logDecrementsMatrixSniff', 'std');
                stderrLogDecrementsSniff = compute_combined_trace(logDecrementsMatrixSniff', 'stderr');
                err95LogDecrementsSniff = compute_combined_trace(logDecrementsMatrixSniff', 'err95');
                lower95LogDecrementsSniff = meanLogDecrementsSniff - err95LogDecrementsSniff;
                upper95LogDecrementsSniff = meanLogDecrementsSniff + err95LogDecrementsSniff;
                
                % Calculate geometric mean of amplitude ratios from mean log decrements
                avgWhiskAmpRatiosSniff = exp(meanLogDecrementsSniff);

                % Save results to the analysis structure
                analysis.whisk.peakAmplitudesBySniff = peakAmplitudesBySniff;
                analysis.whisk.maxDecrementsSniff = maxDecrementsSniff;
                analysis.whisk.logDecrementsBySniff = logDecrementsBySniff;
                analysis.whisk.logDecrementsMatrixSniff = logDecrementsMatrixSniff;
                analysis.whisk.meanLogDecrementsSniff = meanLogDecrementsSniff;
                analysis.whisk.stdLogDecrementsSniff = stdLogDecrementsSniff;
                analysis.whisk.stderrLogDecrementsSniff = stderrLogDecrementsSniff;
                analysis.whisk.err95LogDecrementsSniff = err95LogDecrementsSniff;
                analysis.whisk.lower95LogDecrementsSniff = lower95LogDecrementsSniff;
                analysis.whisk.upper95LogDecrementsSniff = upper95LogDecrementsSniff;
                analysis.whisk.avgWhiskAmpRatiosSniff = avgWhiskAmpRatiosSniff;

                % Format and print the results to the analysis string
                analysisCellStr{end+1} = sprintf('Order#\tmean[ln(A_n+1/A_n)]\tSEM\t\tGeoMean[Ratio]');
                for iDec = 1:maxDecrementsSniff
                    lineStr = sprintf('%d\t%.3f\t\t\t%.3f\t\t%.3f', iDec, ...
                        meanLogDecrementsSniff(iDec), stderrLogDecrementsSniff(iDec), avgWhiskAmpRatiosSniff(iDec));
                    analysisCellStr{end+1} = lineStr;
                end
            else
                post_message('Not enough sequential whisks in sniff start windows to compute decrements!');
                analysisCellStr{end+1} = 'Not enough sequential whisks in sniff start windows to compute decrements.';
            end
        else
            post_message('No sniff start windows found for analysis!');
            analysisCellStr{end+1} = 'No sniff start windows found for analysis.';
        end

        % Ending message
        post_message('Sniff Start Window Logarithmic Decrement Analysis complete!');

    catch ME
        errorStr = sprintf('An error occurred during whisk amplitude analysis: %s', ME.message);
        post_message(errorStr);
        analysisCellStr{end+1} = errorStr;
    end
end

% --- Interval and Phase Analysis ---
if toComputeWhiskInterval
    try
        % Starting message
        post_message('Starting Peak Timing Analysis ...');
        analysisCellStr{end+1} = ''; % Add a blank line
        analysisCellStr{end+1} = '--- Peak Timing Analysis ---';

        % Extract peak times and print
        post_message(sprintf('Pulse#\tCycleStart(ms)\tPeakTimes(ms)\n'));
        analysisCellStr{end+1} = sprintf('Pulse#\tCycleStart(ms)\tPeakTimes(ms)');
        for iPulse = 1:nPulses
            lineStr = sprintf('%d\t%.0f\t\t', iPulse, pulseCycleStartTimes(iPulse));
            lineStr = [lineStr, sprintf('%.0f ', peakTimesByPulse{iPulse})];
            post_message(lineStr);
            analysisCellStr{end+1} = lineStr;
        end

        % Compute Rp-nearest-spike → effector‐peak latencies (signed, relative to peak) 
        %       for each cycle
        %   Only consider spikes within the surrounding peak valley times
        % TODO: It may be more accurate to use the first spike of the nearest burst instead
        RpEffectorLatencies = cell(nPulses, 1);
        chosenSpikeTimesByPulse = cell(nPulses, 1);
        for iPulse = 1:nPulses
            % Get effector peak times for this cycle
            peakTimes = peakTimesByPulse{iPulse};
            preValleyTimes = preValleyTimesByPulse{iPulse};
            postValleyTimes = postValleyTimesByPulse{iPulse};
            nPeaks = nPeaksByPulse(iPulse);
        
            % Compute latencies from nearest Rp spike to effector peak 
            %   for each Rp cell and each peak within this cycle
            latencies = nan(nRp, nPeaks);
            chosenSpikeTimes = nan(nRp, nPeaks);
            for iPeak = 1:nPeaks
                % Get peak times and surrounding valley times
                tPeak = peakTimes(iPeak);
                tPreValley = preValleyTimes(iPeak);
                tPostValley = postValleyTimes(iPeak);

                % Run through each cell
                for iCell = 1:nRp
                    % Get all the spike times for this cell
                    tSpikesThisCell = State.Rp.spikeTimes{iCell};

                    % Restrict to this peak
                    tSpikesThisPeak = tSpikesThisCell(tSpikesThisCell >= tPreValley & ...
                                                        tSpikesThisCell < tPostValley);

                    % If no spikes surrounding this peak, this cell is noncontributory and latency is NaN
                    if isempty(tSpikesThisPeak)
                        continue
                    end

                    % Find the spike whose time is closest to tPeak
                    [~, iSpike] = min(abs(tSpikesThisPeak - tPeak));

                    % Record the chosen spike time and its latency
                    tChosenSpike = tSpikesThisPeak(iSpike);
                    chosenSpikeTimes(iCell, iPeak) = tChosenSpike;

                    % Record signed latency: (neuron_spike - effector_peak)
                    latencies(iCell, iPeak) = tChosenSpike - tPeak;
                end
            end

            % Save the chosen Rp Spikes and latencies for this cycle    
            chosenSpikeTimesByPulse{iPulse} = chosenSpikeTimes;
            RpEffectorLatencies{iPulse} = latencies;
        end

        % Compute std(|latency|) per peak, per cycle
        if supportsOmitNaN
            stdAbsLatencyByPulse = cellfun(@(x) std(abs(x), 0, 1, 'omitnan'), RpEffectorLatencies, 'UniformOutput', false);
        else
            stdAbsLatencyByPulse = cellfun(@(x) nanstd(abs(x), 0, 1), RpEffectorLatencies, 'UniformOutput', false);
        end

        % Print and save results to output
        analysisCellStr{end+1} = ''; % Add a blank line
        post_message(sprintf('vIRt-protraction Nearest Spike to Whisk Peak Latency Jitter:'));
        analysisCellStr{end+1} = '--- vIRt-protraction Nearest Spike to Whisk Peak Latency Jitter ---';
        post_message(sprintf('Pulse#\tstd|latency|(ms)\n'));
        analysisCellStr{end+1} = sprintf('Cycle#\tstd|latency|(ms)');
        for iPulse = 1:nPulses
            lineStr = sprintf('%.1f ', stdAbsLatencyByPulse{iPulse});
            lineStrTable = sprintf('%d\t%s', pulseNumbers(iPulse), lineStr);
            post_message(lineStrTable);
            analysisCellStr{end+1} = lineStrTable;
        end

        % Compute inter-peak intervals for each cycle
        interPeakIntervalsByPulse = cellfun(@diff, peakTimesByPulse, 'UniformOutput', false);

        % Compute the maximum number of inter-peak intervals in a cycle
        maxIntervals = max(cellfun(@numel, interPeakIntervalsByPulse));
        
        % Group "1st interval" across cycles, "2nd interval" across cycles, …
        interPeakIntervalsByOrder = extract_elements(interPeakIntervalsByPulse, 'all');

        % Compute mean and standard deviation for each interval order
        meanPeakIntervalByOrder = cellfun(@nanmean, interPeakIntervalsByOrder);
        stdPeakIntervalByOrder = cellfun(@nanstd, interPeakIntervalsByOrder);

        % Print and save results to output
        post_message('Whisk peak intervals by order within a pulse cycle:');
        analysisCellStr{end+1} = ''; % Add a blank line
        analysisCellStr{end+1} = '--- Whisk peak intervals by order within a pulse cycle ---';
        post_message(sprintf('Order#\tmean(ms)\tstd(ms)\n'));
        analysisCellStr{end+1} = sprintf('Order#\tmean(ms)\tstd(ms)');
        for iInterval = 1:maxIntervals
            lineStr = sprintf('%d->%d\t%.1f\t\t%.1f', ...
                    iInterval, iInterval + 1, ...
                    meanPeakIntervalByOrder(iInterval), stdPeakIntervalByOrder(iInterval));
            post_message(lineStr);
            analysisCellStr{end+1} = lineStr;
        end

        % Save to analysis structure output
        analysis.whisk.chosenRpSpikeTimesByPulse = chosenSpikeTimesByPulse;
        analysis.whisk.stdAbsLatencyByPulse = stdAbsLatencyByPulse;
        analysis.whisk.interPeakIntervalsByPulse = interPeakIntervalsByPulse;
        analysis.whisk.interPeakIntervalsByOrder = interPeakIntervalsByOrder;
        analysis.whisk.meanPeakIntervalByOrder = meanPeakIntervalByOrder;
        analysis.whisk.stdPeakIntervalByOrder = stdPeakIntervalByOrder;

        % --- Basal Respiration Phase Response Curve Analysis ---
        % Starting message
        post_message('Starting Basal Respiration Phase Response Curve  Analysis ...');
        analysisCellStr{end+1} = '';
        analysisCellStr{end+1} = '--- Phase Response Curve ---';
        analysisCellStr{end+1} = sprintf('Pulse#\tPhaseReset (rad)\tPhaseChange (rad)');

        % Extract amplitudes from the table if it exists and is not empty
        if isfield(analysis.whisk, 'basalRespCycleTable') && ~isempty(analysis.whisk.basalRespCycleTable)
            % Extract phase data from the basal cycle table
            pulseNumbersBasal = basalRespCycleTable.pulseNumber;
            phaseResets = basalRespCycleTable.phaseReset;
            phaseChangesWhisk = basalRespCycleTable.phaseChangeWhisk;

            % Print and save results to output
            for iBasal = 1:nBasals
                if ~isnan(phaseResets(iBasal)) && ~isnan(phaseChangesWhisk(iBasal))
                    lineStr = sprintf('%d\t%.2f\t\t%.2f', pulseNumbersBasal(iBasal), phaseResets(iBasal), phaseChangesWhisk(iBasal));
                    analysisCellStr{end+1} = lineStr;
                end
            end

            % Save to analysis structure output
            analysis.whisk.phaseResets = phaseResets;
            analysis.whisk.phaseChangesWhisk = phaseChangesWhisk;
            analysis.whisk.breathOnsetTimes = basalRespCycleTable.breathOnsetTime;
            analysis.whisk.relativeResetTimes = basalRespCycleTable.relativeResetTime;
            analysis.whisk.preIEIsWhiskBefore = basalRespCycleTable.preIEIWhiskBefore;
            analysis.whisk.preIEIsWhiskAfter = basalRespCycleTable.preIEIWhiskAfter;
            analysis.whisk.eventTimesWhiskBefore = basalRespCycleTable.eventTimeWhiskBefore;
            analysis.whisk.eventTimesTwoWhisksBefore = basalRespCycleTable.eventTimeTwoWhisksBefore;
        else
            post_message('No basal respiration cycles found for analysis!');
            analysisCellStr{end+1} = 'No basal respiration cycles found for analysis.';
        end

        % Ending message
        post_message('Peak Timing and Basal Respiration Phase Response Analysis complete!');
    catch ME
        errorStr = sprintf('An error occurred during whisk timing analysis: %s', ME.message);
        post_message(errorStr);
        analysisCellStr{end+1} = errorStr;
    end
end

% Add the cell array of strings into a single string with newlines
analysisString = strjoin(analysisCellStr, '\n');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
