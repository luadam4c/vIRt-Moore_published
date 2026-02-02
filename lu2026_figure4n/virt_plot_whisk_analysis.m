function handles = virt_plot_whisk_analysis (analysis, pPlot, varargin)
%% Plotting whisk analysis results from a vIRt simulation
% Usage: handles = virt_plot_whisk_analysis (analysis, pPlot, varargin)
% Explanation:
%       This function plots the whisk analysis results from virt_analyze_whisk.m.
%       It can create new figures or update existing ones if handles are provided.
%
% Example(s):
%       handles = virt_plot_whisk_analysis(analysis, pPlot, 'Handles', handles);
%
% Outputs:
%       handles     - A structure of handles to the generated or updated figures.
%
% Arguments:
%       analysis    - The analysis structure from virt_analyze_whisk.m.
%                   specified as a structure
%       pPlot       - The plotting parameters with color specified as RGB values
%                   specified as a structure
%       varargin    - 'Handles': A structure of graphics handles. If fields for specific
%                     figures (e.g., handles.figJitter) exist and are valid,
%                     the figures will be updated instead of recreated.
%                   must be a structure
%                   default == []
%                   - 'PlotWhiskDetectionFlag': Whether to plot whisk detection markers.
%                   must be a logical scalar
%                   default == true
%                   - 'PlotSpikeDetectionFlag': Whether to plot spike detection markers.
%                   must be a logical scalar
%                   default == true
%                   - 'PlotJitterFlag': Whether to plot the jitter plot.
%                   must be a logical scalar
%                   default == true
%                   - 'PlotScatterFlag': Whether to plot the scatter plots.
%                   must be a logical scalar
%                   default == true
%                   - 'PlotPRCFlag': Whether to plot the phase response curve.
%                   must be a logical scalar
%                   default == true
%                   - 'AlwaysNew': Whether to always create new figures instead of updating.
%                   must be a logical scalar
%                   default == false
%                   - 'ShowFigure': Whether to show the figure.
%                   must be a logical scalar
%                   default == true
%                   - 'OutDir': Output directory for saving figures.
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'FigTypes': File type(s) for saving figures.
%                   must be a string scalar, character vector, or a cell array of character vectors
%                   default == {'png'}
%                   - 'TimeStamp': Timestamp for naming files.
%                   must be a string scalar or a character vector
%                   default == create_time_stamp
%                   - 'ToPostMessage': Whether to display messages.
%                   must be a logical scalar
%                   default == true
%                   - 'MessageHandler': Function handle for displaying messages.
%                   must be a function handle
%                   default == []
%                   - 'ToSaveOutput': Whether to save the figures.
%                   must be a logical scalar
%                   default == true
%                   - 'NameModel': The base name for the model in saved files.
%                   must be a string scalar or character vector
%                   default == 'virt_moore'
%
% Requires:
%       \Shared\Code\Adams_Functions\create_time_stamp.m
%       \Shared\Code\Adams_Functions\plot_vertical_shade.m
%       \Shared\Code\Adams_Functions\virt_plot_amplitude_correlation.m
%       \Shared\Code\Adams_Functions\virt_plot_jitter.m
%       \Shared\Code\Adams_Functions\virt_plot_phase_response.m
%
% Used by:
%       \Shared\Code\vIRt-Moore\virt_moore.m

% File History:
% 2025-09-24 - Pulled Code from virt_moore.m by Gemini
% 2025-09-25 - Restored update logic, annotations, and message handling.
%              Variable extraction and save_all_figtypes implemented.
% 2025-10-05 - Updated to use basalRespCycleTable for plotting basal shades.
%              Updated jitter plot to use variables with 'Basal' suffix.
% 2025-10-05 - Added vertical shade for pre-analysis time.
% 2025-10-06 - Refactored by Gemini to use standalone plotting functions for creation.
% 2025-10-06 - Refactored by Gemini to move update logic to standalone functions.
% 2025-10-17 - Made 'AlwaysNew' an optional argument, default == false
% 2025-10-17 - Made 'ShowFigure' an optional argument, default == true
% 2026-01-14 - Fixed bug where error bars in PRC plots were not being cleared by Gemini

%% Hard-coded parameters
% File names for saved figures
nameJitterFig = 'Jitter';
nameScatterFig = 'Scatter';
namePRCFig = 'PhaseResponseCurve';

% Figure titles
figTitleJitter = 'Whisk Logarithmic Decrements';
figTitleScatter = 'Successive Whisk Amplitude Correlations';
figTitlePRC = 'Whisk Phase Response Curve';

%% Default values for optional arguments
handlesDefault = [];
plotWhiskDetectionFlagDefault = true;
plotSpikeDetectionFlagDefault = true;
plotJitterFlagDefault = true;
plotScatterFlagDefault = true;
plotPRCFlagDefault = true;
alwaysNewDefault = false;
showFigureDefault = true;
dirOutDefault = pwd;
figTypesDefault = {'png'};
timeStampDefault = create_time_stamp;
toPostMessageDefault = true;
messageHandlerDefault = [];
toSaveOutputDefault = true;
nameModelDefault = 'virt_moore';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Functions

% Function to post messages to command window and/or GUI
function post_message(msg)
    if toPostMessage
        % Always print to command window
        fprintf('%s\n', msg);
        
        % Also send to GUI if handler exists
        if ~isempty(messageHandler)
            messageHandler(msg);
        end
    end
end

% Function to plot chosen spikes on the vIRt-protraction sample trace
function plot_chosen_spikes_sample(axVRp, analysis, pPlot)
    if ~isfield(analysis.whisk, 'chosenRpSpikeTimesByPulse')
        post_message('No chosen spikes to plot!');
        return
    end

    % Extract plotting parameters
    vMarkerChosenSpikes = pPlot.vMarkerChosenSpikes;
    markerSizeChosenSpikes = pPlot.markerSizeChosenSpikes;
    lineWidthChosenSpikes = pPlot.lineWidthChosenSpikes;

    % Extract chosen spike times from analysis structure
    chosenSpikeTimesByPulse = analysis.whisk.chosenRpSpikeTimesByPulse;
    
    % Extract all chosen spike times for the first Rp neuron
    allChosenSpikesCell1 = cellfun(@(x) x(1, :), ...
        chosenSpikeTimesByPulse, 'UniformOutput', false);
    allChosenSpikesCell1 = horzcat(allChosenSpikesCell1{:});
    allChosenSpikesCell1Sec = allChosenSpikesCell1(~isnan(allChosenSpikesCell1)) / 1000;

    % Add red 'x' markers to the Rp sample trace plot
    if ~isempty(allChosenSpikesCell1Sec)
        post_message('Plotting chosen spikes on Rp sample trace...');
        hold(axVRp, 'on');
        plot(axVRp, allChosenSpikesCell1Sec, vMarkerChosenSpikes, 'rx', ...
                'MarkerSize', markerSizeChosenSpikes, ...
                'LineWidth', lineWidthChosenSpikes, ...
                'Tag', 'analysis_points');
        hold(axVRp, 'off');
        post_message('Finished plotting chosen spikes.');
    end
end

% Function to plot detection on effector trace
function plot_detection_on_effector(axEffector, analysis, pPlot)
    % This function plots whisk analysis results (peaks, valleys, phases)
    % onto an existing effector trace axes.

    % Extract plotting parameters from parameters struct for convenience
    markerTypePeak = pPlot.markerTypePeak;
    colorFirstPeakUsed = pPlot.colorFirstPeakUsed;
    colorFirstPeakNotUsed = pPlot.colorFirstPeakNotUsed;
    colorOtherPeak = pPlot.colorOtherPeak;
    colorValley = pPlot.colorValley;
    lineStylePeakAmplitude = pPlot.lineStylePeakAmplitude;
    colorPeakAmplitude = pPlot.colorPeakAmplitude;
    markerSizeFirstPeak = pPlot.markerSizeFirstPeak;
    markerSizeOtherPeak = pPlot.markerSizeOtherPeak;
    markerSizeValley = pPlot.markerSizeValley;
    colorBasalCycleShade = pPlot.colorBasalCycleShade;
    colorSniffStartWin = pPlot.colorSniffStartWin;
    colorPreAnalysisShade = pPlot.colorPreAnalysisShade;
    transparencyBasalCycle = pPlot.transparencyBasalCycle;
    transparencySniffStartWin = pPlot.transparencySniffStartWin;
    transparencyPreAnalysisShade = pPlot.transparencyPreAnalysisShade;

    colorT0 = pPlot.colorT0;
    colorT1 = pPlot.colorT1;
    colorReset = pPlot.colorReset;
    barFaceAlpha = pPlot.barFaceAlpha;

    % Extract detection results from analysis struct for convenience
    tAnalysisStartSec = analysis.whisk.tAnalysisStartSec;
    pulseNumbers = analysis.whisk.pulseNumbers;
    firstPeakAmplitudesByPulse = analysis.whisk.firstPeakAmplitudesByPulse;
    sniffStartWinTable = analysis.whisk.sniffStartWinTable;
    if isfield(analysis.whisk, 'basalRespCycleTable') && ~isempty(analysis.whisk.basalRespCycleTable)
        basalRespCycleStartTimes = analysis.whisk.basalRespCycleTable.basalRespCycleStartTime;
        basalRespCycleEndTimes = analysis.whisk.basalRespCycleTable.basalRespCycleEndTime;
    else
        basalRespCycleStartTimes = [];
        basalRespCycleEndTimes = [];
    end
    allPeakTimes = analysis.whisk.allPeakTimes;
    allPeakValues = analysis.whisk.allPeakValues;
    allValleyTimes = analysis.whisk.allValleyTimes;
    allValleyValues = analysis.whisk.allValleyValues;
    peakTimesByPulse = analysis.whisk.peakTimesByPulse;
    peakAmplitudesByPulse = analysis.whisk.peakAmplitudesByPulse;

    % Hold on to the axis
    hold(axEffector, 'on');
    
    % Plot a gray vertical shade for the pre-analysis time
    if tAnalysisStartSec > 0
        plot_vertical_shade([0, tAnalysisStartSec], ...
                            'Color', colorPreAnalysisShade, ...
                            'FaceAlpha', transparencyPreAnalysisShade, ...
                            'AxesHandle', axEffector, 'Tag', 'analysis_start_shade');
    end

    % Plot sniff start windows as vertical shades
    if isfield(analysis.whisk, 'sniffStartWinTable') && ~isempty(sniffStartWinTable)
        % Extract the boundaries
        sniffStartWinStartTime = sniffStartWinTable.sniffStartWinStartTime;
        sniffStartWinEndTime = sniffStartWinTable.sniffStartWinEndTime;
        sniffBoundaries = [sniffStartWinStartTime, sniffStartWinEndTime];

        % Plot windows
        plot_vertical_shade(sniffBoundaries', 'Color', colorSniffStartWin, ...
                            'FaceAlpha', transparencySniffStartWin, ...
                            'AxesHandle', axEffector, ...
                            'Tag', 'analysis_sniff_shades');
    end

    % Put basal start and end times together
    basalBoundaries = [basalRespCycleStartTimes, basalRespCycleEndTimes];

    % Plot the basal respiration cycles as vertical shades
    if ~isempty(basalBoundaries)
        plot_vertical_shade(basalBoundaries', 'Color', colorBasalCycleShade, ...
                            'FaceAlpha', transparencyBasalCycle, ...
                            'AxesHandle', axEffector, ...
                            'Tag', 'analysis_basal_shades');
    end

    % Plot all detected valleys
    plot(axEffector, allValleyTimes, allValleyValues, ...
         'o', 'MarkerFaceColor', colorValley, ...
         'MarkerEdgeColor', colorValley, ...
         'Tag', 'analysis_points', 'MarkerSize', markerSizeValley);

    % Plot peaks, coloring the first in each cycle differently
    nPulses = numel(pulseNumbers);
    for iPulse = 1:nPulses
        peakTimesThisPulse = peakTimesByPulse{iPulse};
        peakAmpsThisPulse = peakAmplitudesByPulse{iPulse};
        if ~isempty(peakTimesThisPulse)
            % Obtain first peak amplitude for this cycle
            firstPeakAmplitude = firstPeakAmplitudesByPulse(iPulse);

            % If amplitude is NaN, the peak is not used for analysis
            firstPeakUsed = ~isnan(firstPeakAmplitude);

            % Decide on the color of the first peak
            if firstPeakUsed
                colorFirstPeak = colorFirstPeakUsed;
            else
                colorFirstPeak = colorFirstPeakNotUsed;
            end

            % Plot the first peak
            firstPeakTime = peakTimesThisPulse(1);
            firstPeakValue = allPeakValues(allPeakTimes == firstPeakTime);
            plot(axEffector, firstPeakTime, firstPeakValue, ...
                 markerTypePeak, 'MarkerFaceColor', colorFirstPeak, ...
                 'MarkerEdgeColor', colorFirstPeak, ...
                 'Tag', 'analysis_points', 'MarkerSize', markerSizeFirstPeak);

            % Plot the first peak amplitude as a red dashed line if used
            if firstPeakUsed
                plot(axEffector, [firstPeakTime, firstPeakTime], ...
                     [firstPeakValue - firstPeakAmplitude, ...
                        firstPeakValue], lineStylePeakAmplitude, 'Color', colorPeakAmplitude, ...
                     'Tag', 'analysis_points');
            end

            % Plot the rest of the peaks in cyan
            if numel(peakTimesThisPulse) > 1
                otherPeakTimes = peakTimesThisPulse(2:end);
                otherPeakValues = allPeakValues(ismember(allPeakTimes, otherPeakTimes));
                plot(axEffector, otherPeakTimes, otherPeakValues, ...
                     markerTypePeak, 'MarkerFaceColor', colorOtherPeak, ...
                     'MarkerEdgeColor', colorOtherPeak, ...
                     'Tag', 'analysis_points', 'MarkerSize', markerSizeOtherPeak);

                % Get the amplitudes for these other peaks from the cycle-specific table
                otherPeakAmplitudes = peakAmpsThisPulse(2:end);
                
                % Plot the subsequent peak amplitudes as a pink dashed line if used
                for iOtherPeak = 1:numel(otherPeakTimes)
                    peakTime = otherPeakTimes(iOtherPeak);
                    peakValue = otherPeakValues(iOtherPeak);
                    peakAmplitude = otherPeakAmplitudes(iOtherPeak);
                    
                    % Only plot if the amplitude is a valid number
                    if ~isnan(peakAmplitude)
                        plot(axEffector, [peakTime, peakTime], ...
                             [peakValue - peakAmplitude, peakValue], ...
                             lineStylePeakAmplitude, 'Color', colorPeakAmplitude, ... 
                             'Tag', 'analysis_points');
                    end
                end
            end
        end
    end

    % Plot phase response bars for each basal cycle
    if isfield(analysis.whisk, 'nBasals') && analysis.whisk.nBasals > 0
        % Extract necessary phase data
        breathOnsetTimes = analysis.whisk.breathOnsetTimes;
        eventTimesWhiskBefore = analysis.whisk.eventTimesWhiskBefore;
        eventTimesTwoWhisksBefore = analysis.whisk.eventTimesTwoWhisksBefore;
        preIEIsWhiskAfter = analysis.whisk.preIEIsWhiskAfter;
        nBasals = analysis.whisk.nBasals;

        % Loop through each basal cycle with valid phase data
        for iBasal = 1:nBasals
            % Get data for this specific cycle
            eventTimeWhiskBefore = eventTimesWhiskBefore(iBasal);
            eventTimeTwoWhisksBefore = eventTimesTwoWhisksBefore(iBasal);
            preIEIWhiskAfter = preIEIsWhiskAfter(iBasal);
            breathOnsetTime = breathOnsetTimes(iBasal);

            % Check if all necessary data for this cycle is valid (not NaN)
            if ~isnan(eventTimeWhiskBefore) && ~isnan(eventTimeTwoWhisksBefore) && ...
               ~isnan(preIEIWhiskAfter) && ~isnan(breathOnsetTime)
               
                % Set y-coordinates for the bars relative to the effector plot range
                yLimits = get(axEffector, 'YLim');
                barYLowEffector = yLimits(1) + 0.05 * diff(yLimits); % Position at 5% of y-range
                barHeightEffector = 0.05 * diff(yLimits); % Bar height is 5% of y-range

                yT0T1Rect = barYLowEffector + barHeightEffector * [-1, 1, 1, -1] / 2;
                yResetRect = barYLowEffector + barHeightEffector + barHeightEffector * [-1, 1, 1, -1] / 2;

                % 1. Plot T0 (unperturbed inter-whisk-interval)
                xT0Rect = [eventTimeTwoWhisksBefore, eventTimeTwoWhisksBefore, eventTimeWhiskBefore, eventTimeWhiskBefore];
                fill(axEffector, xT0Rect, yT0T1Rect, colorT0, 'FaceAlpha', barFaceAlpha, 'Tag', 'analysis_points', 'EdgeColor', 'none');

                % 2. Plot T1 (perturbed inter-whisk-interval)
                peakTimeAfter = eventTimeWhiskBefore + preIEIWhiskAfter;
                xT1Rect = [eventTimeWhiskBefore, eventTimeWhiskBefore, peakTimeAfter, peakTimeAfter];
                fill(axEffector, xT1Rect, yT0T1Rect, colorT1, 'FaceAlpha', barFaceAlpha, 'Tag', 'analysis_points', 'EdgeColor', 'none');

                % 3. Plot reset time (breath onset)
                xResetRect = [eventTimeWhiskBefore, eventTimeWhiskBefore, breathOnsetTime, breathOnsetTime];
                fill(axEffector, xResetRect, yResetRect, colorReset, 'FaceAlpha', barFaceAlpha, 'Tag', 'analysis_points', 'EdgeColor', 'none');
            end
        end
    end
    
    % Hold off
    hold(axEffector, 'off');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Handles', handlesDefault);
addParameter(iP, 'PlotWhiskDetectionFlag', plotWhiskDetectionFlagDefault, ...
    @(x) validateattributes(x, {'logical'}, {'scalar'}));
addParameter(iP, 'PlotSpikeDetectionFlag', plotSpikeDetectionFlagDefault, ...
    @(x) validateattributes(x, {'logical'}, {'scalar'}));
addParameter(iP, 'PlotJitterFlag', plotJitterFlagDefault, ...
    @(x) validateattributes(x, {'logical'}, {'scalar'}));
addParameter(iP, 'PlotScatterFlag', plotScatterFlagDefault, ...
    @(x) validateattributes(x, {'logical'}, {'scalar'}));
addParameter(iP, 'PlotPRCFlag', plotPRCFlagDefault, ...
    @(x) validateattributes(x, {'logical'}, {'scalar'}));
addParameter(iP, 'AlwaysNew', alwaysNewDefault, ...
    @(x) validateattributes(x, {'logical'}, {'scalar'}));
addParameter(iP, 'ShowFigure', showFigureDefault, ...
    @(x) validateattributes(x, {'logical'}, {'scalar'}));
addParameter(iP, 'OutDir', dirOutDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FigTypes', figTypesDefault, ...
    @(x) ischar(x) || isstring(x) || iscell(x));
addParameter(iP, 'TimeStamp', timeStampDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'ToPostMessage', toPostMessageDefault, ...
    @(x) validateattributes(x, {'logical'}, {'scalar'}));
addParameter(iP, 'MessageHandler', messageHandlerDefault);
addParameter(iP, 'ToSaveOutput', toSaveOutputDefault, ...
    @(x) validateattributes(x, {'logical'}, {'scalar'}));
addParameter(iP, 'NameModel', nameModelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, varargin{:});
handles = iP.Results.Handles;
plotWhiskDetectionFlag = iP.Results.PlotWhiskDetectionFlag;
plotSpikeDetectionFlag = iP.Results.PlotSpikeDetectionFlag;
plotJitterFlag = iP.Results.PlotJitterFlag;
plotScatterFlag = iP.Results.PlotScatterFlag;
plotPRCFlag = iP.Results.PlotPRCFlag;
alwaysNew = iP.Results.AlwaysNew;
showFigure = iP.Results.ShowFigure;
dirOut = iP.Results.OutDir;
figTypes = iP.Results.FigTypes;
timeStamp = iP.Results.TimeStamp;
toPostMessage = iP.Results.ToPostMessage;
messageHandler = iP.Results.MessageHandler;
toSaveOutput = iP.Results.ToSaveOutput;
nameModel = iP.Results.NameModel;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preparation
% If handles is empty, create it
if isempty(handles)
    handles = struct;
end

% Create figure names
figNameJitter = sprintf('%s_%s_%s', nameModel, nameJitterFig, timeStamp);
figNameScatter = sprintf('%s_%s_%s', nameModel, nameScatterFig, timeStamp);
figNamePRC = sprintf('%s_%s_%s', nameModel, namePRCFig, timeStamp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Do the job
% Update FIGURE 1
if plotSpikeDetectionFlag
    % Plot chosen spikes on the Rp sample trace
    if isfield(handles, 'ax1') && numel(handles.ax1) >= 4 && isgraphics(handles.ax1(4))
        plot_chosen_spikes_sample(handles.ax1(4), analysis, pPlot);
    else
        post_message('Warning: Cannot plot chosen spikes because Rp sample trace axis handle is invalid.');
    end
end

% Update FIGURE 2
if plotWhiskDetectionFlag
    % Plot detection results on the main raster plot figure
    if isfield(handles, 'ax2') && numel(handles.ax2) >= 2 && isgraphics(handles.ax2(2))
        post_message('Plotting whisk analysis results on effector trace...');
        plot_detection_on_effector(handles.ax2(2), analysis, pPlot);
        post_message('Finished plotting whisk analysis detection.');
    else
        post_message('Warning: Cannot plot detections because effector axis handle is invalid.');
    end
end

% FIGURE 3: Logarithmic Decrement Jitter Plot
if plotJitterFlag && isfield(analysis.whisk, 'basalRespCycleTable') && ...
        isfield(analysis.whisk, 'logDecrementsMatrixBasal')
    % Extract variables
    logDecrementsMatrix = analysis.whisk.logDecrementsMatrixBasal;
    [nBasals, maxDecrements] = size(logDecrementsMatrix);
    
    if maxDecrements > 0 && nBasals > 0
        dataTable = analysis.whisk.basalRespCycleTable;
        if ~isempty(dataTable)
            % Add a dummy grouping column for single runs
            if ~isfield(dataTable, 'repetitionNumber')
                dataTable.repetitionNumber = ones(height(dataTable), 1);
            end

            % Decide whether to create or update
            if ~alwaysNew && isfield(handles, 'figJitter') && ...
                    isgraphics(handles.figJitter)
                post_message('Updating logarithmic decrement jitter plot...');

                % Package specific handles for the update function
                handlesToUpdate.fig = handles.figJitter;
                handlesToUpdate.axJitter = handles.axJitter;
                handlesToUpdate.hJitter = handles.hJitter;
                handlesToUpdate.hErrorBars = handles.hErrorBars;
                handlesToUpdate.hMeans = handles.hMeans;
                handlesToUpdate.pTextJitter = handles.pTextJitter;
                handlesToUpdate.hNull = handles.hNull;
                handlesToUpdate.sigMarkerJitter = handles.sigMarkerJitter;
                handlesToUpdate.transformedLabels = handles.transformedLabels;
                updateOrCreate = 'updated';
            else
                post_message('Plotting logarithmic decrement jitter plot...');
                handlesToUpdate = [];
                updateOrCreate = 'created';
            end
            
            % Call the unified plotting function
            [~, hJitterOut] = virt_plot_jitter(dataTable, pPlot, ...
                                    'DataColumn', 'whiskLogDecrements', ...
                                    'DataMode', 'LogDecrement', ...
                                    'Handles', handlesToUpdate, ...
                                    'MaxOrders', maxDecrements, ...
                                    'FigTitle', figTitleJitter, ...
                                    'FigName', figNameJitter, ...
                                    'OutDir', dirOut, 'FigTypes', figTypes, ...
                                    'ToSaveOutput', toSaveOutput, ...
                                    'ShowFigure', showFigure);

            % Store/update the returned handles
            handles.figJitter = hJitterOut.fig;
            handles.axJitter = hJitterOut.axJitter;
            handles.hJitter = hJitterOut.hJitter;
            handles.hErrorBars = hJitterOut.hErrorBars;
            handles.hMeans = hJitterOut.hMeans;
            handles.hNull = hJitterOut.hNull;
            handles.pTextJitter = hJitterOut.pTextJitter;
            handles.sigMarkerJitter = hJitterOut.sigMarkerJitter;
            handles.transformedLabels = hJitterOut.transformedLabels;

            % Post message
            post_message(['Jitter plot ', updateOrCreate, '.']);
        end
    end
end


% FIGURE 4: Successive Amplitude Scatter Plots
if plotScatterFlag && isfield(analysis.whisk, 'basalRespCycleTable') && ...
    isfield(analysis.whisk, 'whiskAmpMatrixBasal')
    % Extract variables for clarity
    whiskAmplitudesMatrix = analysis.whisk.whiskAmpMatrixBasal;
    nCorrToAnalyze = analysis.whisk.nCorrToAnalyze;
    [~, maxPeaks] = size(whiskAmplitudesMatrix);
    nCorrelations = min(nCorrToAnalyze, maxPeaks - 1);
    
    if nCorrelations > 0
        dataTable = analysis.whisk.basalRespCycleTable;
        if ~isempty(dataTable)
            % Add a dummy grouping column for single runs
            if ~isfield(dataTable, 'repetitionNumber')
                dataTable.repetitionNumber = ones(height(dataTable), 1);
            end
            
            % Decide whether to create or update
            if ~alwaysNew && isfield(handles, 'figScatter') && ...
                    isgraphics(handles.figScatter) && ...
                    isfield(handles, 'hScatters') && ...
                    numel(handles.hScatters) == nCorrelations
                post_message('Updating successive whisk amplitude scatter plots...');
                handlesToUpdate.fig = handles.figScatter;
                handlesToUpdate.axScatter = handles.axScatter;
                handlesToUpdate.hScatters = handles.hScatters;
                handlesToUpdate.hBestFits = handles.hBestFits;
                handlesToUpdate.hTextBestFit = handles.hTextBestFit;
                handlesToUpdate.hThrOrigs = handles.hThrOrigs;
                handlesToUpdate.hTextThrOrig = handles.hTextThrOrig;
                updateOrCreate = 'updated';
            else
                post_message('Plotting successive whisk amplitude scatter plots...');
                handlesToUpdate = [];
                updateOrCreate = 'created';
            end

            % Call the unified plotting function
            [~, hScatterOut] = virt_plot_amplitude_correlation(dataTable, pPlot, ...
                                    'Handles', handlesToUpdate, ...
                                    'NCorrelations', nCorrelations, ...
                                    'FigTitle', figTitleScatter, ...
                                    'FigName', figNameScatter, ...
                                    'OutDir', dirOut, 'FigTypes', figTypes, ...
                                    'ToSaveOutput', toSaveOutput, ...
                                    'ShowFigure', showFigure);

            % Store handles
            handles.figScatter = hScatterOut.fig;
            handles.axScatter = hScatterOut.axScatter;
            handles.hScatters = hScatterOut.hScatters;
            handles.hBestFits = hScatterOut.hBestFits;
            handles.hTextBestFit = hScatterOut.hTextBestFit;
            handles.hThrOrigs = hScatterOut.hThrOrigs;
            handles.hTextThrOrig = hScatterOut.hTextThrOrig;

            % Post message
            post_message(['Scatter plots ', updateOrCreate, '.']);
        end
    end
end

% FIGURE 5: Phase Response Curve
if plotPRCFlag && isfield(analysis.whisk, 'basalRespCycleTable')
    basalTable = analysis.whisk.basalRespCycleTable;
    toKeep = ~isnan(basalTable.phaseReset) & ~isnan(basalTable.phaseChangeWhisk);
    
    if sum(toKeep) > 0
        dataTable = analysis.whisk.basalRespCycleTable;
        if ~isempty(dataTable)
            % Add a dummy grouping column for single runs
            if ~isfield(dataTable, 'repetitionNumber')
                dataTable.repetitionNumber = ones(height(dataTable), 1);
            end

            % Decide whether to create or update
            if ~alwaysNew && isfield(handles, 'figPRC') && ...
                    isgraphics(handles.figPRC)
                post_message('Updating phase response curve...');
                handlesToUpdate.fig = handles.figPRC;
                handlesToUpdate.axPRC = handles.axPRC;
                handlesToUpdate.hPRCScatter = handles.hPRCScatter;
                handlesToUpdate.hPRCRegLine = handles.hPRCRegLine;
                handlesToUpdate.hPRCRegText = handles.hPRCRegText;
                handlesToUpdate.hBinError = handles.hPRCBinError;
                handlesToUpdate.hBinSig = handles.hPRCBinSig;
                updateOrCreate = 'updated';
            else
                post_message('Plotting phase response curve...');
                handlesToUpdate = [];
                updateOrCreate = 'created';
            end

            % Plot the phase response
            [~, hPrcOut] = virt_plot_phase_response(dataTable, pPlot, ...
                            'Handles', handlesToUpdate, ...
                            'WhiskDir', analysis.whisk.whiskDirForPhase, ...
                            'FigTitle', figTitlePRC, 'FigName', figNamePRC, ...
                            'OutDir', dirOut, 'FigTypes', figTypes, ...
                            'ToSaveOutput', toSaveOutput, ...
                            'ShowFigure', showFigure);

            % Store handles
            handles.figPRC = hPrcOut.fig;
            handles.axPRC = hPrcOut.axPRC;
            handles.hPRCScatter = hPrcOut.hPRCScatter;
            handles.hPRCRegLine = hPrcOut.hPRCRegLine;
            handles.hPRCRegText = hPrcOut.hPRCRegText;
            handles.hPRCBinError = hPrcOut.hBinError;
            handles.hPRCBinSig = hPrcOut.hBinSig;

            % Post message
            post_message(['Phase Response Curve plot ', updateOrCreate, '.']);
        end
    end
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
