function [handles, P, simData] = run_rlc_simulations (varargin)
%% Simulate parallel RLC circuit driven by sinusoidal current source
% Usage: [handles, P, simData] = run_rlc_simulations (varargin)
% Explanation:
%       This function simulates a parallel RLC circuit consisting of a 
%       current source, resistor (R), inductor (L), and capacitor (C).
%
%       The state-space equations are derived from Kirchhoff's Current Law 
%       (KCL) at the top node.
%
%       It performs two types of analyses:
%       1. Frequency Sweep: Varying input frequency with fixed C.
%       2. Capacitance Sweep: Varying C (and thus natural frequency) with 
%          fixed input frequency.
%
% Example(s):
%       run_rlc_simulations;
%       run_rlc_simulations('SaveOutput', false);
%       [h, P, data] = run_rlc_simulations('InputAmp', 2, 'Resistor', 20);
%
% Outputs:
%       handles     - handles to all created graphics objects.
%                   specified as a structure
%       P           - parameter structure used for the simulation.
%                   specified as a structure
%       simData     - gathered data from the sweeps.
%                   specified as a structure with fields:
%                       .freqSweep: results from frequency sweep
%                       .capSweep: results from capacitance sweep
%
% Arguments:
%       varargin    - 'OutDir': output directory
%                   must be a string scalar or a character vector
%                   default == fullfile(pwd, sprintf('RLC_Test-%s', timeStamp))
%                   - 'InputAmp': Amplitude of current source (Amps)
%                   must be a numeric scalar
%                   default == 1
%                   - 'PulseDur': Pulse duration (seconds)
%                   must be a numeric scalar
%                   default == 0.05
%                   - 'Resistor': Resistance (Ohms)
%                   must be a numeric scalar
%                   default == 10.0
%                   - 'CapacitanceRef': Reference Capacitance (Farads)
%                   must be a numeric scalar
%                   default == 10e-3
%                   - 'NaturalFreqRef': Target natural frequency for L calc (Hz)
%                   must be a numeric scalar
%                   default == 6.0
%                   - 'SaveOutput': Whether to save figures.
%                   must be a logical scalar
%                   default == true
%                   - 'ShowFigure': Whether to show figures.
%                   must be a logical scalar
%                   default == true
%
% Requires:
%       (Standard MATLAB Toolboxes)
%

% File History:
% 2025-11-20 rlc_sim_genfigs_C.m created by Dr. Jeff Moore with the help of Gemini
% 2026-01-09 Adapted from rlc_sim_genfigs_C.m to function format
% 2025-01-09 Enforced camelCase for variables and snake_case for functions

%% Compatibility Check
% Check MATLAB version for function flag support
v = ver('matlab');
versionStr = regexp(v.Version, '^\d+\.\d+', 'match', 'once');
matlabVersion = str2double(versionStr);

%% Hard-Coded Information For This Script Only
% timestamp generator
timeStamp = datestr(now, 'yyyymmdd_HHMMSS'); 
nameModel = 'RLC_Sim';

%% Default values for optional arguments
dirOutDefault = ''; 
inputAmpDefault = 1;
pulseDurDefault = 0.05;
resistorDefault = 10.0;
capRefDefault = 10e-3;
natFreqRefDefault = 6.0;
toSaveOutputDefault = true;
showFigureDefault = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Input Parser implementation
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;

addParameter(iP, 'OutDir', dirOutDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'InputAmp', inputAmpDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'PulseDur', pulseDurDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'Resistor', resistorDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'CapacitanceRef', capRefDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'NaturalFreqRef', natFreqRefDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'SaveOutput', toSaveOutputDefault, ...
    @(x) validateattributes(x, {'logical'}, {'scalar'}));
addParameter(iP, 'ShowFigure', showFigureDefault, ...
    @(x) validateattributes(x, {'logical'}, {'scalar'}));

parse(iP, varargin{:});

% Extract parameters into camelCase variables
dirOut = iP.Results.OutDir;
iAmp = iP.Results.InputAmp;
pulseDur = iP.Results.PulseDur;
resistanceVal = iP.Results.Resistor;
capRef = iP.Results.CapacitanceRef;
natFreqRef = iP.Results.NaturalFreqRef;
toSaveOutput = iP.Results.SaveOutput;
showFigure = iP.Results.ShowFigure;

%% I. INITIALIZATION AND SETUP
% =========================================================================

% Handle output directory
if isempty(dirOut)
    dirOut = fullfile(pwd, sprintf('%s_Test-%s', nameModel, timeStamp));
end

if toSaveOutput
    if ~exist(dirOut, 'dir')
        mkdir(dirOut);
    end
    
    % Initialize diary
    diaryFile = fullfile(dirOut, 'simulation_log.txt');
    if exist(diaryFile, 'file')
        delete(diaryFile);
    end
    diary(diaryFile);
end

% Calculate Inductance based on reference values
% w_n^2 = 1 / (LC) => L = 1 / (w_n^2 * C)
wNRef = 2 * pi * natFreqRef;
inductanceVal = 1 / (wNRef^2 * capRef);

% Print calculations to standard output
fprintf('--- Circuit Parameter Calculations ---\n');
fprintf('Target Natural Frequency: %.2f Hz\n', natFreqRef);
fprintf('Reference Capacitance: %.2e F\n', capRef);
fprintf('Calculated Angular Frequency (w_n): %.4f rad/s\n', wNRef);
fprintf('Calculated Inductance (L): %.6e H\n', inductanceVal);
fprintf('--------------------------------------\n');

% Store parameters in P structure for output
P.iAmp = iAmp;
P.pulseDur = pulseDur;
P.resistanceVal = resistanceVal;
P.capRef = capRef;
P.natFreqRef = natFreqRef;
P.inductanceVal = inductanceVal;

handles = struct();
simData = struct();

%% II. FREQUENCY SWEEP SIMULATION
% =========================================================================

% Refactored Frequency Sweep
fprintf('Starting Frequency Sweep...\n');

fInSweepVector = 1:0.1:10;
cStatic = capRef; % Use the parameter passed in
plotTraceSweep = (fInSweepVector - floor(fInSweepVector) == 0); % Plot integer frequencies

% Pre-allocate
vSwing = zeros(size(fInSweepVector));
fNOut = zeros(size(fInSweepVector));

for i = 1:length(fInSweepVector)
    % Call the local simulation function
    % Note: passing 'toSaveOutput' to the subfunction controls individual trace saving
    [vSwing(i), fNOut(i)] = rlc_simulation_core(iAmp, fInSweepVector(i), pulseDur, ...
                                            cStatic, inductanceVal, resistanceVal, ...
                                            plotTraceSweep(i), toSaveOutput, ...
                                            dirOut, showFigure);
end

simData.freqSweep.fIn = fInSweepVector;
simData.freqSweep.vSwing = vSwing;
simData.freqSweep.fN = fNOut;

% Plotting Frequency Response
if showFigure
    handles.figFreqResp = figure('Visible', 'on');
else
    handles.figFreqResp = figure('Visible', 'off');
end

plot(fInSweepVector, vSwing);
xlabel('Drive frequency (Hz)');
ylabel('Vc amplitude (V)');
title('RLC circuit response to current pulse trains')
line([1 1]*fNOut(1), ylim, 'Color', 'r', 'LineStyle', '--');
grid on;

% Saving
if toSaveOutput
    saveas(handles.figFreqResp, fullfile(dirOut, 'Frequency_response_fin.fig'));
    saveas(handles.figFreqResp, fullfile(dirOut, 'Frequency_response_fin.png')); 
end


%% III. CAPACITANCE SWEEP SIMULATION
% =========================================================================

% Refactored Capacitance Sweep
fprintf('Starting Capacitance Sweep...\n');

cTuneVector = logspace(log10(4), log10(50), 100) * 1e-3;
fInStatic = 6; % Fixed frequency for this test
plotTraceTune = zeros(size(cTuneVector));
plotTraceTune(1:10:end) = 1; % Plot every 10th
saveIndividualPlots = false; % Override saving for these traces to avoid clutter

% Pre-allocate
vSwingC = zeros(size(cTuneVector));
fNC = zeros(size(cTuneVector));

for i = 1:length(cTuneVector)
    % Call local simulation function
    % Recalculate Inductance is NOT done here; L is fixed from Section I.
    % Only C changes, which changes the natural frequency of the circuit.
    [vSwingC(i), fNC(i)] = rlc_simulation_core(iAmp, fInStatic, pulseDur, ...
                                             cTuneVector(i), inductanceVal, resistanceVal, ...
                                             plotTraceTune(i), saveIndividualPlots, ...
                                             dirOut, showFigure);
end

simData.capSweep.cTune = cTuneVector;
simData.capSweep.vSwing = vSwingC;
simData.capSweep.fN = fNC;

% Plotting Capacitance Sweep Response
if showFigure
    handles.figCapResp = figure('Visible', 'on');
else
    handles.figCapResp = figure('Visible', 'off');
end

plot(fNC, vSwingC);
xlabel('RLC Intrinsic Frequency (Hz)');
ylabel('Vc amplitude (V)');
title('RLC circuit response to current pulse trains (Varying C)')
line([1 1]*fInStatic, ylim, 'Color', 'r', 'LineStyle', '--');

% Add labels
cTunePlot = round(cTuneVector(logical(plotTraceTune))*1e3, 0); %mF
text(fNC(logical(plotTraceTune)), vSwingC(logical(plotTraceTune)), num2cell(cTunePlot));
hold on
plot(fNC(logical(plotTraceTune)), vSwingC(logical(plotTraceTune)), 'r.');
xlim([1 10]);
grid on;

% Saving
if toSaveOutput
    saveas(handles.figCapResp, fullfile(dirOut, 'Frequency_response_C.fig'));
    saveas(handles.figCapResp, fullfile(dirOut, 'Frequency_response_C.png'));
    
    % Save P and simData structs to .mat file
    saveFileName = fullfile(dirOut, 'simulation_data.mat');
    save(saveFileName, 'P', 'simData');
    fprintf('Saved simulation parameters and data to: %s\n', saveFileName);
end

fprintf('Simulations Complete.\n');

if toSaveOutput
    diary off;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Local Functions

function [vSwing, fNCalc] = rlc_simulation_core(iAmp, fIn, pulseDur, capacitanceVal, inductanceVal, resistanceVal, plotTrace, toSavePlots, dirOut, showFigure)
% Note: This function is adapted from the original `rlc_simulation` but now accepts
% L and R as arguments rather than calculating them internally or hardcoding them.

    % --- Calculate Pulse Parameters ---
    period = 1 / fIn; % Period of the square wave
    timeOn = pulseDur; % Time the pulse is "on"
    
    % New Code:
    wNCalc = 1/sqrt(inductanceVal * capacitanceVal);
    fNCalc = wNCalc/(2*pi);

    % --- Simulation Parameters ---
    tEnd = 3.0; % Simulate for 3 seconds
    tSpan = [0 tEnd];
    
    % Initial conditions: [v(0), I_L(0)]
    x0 = [0; 0]; % Start with zero voltage and zero inductor current

    % --- Define the ODE System ---
    % x is a vector [x1; x2] = [v; I_L]
    % t is time
    
    % Define the input current function (square wave)
    % iInFunc = iAmp if mod(t, period) < timeOn, else 0
    iInFunc = @(t) iAmp * (mod(t, period) < timeOn);
    
    odeSystem = @(t, x) [ (iInFunc(t) - x(1)/resistanceVal - x(2)) / capacitanceVal; ...
                           x(1) / inductanceVal ];

    % --- Run the Simulation ---
    % Use ode45, a standard ODE solver in MATLAB
    [tVec, xMat] = ode45(odeSystem, tSpan, x0);

    % --- Extract Results ---
    vTime = xMat(:, 1);     % Voltage across components (x1)
    iLTime = xMat(:, 2);    % Current through inductor (x2)
    
    % Calculate other currents for plotting
    iInTime = iInFunc(tVec);            % Input current (using the same function)
    iRTime = vTime / resistanceVal;     % Current through resistor
    
    % Use gradient or diff for Capacitor current
    iCTime = capacitanceVal * gradient(vTime, tVec);      % Current through capacitor (numerical derivative)
    
    % Verify KCL (for checking): I_total = I_R + I_L + I_C
    iTotalTime = iRTime + iLTime + iCTime;

    % --- Plot the Results for integer values of input frequency---
    if plotTrace
        if showFigure
            figTrace = figure('Visible', 'on');
        else
            figTrace = figure('Visible', 'off');
        end

        % Plot 1: Voltage
        subplot(2, 1, 1);
        plot(tVec, vTime, 'b-', 'LineWidth', 1.5);
        title(['Parallel RLC Circuit Simulation (f_{in} = ', num2str(fIn), ' Hz; f_n = ',num2str(fNCalc),')']);
        ylabel('Voltage (V)');
        xlabel('Time (s)');
        grid on;
        legend('Voltage (v)');
        ylim([-10 10])

        % Plot 2: Currents
        subplot(2, 1, 2);
        plot(tVec, iInTime, 'k-', 'LineWidth', 2.0);
        hold on;
        plot(tVec, iRTime, 'r--');
        plot(tVec, iLTime, 'g--');
        plot(tVec, iCTime, 'm--');
        % Plot KCL check
        % plot(tVec, iTotalTime, 'c:', 'LineWidth', 2.0);
        hold off;
        ylabel('Current (A)');
        xlabel('Time (s)');
        grid on;
        legend('I_{source}', 'I_R', 'I_L', 'I_C');
        % legend('I_{source}', 'I_R', 'I_L', 'I_C', 'I_R+I_L+I_C');
        ylim([-5 5])
        
        if toSavePlots
            % Use fullfile and dirOut
            fname = fullfile(dirOut, ['Circuit_response_fin=', num2str(fIn)]);
            saveas(figTrace, [fname, '.fig']);
            saveas(figTrace, [fname, '.png']);
        end
        
        % Close the trace figure to prevent memory buildup during loops
        if ~showFigure || toSavePlots
            % Optional: close(figTrace); 
        end
    end

    % Calculate maximum voltage swing at steady state
    % ensure we have enough data points
    if length(vTime) > 10
        vSwing = range(vTime(round(end/2):end));
    else
        vSwing = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%