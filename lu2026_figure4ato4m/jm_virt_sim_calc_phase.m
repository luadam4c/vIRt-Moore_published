function [phaseResetIndices, phase, effectorDownsampled, fsDownsampled] = jm_virt_sim_calc_phase (effectorPosition, samplingPeriod, downsampleFactor, filterOrder, filterCutoffsBp)
%% Calculate Instantaneous Phase and Phase Resets
% Usage: [phaseResetIndices, phase, effectorDownsampled, fsDownsampled] = jm_virt_sim_calc_phase (effectorPosition, samplingPeriod, downsampleFactor, filterOrder, filterCutoffsBp)
% Explanation:
%       Downsamples the effector position data, filters it using a bandpass
%       Butterworth filter, and calculates the instantaneous phase using the
%       Hilbert transform. It also identifies phase reset times (transitions
%       from pi to -pi).
%
% Outputs:
%       phaseResetIndices   - Indices where the phase resets (diff < -pi)
%       phase               - Instantaneous phase (radians)
%       effectorDownsampled - Downsampled effector position data
%       fsDownsampled       - Downsampled sampling frequency (Hz)
%
% Arguments:
%       effectorPosition    - Vector of effector positions
%                           (numeric vector) -- required
%       samplingPeriod      - Sampling period of the original data (s)
%                           (numeric scalar) -- required
%       downsampleFactor    - Factor by which to downsample the data
%                           (integer scalar) -- required
%       filterOrder         - Order of the Butterworth filter
%                           (integer scalar) -- required
%       filterCutoffsBp     - Bandpass cutoffs for the filter [low, high] (Hz)
%                           (numeric vector) -- required
%
% Requires:
%       Signal Processing Toolbox (butter, filtfilt, hilbert)
%
% Used by:
%       \Shared\Code\vIRt-Moore\jm_virt_sim_calc_virt_spike_coherence.m
%
% File History:
% 2025-11-24 Created by Jeff Moore
% 2026-01-13 Reorganized and reannotated by Gemini

%% Deal with arguments
if nargin < 5
    error('Not enough input arguments.');
end

%% Preparation
% Calculate frequencies
fs = 1 / samplingPeriod;                 % Original frequency
fsDownsampled = fs / downsampleFactor;   % Downsampled frequency

% Downsample data
effectorDownsampled = effectorPosition(1:downsampleFactor:end); 

%% Do the job
% 1. Filter the data
% Calculate Nyquist-normalized cutoffs
nyquistFreq = fsDownsampled / 2;
normalizedCutoffs = filterCutoffsBp / nyquistFreq;

% Design Butterworth filter
[bCoeff, aCoeff] = butter(filterOrder, normalizedCutoffs); 

% Apply zero-phase filtering
filteredEffector = filtfilt(bCoeff, aCoeff, effectorDownsampled);

% 2. Calculate Instantaneous Phase
% Use Hilbert transform
% Note: Peaks will have phase = 0 and valleys will have phase = +/- pi
analyticSignal = hilbert(filteredEffector);
phase = angle(analyticSignal); 

% 3. Detect Phase Resets
% Detect where phase jumps from ~pi to ~-pi
phaseDiff = diff(phase);
phaseResetLogical = phaseDiff < -pi;
phaseResetIndices = find(phaseResetLogical);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%