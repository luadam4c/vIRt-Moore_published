function [intrinsicFreq, intrinsicAmp, intrinsicSetPt, frequencies, spectrum] = ...
    jm_virt_sim_calc_dominant_freq_spectral (effectorSteadyState, samplingPeriod, downsampleFactor, chronuxParams)
%% Calculates dominant frequency and amplitude of effector motion using spectral estimation
% Usage: [intrinsicFreq, intrinsicAmp, intrinsicSetPt, frequencies, spectrum] = ...
%           jm_virt_sim_calc_dominant_freq_spectral (effectorSteadyState, samplingPeriod, downsampleFactor, chronuxParams)
% Explanation:
%       Calculates the dominant frequency, signal amplitude, and set point
%       of the effector motion using multi-taper spectral estimation 
%       (Chronux).
%
% Outputs:
%       intrinsicFreq       - Dominant frequency (Hz)
%       intrinsicAmp        - Signal amplitude (derived from spectral power)
%       intrinsicSetPt      - Set point (mean of the signal)
%       frequencies         - Frequency vector from spectrum
%       spectrum            - Power spectrum
%
% Arguments:
%       effectorSteadyState - Vector of effector positions (steady state)
%                           (numeric vector) -- required
%       samplingPeriod      - Sampling period (s)
%                           (numeric scalar) -- required
%       downsampleFactor    - Downsampling factor
%                           (integer scalar) -- required
%       chronuxParams       - Structure of Chronux parameters
%                           (structure) -- required
%
% Requires:
%       mtspectrumc.m and dependent files from the Chronux toolbox
%
% Used by:
%       \Shared\Code\vIRt-Moore\jm_postprocess_virt_sim_amplitude_analysis.m
%       \Shared\Code\vIRt-Moore\jm_postprocess_virt_sim_amplitude_analysis_ExtCurrentRuns.m
%       \Shared\Code\vIRt-Moore\jm_postprocess_virt_sim_amplitude_analysis_reps.m
%       \Shared\Code\vIRt-Moore\jm_postprocess_virt_sim_ExtCurrent_analysis.m
%
% File History:
% 2025-12-02 Created by Jeff Moore
% 2026-01-13 Reorganized and reannotated by Gemini

%% Deal with arguments
if nargin < 4
    error('Not enough input arguments.');
end

%% Preparation
% Downsample the data
effectorDownsampled = effectorSteadyState(1:downsampleFactor:end);

% Remove mean for spectral calculation
intrinsicSetPt = mean(effectorDownsampled);
effectorCentered = effectorDownsampled - intrinsicSetPt;

% Calculate downsampled sampling frequency
fsDownsampled = (1 ./ samplingPeriod) / downsampleFactor;

% Update Chronux parameters with new Fs
chronuxParams.Fs = fsDownsampled;

%% Do the job
% 1. Compute Spectrum using Chronux
[spectrum, frequencies] = mtspectrumc(effectorCentered, chronuxParams);

% 2. Find Peak Frequency
% Find the index of the maximum power in the spectrum
peakFreqIndex = find(spectrum == max(spectrum), 1);
intrinsicFreq = frequencies(peakFreqIndex);

% 3. Calculate Signal Amplitude
% Calculate P-P signal amplitude approximation based on spectral power
% (Normalization tested with sinusoid of known amplitude)
intrinsicAmp = sqrt(sum(spectrum));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%