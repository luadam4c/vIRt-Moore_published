function coherenceResults = jm_virt_sim_calc_virt_spike_coherence (effectorPosition, timeVector, spikeTimesRr, spikeTimesRp, chronuxParams, downsampleFactor, filterCutoffsBp, filterOrder)
%% Calculate Coherence Statistics for vIRt Simulation
% Usage: coherenceResults = jm_virt_sim_calc_virt_spike_coherence (effectorPosition, timeVector, spikeTimesRr, spikeTimesRp, chronuxParams, downsampleFactor, filterCutoffsBp, filterOrder)
% Explanation:
%       Calculates coherence statistics for simulated spiking activity 
%       (vIRt-Retraction and vIRt-Protraction pools) against the simulated 
%       whisker motion.
%
% Outputs:
%       coherenceResults    - Structure containing coherence analysis:
%           .Rr             - Stats for Retraction pool (C, phi, etc.)
%           .Rp             - Stats for Protraction pool (C, phi, etc.)
%           .phaseRef       - Reference phases (Retracted, Mid-protraction)
%           .effector       - Effector spectral stats
%           .Fsdwn          - Downsampled frequency
%           .timeSec        - Total simulation time
%           .ChronuxParams  - Parameters used
%
% Arguments:
%       effectorPosition    - Vector of effector position from simulation
%       timeVector          - Sample times corresponding to effectorPosition
%       spikeTimesRr        - Matrix of spike times (Time x Neuron) for Rr
%       spikeTimesRp        - Matrix of spike times (Time x Neuron) for Rp
%       chronuxParams       - Chronux parameters structure
%       downsampleFactor    - Integer downsampling factor
%       filterCutoffsBp     - Bandpass filter cutoffs (Hz)
%       filterOrder         - Filter order
%
% Requires:
%       mtspectrumc.m and dependent files from the Chronux toolbox
%       coherencycpt.m and dependent files from the Chronux toolbox
%       \Shared\Code\vIRt-Moore\jm_virt_sim_calc_phase.m
%
% Used by:
%       \Shared\Code\vIRt-Moore\jm_postprocess_virt_sim_spikeTrainStats.m
%
% File History:
% 2025-11-28 Created by Jeff Moore
% 2026-01-13 Reorganized and reannotated by Gemini
% 2026-01-17 Updated to use parfor for coherence loops by Gemini
% 2026-01-22 Modularized repetitive logic into subfunctions by Gemini

%% Deal with arguments
if nargin < 8
    error('Not enough input arguments.');
end

%% Preparation: Data Slicing
% Find the steady state indices (2nd half of simulation)
nSamples = length(effectorPosition);
steadyStateIndices = ceil(nSamples/2):nSamples;

% Restrict data to steady state (2nd half of simulation)
effectorSteadyState = effectorPosition(steadyStateIndices);
timeSteadyState = timeVector(steadyStateIndices);

% Downsample steady state data
effectorDownsampled = effectorSteadyState(1:downsampleFactor:end);
timeDownsampled = timeSteadyState(1:downsampleFactor:end);

% Compute downsampled sampling frequency
samplingInterval = mean(diff(timeVector)); 
fsDownsampled = (1/samplingInterval) / downsampleFactor;

% Remove spikes outside steady state window (replace with NaN)
tMin = timeSteadyState(1);
tMax = timeSteadyState(end);
spikeTimesRr(spikeTimesRr < tMin | spikeTimesRr > tMax) = NaN;
spikeTimesRp(spikeTimesRp < tMin | spikeTimesRp > tMax) = NaN;

% Shift time vector relative to steady state start so it starts at 0
timeDownsampledZeroed = timeDownsampled - tMin;

% Align spike times to steady state start
spikeTimesRr = spikeTimesRr - tMin;
spikeTimesRp = spikeTimesRp - tMin;

%% Do the job
% --- 1. Spectral Calculation of Effector ---
% Update sampling frequency
chronuxParams.Fs = fsDownsampled;

% Mean-substract the downsampled steady state effector trace
effectorZeroMean = effectorDownsampled - mean(effectorDownsampled);

% Compute the multi-taper spectrum of 
%   the mean-subtracted, downsampled steady state effector trace (but unfiltered)
[spectrumEffector, freqsEffector] = mtspectrumc(effectorZeroMean, chronuxParams);

% Find index of the peak (dominant) frequency in the multi-taper spectrum
idxDominantFreq = find(spectrumEffector == max(spectrumEffector), 1);
dominantFreq = freqsEffector(idxDominantFreq);

% --- 2. Phase Calculation ---
% Band-pass filter the downsampled effector trace at steady state,
%   then calculate the instantaneous phase using the Hilbert transform. 
% Note: protraction with phase -pi to 0 and retraction with phase 0 to pi
%   Also identifies phase reset indices (transitions from pi to -pi),
%       which are the protraction onset times
% Note: the indices are relative to the downsampled but non-zeroed time vector
[phaseResetIndices, instantaneousPhase, ~, ~] = ...
    jm_virt_sim_calc_phase(effectorSteadyState, samplingInterval, ...
                            downsampleFactor, filterOrder, filterCutoffsBp);

% --- 3. Coherence for Rp (Protraction) Units ---
disp('Calculating Rp coherence values (Parallel)...')

% Calculate the coherency of each Rp neuron with the steady state effector trace
%   and the phases of the spikes relative to the protraction onset of the steady state effector trace
[cohereRp, phiRp, confCRp, nSpikesRp, spikePhaseRp] = ...
    calc_coherence_by_neuron(spikeTimesRp, effectorZeroMean, chronuxParams, ...
                              idxDominantFreq, timeDownsampledZeroed, instantaneousPhase);


% --- 4. Coherence for Rr (Retraction) Units ---
disp('Calculating Rr coherence values (Parallel)...')

% Calculate the coherency of each Rr neuron with the steady state effector trace
%   and the phases of the spikes relative to the protraction onset of the steady state effector trace
[cohereRr, phiRr, confCRr, nSpikesRr, spikePhaseRr] = ...
    calc_coherence_by_neuron(spikeTimesRr, effectorZeroMean, chronuxParams, ...
                              idxDominantFreq, timeDownsampledZeroed, instantaneousPhase);

% --- 5. Aggregate Coherence (Network Level) ---
% Flatten spike matrices to single vectors and only retain valid spikes
allSpikesRp = spikeTimesRp(:);
allSpikesRp = allSpikesRp(~isnan(allSpikesRp));

allSpikesRr = spikeTimesRr(:);
allSpikesRr = allSpikesRr(~isnan(allSpikesRr));

% Compute Aggregate Coherency of Rp neurons with the steady state effector trace
[C, phi, ~, ~, ~, ~, ~, confC] = coherencycpt(effectorZeroMean, allSpikesRp, chronuxParams);
cohereRpAgg = C(idxDominantFreq); 
phiRpAgg = -phi(idxDominantFreq); 
confCRpAgg = confC;

% Compute Aggregate Coherency of Rr neurons with the steady state effector trace
[C, phi, ~, ~, ~, ~, ~, confC] = coherencycpt(effectorZeroMean, allSpikesRr, chronuxParams);
cohereRrAgg = C(idxDominantFreq); 
phiRrAgg = -phi(idxDominantFreq); 
confCRrAgg = confC;

% Flatten spike phases
allSpikesRpPhase = spikePhaseRp(:);
allSpikesRpPhase = allSpikesRpPhase(~isnan(allSpikesRpPhase));

allSpikesRrPhase = spikePhaseRr(:);
allSpikesRrPhase = allSpikesRrPhase(~isnan(allSpikesRrPhase));

% --- 6. Phase Reference Calculation ---
% Determine protraction onset times relative to the zeroed downsampled time vector
protractOnsetTimes = timeDownsampled(phaseResetIndices) - tMin;

% Calculate coherency of the protraction onset times with the steady state effector trace
[C, phi] = coherencycpt(effectorZeroMean, protractOnsetTimes, chronuxParams);
cohereProtractOnset = C(idxDominantFreq);
phiProtractOnset = -phi(idxDominantFreq);

% Compute quarter cycle times (mid-protraction) relative to the zeroed downsampled time vector
dominantPeriod = 1 / dominantFreq;
dominantQuarterCycle = dominantPeriod / 4;
quarterCycleTime = protractOnsetTimes + dominantQuarterCycle;

% Calculate coherency of the quarter cycle times with the steady state effector trace
[C, phi] = coherencycpt(effectorZeroMean, quarterCycleTime, chronuxParams);
cohereQuarterCycle = C(idxDominantFreq);
phiQuarterCycle = -phi(idxDominantFreq);

% --- 7. ISI and CV2 Calculation ---

% Calculate interspike intervals and CV2 for Rr
[isiRrAgg, cv2Rr] = calc_isi_stats(spikeTimesRr, 2 * dominantQuarterCycle);

% Calculate interspike intervals and CV2 for Rp
[isiRpAgg, cv2Rp] = calc_isi_stats(spikeTimesRp, 2 * dominantQuarterCycle);

%% Output results
% Collect output structure
% Note: Define protraction from phase 0 to pi instead of -pi to 0
%               retraction from phase pi to 2pi instead of 0 to pi
%       (add pi to phase outputs)

% Rr Data
coherenceResults.Rr.C = cohereRr;
coherenceResults.Rr.phi = phiRr + pi;
coherenceResults.Rr.confC = confCRr;
coherenceResults.Rr.CAgg = cohereRrAgg;
coherenceResults.Rr.phiAgg = phiRrAgg + pi;
coherenceResults.Rr.confCAgg = confCRrAgg;
coherenceResults.Rr.spikePhase = allSpikesRrPhase + pi;
coherenceResults.Rr.Nspikes = nSpikesRr;
coherenceResults.Rr.isiAgg = isiRrAgg;
coherenceResults.Rr.cv2 = cv2Rr;

% Rp Data
coherenceResults.Rp.C = cohereRp;
coherenceResults.Rp.phi = phiRp + pi;
coherenceResults.Rp.confC = confCRp;
coherenceResults.Rp.CAgg = cohereRpAgg;
coherenceResults.Rp.phiAgg = phiRpAgg + pi;
coherenceResults.Rp.confCAgg = confCRpAgg;
coherenceResults.Rp.spikePhase = allSpikesRpPhase + pi;
coherenceResults.Rp.Nspikes = nSpikesRp;
coherenceResults.Rp.isiAgg = isiRpAgg;
coherenceResults.Rp.cv2 = cv2Rp;

% Phase References
coherenceResults.phaseRef.phiProtractOnset = phiProtractOnset + pi;
coherenceResults.phaseRef.CProtractOnset = cohereProtractOnset;
coherenceResults.phaseRef.phiQuarterCycle = phiQuarterCycle + pi;
coherenceResults.phaseRef.CQuarterCycle = cohereQuarterCycle;

% Effector Stats
coherenceResults.effector.dominantFreq = dominantFreq;
coherenceResults.effector.allPhase = instantaneousPhase + pi;
coherenceResults.effector.S = spectrumEffector;
coherenceResults.effector.f = freqsEffector;
coherenceResults.filtCutoffsBP = filterCutoffsBp;
coherenceResults.filtOrder = filterOrder;

% General Stats
coherenceResults.Fsdwn = fsDownsampled;
coherenceResults.timeSec = timeDownsampledZeroed(end);
coherenceResults.ChronuxParams = chronuxParams;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cohere, phi, confC, nSpikes, spikePhase] = ...
    calc_coherence_by_neuron(spikeTimes, effectorZeroMean, chronuxParams, ...
                    idxDominantFreq, timeDownsampledZeroed, instantaneousPhase)
% Calculates coherence for a matrix of spike times against the effector.
% Uses parfor for parallel processing of neurons.

% Count the number of neurons
nNeurons = size(spikeTimes, 2);

% Initialize arrays
cohere = zeros(1, nNeurons); 
phi = zeros(1, nNeurons); 
confC = zeros(1, nNeurons); 
nSpikes = zeros(1, nNeurons);

% Use cell array for parfor variable length output
spikePhaseCell = cell(1, nNeurons); 

parfor iNeuron = 1:nNeurons
    % Extract valid spikes for this neuron
    spikesThisNeuron = spikeTimes(~isnan(spikeTimes(:, iNeuron)), iNeuron);

    % Store the number of spikes for this neuron
    nSpikes(iNeuron) = length(spikesThisNeuron);
    
    % Calculate Coherency (Continuous vs Point Process)
    [C, phiOut, ~, ~, ~, ~, ~, confCOut] = ...
        coherencycpt(effectorZeroMean, spikesThisNeuron, chronuxParams);
    
    % Store coherency, its phase, and its confidence level at the dominant frequency
    % Note: phase returned by coherencycpt() is sign-inverted due to 
    %       the computation of the cross-spectrum having the -i*omega*t term
    cohere(iNeuron) = C(idxDominantFreq);
    phi(iNeuron) = -phiOut(idxDominantFreq);
    confC(iNeuron) = confCOut;
    
    % Compute the phases of the spike times by interpolation
    if ~isempty(spikesThisNeuron)
        spikePhaseCell{iNeuron} = ...
            interp1(timeDownsampledZeroed, instantaneousPhase, spikesThisNeuron);
    else
        spikePhaseCell{iNeuron} = [];
    end
end

% Reconstruct matrix from cell array
spikePhase = NaN(size(spikeTimes));
for iNeuron = 1:nNeurons
    phases = spikePhaseCell{iNeuron};
    if ~isempty(phases)
        spikePhase(1:length(phases), iNeuron) = phases;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [isiAgg, cv2] = calc_isi_stats(spikeTimes, dominantHalfCycle)
% Calculates aggregated ISI and CV2 (local variation) for a spike matrix.

% Compute ISIs
isiMatrix = diff(spikeTimes);

% Aggregated ISI
isiAll = isiMatrix(:);
isiAgg = isiAll(~isnan(isiAll));

% Exclude long pauses greater than half the dominant period
isiCycle = isiMatrix;
isiCycle(isiMatrix > dominantHalfCycle) = NaN;

% Calculate CV2 (local variation)
isi1 = isiCycle(1:(end-1), :);
isi2 = isiCycle(2:end, :);
cv2 = nanmean(2 * abs(isi2 - isi1) ./ (isi2 + isi1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%