function spikeWhiskCohereStats = jm_virt_sim_plot_virt_spike_coherence (coherenceResults)
%% Plot vIRt Simulation Spike Coherence
% Usage: spikeWhiskCohereStats = jm_virt_sim_plot_virt_spike_coherence (coherenceResults)
% Explanation:
%       Plots coherence statistics between simulated unit spiking and effector
%       motion. Includes power spectra, polar plots of coherence, phase
%       histograms, firing rate distributions, and ISI distributions.
%
% Outputs:
%       spikeWhiskCohereStats - Unit statistics (currently empty/undefined in original code logic)
%
% Arguments:
%       coherenceResults      - Structure created by jm_virt_sim_calc_virt_spike_coherence
%                               containing fields Rr, Rp, phaseRef, effector, etc.
%
% Requires:
%       Matlab Plotting Toolbox
%
% Used by:
%       \Shared\Code\vIRt-Moore\jm_postprocess_virt_sim_spikeTrainStats.m
%
% File History:
% 2025-12-01 Created by Jeff Moore
% 2026-01-13 Reorganized and reannotated by Gemini
% 2026-01-22 Fixed bug in rate calculations

%% Hard-coded parameters
nBinsPhase = 16;        % Number of phase bins
nBinsCoherence = 16;    
nBinsFiringRate = 16;
nBinsIsi = 20;
nBinsCv2 = 20;

%% Deal with arguments
if nargin < 1
    error('Need coherenceResults structure');
end

%% Preparation
% Extract Data for readability
effFreqs = coherenceResults.effector.f;
effSpectrum = coherenceResults.effector.S;
allPhase = coherenceResults.effector.allPhase;

% Rr Data
rrPhi = coherenceResults.Rr.phi;
rrC = coherenceResults.Rr.C;
rrConf = coherenceResults.Rr.confC;
rrSpikePhase = coherenceResults.Rr.spikePhase;
rrNSpikes = coherenceResults.Rr.Nspikes;
rrIsiAgg = coherenceResults.Rr.isiAgg;
rrCv2 = coherenceResults.Rr.cv2;

% Rp Data
rpPhi = coherenceResults.Rp.phi;
rpC = coherenceResults.Rp.C;
rpConf = coherenceResults.Rp.confC;
rpSpikePhase = coherenceResults.Rp.spikePhase;
rpNSpikes = coherenceResults.Rp.Nspikes;
rpIsiAgg = coherenceResults.Rp.isiAgg;
rpCv2 = coherenceResults.Rp.cv2;

% Reference Data
refPhiProtractOnset = coherenceResults.phaseRef.phiProtractOnset;
refCProtractOnset = coherenceResults.phaseRef.CProtractOnset;
refPhiQuarter = coherenceResults.phaseRef.phiQuarterCycle;
refCQuarter = coherenceResults.phaseRef.CQuarterCycle;

%% Do the job
% Calculate histogram bins for phase
[nSamplesEachPhase, binCentersPhase] = hist(allPhase, nBinsPhase);

% Calculate histogram bins for coherence (not used)
[~, binCentersCoherence] = hist(0:.01:1, nBinsCoherence); 

% Calculate histogram bins for mean spike rate
[~, binCentersFr] = hist(10:.1:25, nBinsFiringRate);

% Calculate histogram bins for ISIs (not used)
[~, binCentersIsi] = hist(0:1:20, nBinsIsi);

% Calculate histogram bins for CV2
[~, binCentersCv2] = hist(0:0.01:2, nBinsCv2);

% If data changed, use the following for full range:
% [~, binCentersCoherence] = hist(0:1/nBinsCoherence:1, nBinsCoherence); % Full 0-1 range
% [~, binCentersFr] = hist(rrNSpikes / coherenceResults.timeSec, nBinsFiringRate); % Dynamic based on data
% [~, binCentersIsi] = hist(rrIsiAgg * 1000, nBinsIsi); % Dynamic based on data
% [~, binCentersCv2] = hist(rrCv2, nBinsCv2); % Dynamic based on data

% Calculate Distributions
% 1. Phase histograms (tuning curves)
% nRrPhi = hist(rrPhi, binCentersPhase);
% nRpPhi = hist(rpPhi, binCentersPhase);

nRrSpikePhase = hist(rrSpikePhase, binCentersPhase);
nRpSpikePhase = hist(rpSpikePhase, binCentersPhase);

% 2. Compute Phase Rates normalize by time spent in bin
nNeuronsRr = length(rrNSpikes);
nNeuronsRp = length(rpNSpikes);

% Calculate occupancy (time spent in each phase bin)
timeEachPhase = nSamplesEachPhase / coherenceResults.Fsdwn;

% Mean Spike Rate = Spikes / Time (Hz) by Phase
ratePhaseRr = (nRrSpikePhase ./ nNeuronsRr) ./ timeEachPhase;
ratePhaseRp = (nRpSpikePhase ./ nNeuronsRp) ./ timeEachPhase;

% 3. Coherence Histograms
% nRrC = hist(rrC, binCentersCoherence);
% nRpC = hist(rpC, binCentersCoherence);

% 4. Mean Spike Rate Histograms
histFrRr = hist(rrNSpikes / coherenceResults.timeSec, binCentersFr);
histFrRp = hist(rpNSpikes / coherenceResults.timeSec, binCentersFr);

% 5. ISI Histograms
% nIsiRr = hist(rrIsiAgg * 1000, binCentersIsi);
% nIsiRp = hist(rpIsiAgg * 1000, binCentersIsi);

% 6. CV2 Histograms
histCv2Rr = hist(rrCv2, binCentersCv2);
histCv2Rp = hist(rpCv2, binCentersCv2);

%% Plotting
figure;

% Subplot 1: Effector Spectrum
subplot(1, 9, 3);
diffF = mean(diff(effFreqs));
plot(effFreqs, effSpectrum / sum(effSpectrum) / diffF);
xlabel('Frequency (Hz)');
ylabel('Fraction of Signal Power');

% Subplot 2: Polar Plot of Coherence
subplot(1, 9, [4, 5]);

% Plot Rr Units (Black), but only fill in the significant units
sigIdxRr = rrC > rrConf;
polarplot(rrPhi(sigIdxRr), rrC(sigIdxRr), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 2);
hold on;
polarplot(rrPhi(~sigIdxRr), rrC(~sigIdxRr), 'ko', 'MarkerSize', 2);

% Plot Rp Units (Green), but only fill in the significant units
sigIdxRp = rpC > rpConf;
polarplot(rpPhi(sigIdxRp), rpC(sigIdxRp), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 2);
polarplot(rpPhi(~sigIdxRp), rpC(~sigIdxRp), 'go', 'MarkerSize', 2);

% Plot Phase References
polarplot(refPhiProtractOnset, refCProtractOnset, 'b.'); 
polarplot(refPhiQuarter, refCQuarter, 'r.'); 
text(refPhiProtractOnset, refCProtractOnset, ' Protraction-Onset'); 
text(refPhiQuarter, refCQuarter, ' Mid-Protraction'); 

% Plot Aggregate Coherence
polarplot(coherenceResults.Rr.phiAgg, coherenceResults.Rr.CAgg, 'kx');
polarplot(coherenceResults.Rp.phiAgg, coherenceResults.Rp.CAgg, 'gx');
text(coherenceResults.Rr.phiAgg, coherenceResults.Rr.CAgg, ' Aggregate Rr spikes');
text(coherenceResults.Rp.phiAgg, coherenceResults.Rp.CAgg, 'Aggregate Rp spikes');

hold off;

% Subplot 3: Phase Rate Tuning Curves (Polar Line)
subplot(1, 9, [6, 7]);
polarplot([binCentersPhase, binCentersPhase(1)], [ratePhaseRr, ratePhaseRr(1)], 'k.-');
hold on;
polarplot([binCentersPhase, binCentersPhase(1)], [ratePhaseRp, ratePhaseRp(1)], 'g.-');
hold off;

% Subplot 4: Firing Rate Distribution
subplot(1, 9, 8); 
hold on;
bar(binCentersFr, histFrRr / sum(histFrRr), 1, 'FaceColor', 'none', 'EdgeColor', 'k');
bar(binCentersFr, histFrRp / sum(histFrRp), 1, 'FaceColor', 'none', 'EdgeColor', 'g');
xlabel('Mean spike rate (Hz)');
ylabel('Fraction of neurons');
hold off;

% Subplot 5: CV2 Distribution
subplot(1, 9, 9); 
hold on;
bar(binCentersCv2, histCv2Rr / sum(histCv2Rr), 1, 'FaceColor', 'none', 'EdgeColor', 'k');
bar(binCentersCv2, histCv2Rp / sum(histCv2Rp), 1, 'FaceColor', 'none', 'EdgeColor', 'g');
xlabel('Mean CV2');
ylabel('Fraction of neurons');
hold off;

% Set output (Placeholder as logic was not defined in original)
spikeWhiskCohereStats = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%