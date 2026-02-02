function [r, p, ci95] = jm_virt_sim_corr_bootci(x, y, numIter)
%% Computes the 95% confidence interval using a bootstrap method
% Usage: [r, p, ci95] = jm_virt_sim_corr_bootci(x, y, numIter)
% Explanation:
%       Calculates the Pearson correlation coefficient (r) between two column 
%       vectors x and y.
%       Computes the p-value using a shuffle test (permutation).
%       Computes the 95% confidence interval using a bootstrap method.
%
% Inputs:
%       x           - Column vector of data (Double)
%       y           - Column vector of data (Double)
%       numIter     - (Optional) Number of iterations for bootstrap/shuffle 
%                     (Default: 10000)
%
% Outputs:
%       r           - Pearson correlation coefficient
%       p           - p-value derived from shuffle test
%       ci95        - 95% confidence interval [lower, upper] from bootstrap
%
% Requires:
%
% Used by:
%       \Shared\Code\vIRt-Moore\jm_postprocess_virt_sim_ExtCurrent_analysis.m
%
% File History:
% 2026-01-26 Created by Jeff Moore
% 2026-01-26 Annotated by Gemini

if nargin < 3
    numIter = 10000;
end

% Ensure x, y are column vars
% (Logic assumes alignment, input validation usually handled by caller)
N = size(y, 1);

% Calculate actual correlation
r = corr(x, y);

% --- Shuffle Test (Permutation) ---
nullRs = zeros(numIter, 1);
for i = 1:numIter
    % Shuffle ONLY the labels of y
    shuffledY = y(randperm(N));

    % Calculate correlation on the shuffled data
    nullRs(i) = corr(x, shuffledY);
end

% Calculate the p-value: what fraction of shuffles beat our real correlation?
p = mean(abs(nullRs) >= abs(r));

% --- Bootstrap Confidence Interval ---
bootStat = bootstrp(numIter, @corr, x, y);
ci95 = prctile(bootStat, [2.5, 97.5]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%