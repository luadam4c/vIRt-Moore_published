function jm_annotate_imageplot_pixel_values_text (xVals, yVals, maskMatrix)
%% Displays text of values in matrix M on corresponding pixels in an image plot
% Usage: jm_annotate_imageplot_pixel_values_text (xVals, yVals, maskMatrix)
% Explanation:
%       Iterates through the provided x and y values and places text 
%       annotations on the current plot corresponding to the values in 
%       the matrix M.
%
% Arguments:
%       xVals       - X-axis values corresponding to the columns of M
%                   (numeric vector) -- required
%       yVals       - Y-axis values corresponding to the rows of M
%                   (numeric vector) -- required
%       maskMatrix  - Matrix containing values to display. 
%                     (numeric matrix) -- required
%
% Requires:
%
%
% Used by:
%       \Shared\Code\vIRt-Moore\jm_postprocess_virt_sim_ExtCurrent_analysis.m
%
% File History:
% 2025-12-06 Created by Jeff Moore
% 2026-01-13 Reorganized and reannotated by Gemini

%% Hard-coded parameters
precision = 3;      % Significant figures for text display
textColor = 'white';

%% Do the job
% Loop through x and y values
for iX = 1:length(xVals)
    for iY = 1:length(yVals)
        % Convert value to string with specified precision
        textToPlot = num2str(maskMatrix(iY, iX), precision);
        
        % Place text on the plot
        text(xVals(iX), yVals(iY), textToPlot, 'Color', textColor);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%