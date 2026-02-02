function jm_annotate_imageplot_pixel_outlines (xVals, yVals, maskMatrix)
%% Plots color outlines on current axes for data in matrix M
% Usage: jm_annotate_imageplot_pixel_outlines (xVals, yVals, maskMatrix)
% Explanation:
%       Draws white rectangular outlines around pixels where maskMatrix is true (1).
%       This is used to highlight specific regions in a heatmap or image plot.
%
% Arguments:
%       xVals       - X-axis values corresponding to the columns of M
%                   (numeric vector) -- required
%       yVals       - Y-axis values corresponding to the rows of M
%                   (numeric vector) -- required
%       maskMatrix  - Binary matrix where 1 corresponds to pixels to outline 
%                     and 0 corresponds to no outline.
%                     (logical/numeric matrix) -- required
%
% Requires:
%
%
% Used by:
%       \Shared\Code\vIRt-Moore\jm_postprocess_virt_sim_ExtCurrent_analysis.m
%
% File History:
% 2025-12-08 Created by Jeff Moore
% 2026-01-13 Reorganized and reannotated by Gemini

%% Hard-coded parameters
lineColor = 'white';

%% Preparation
% Calculate half-width of pixels
dx = mean(diff(xVals)) / 2;
dy = mean(diff(yVals)) / 2;

%% Do the job
% Loop through x and y values
for iX = 1:length(xVals)
    for iY = 1:length(yVals)
        % If the mask is true for this pixel, draw the outline
        if maskMatrix(iY, iX)
            hold on;
            
            % Draw top edge
            line(xVals(iX) + [-dx  dx], yVals(iY) + [-dy -dy], 'Color', lineColor);
            
            % Draw bottom edge
            line(xVals(iX) + [-dx  dx], yVals(iY) + [ dy  dy], 'Color', lineColor);
            
            % Draw left edge
            line(xVals(iX) + [-dx -dx], yVals(iY) + [-dy  dy], 'Color', lineColor);
            
            % Draw right edge
            line(xVals(iX) + [ dx  dx], yVals(iY) + [-dy  dy], 'Color', lineColor);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%