function newName = virt_updated_test_name_for_rep (oldName, cleanParamName, strParamVaried)
%% Updates a test name by replacing or appending a parameter-value string.
% Usage: newName = virt_updated_test_name_for_rep (oldName, cleanParamName, strParamVaried)
%
% Explanation:
%       This function updates a test name string by replacing an existing
%       parameter-value pair or appending a new one. It relies on the specific
%       hyphen-separated naming convention used in virt_moore simulations.
%       It also ensures 'seedNumber' and 'nReps' tags are removed or managed
%       appropriately to keep the filename clean.
%
% Outputs:
%       newName         - The updated test name string.
%                       specified as a character vector
%
% Arguments:
%       oldName         - The original test name string.
%                       must be a string scalar or a character vector
%       cleanParamName  - The simplified name of the parameter being varied
%                         (e.g., 'Pb-period1').
%                       must be a string scalar or a character vector
%       strParamVaried  - The formatted string containing the parameter and
%                         its specific value (e.g., 'Pb-period1-200').
%                       must be a string scalar or a character vector
%
% Used by:
%       \Shared\Code\vIRt-Moore\jm_run_virt_sim_replicates.m
%       \Shared\Code\vIRt-Moore\jm_run_virt_sim_vary_ext_current.m
%       \Shared\Code\vIRt-Moore\virt_moore_multiple_reps.m
%

% File History:
% 2026-01-18 Extracted from virt_moore_multiple_reps.m by Gemini.

%% Main Function

% Add leading '-' to get the string to be replaced or appended
addition = ['-', strParamVaried];

% Add leading '-' to create parameter key
paramKey = ['-', cleanParamName];

% Create a regular expression to find an existing entry for this parameter key.
% It matches "-<key>-" followed by one or more characters that are not a hyphen.
pattern = [paramKey, '-[^-]+'];

if regexp(oldName, pattern)
    % If the key already exists, replace the old segment
    newName = regexprep(oldName, pattern, addition);
else
    % If the key doesn't exist, append the new segment
    newName = [oldName, addition];
end

% Remove the seednumber string if it is not the parameter varied
if ~strcmp(cleanParamName, 'seedNumber')
    seedNumPattern = ['-seedNumber', '-[^-]+'];
    newName = regexprep(newName, seedNumPattern, '');    
end

% Remove any nReps parameter-value string
nRepsPattern = ['-nReps', '-[^-]+'];
newName = regexprep(newName, nRepsPattern, '');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%