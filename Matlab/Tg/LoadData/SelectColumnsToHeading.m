function [Tout, missingCols] = SelectColumnsToHeading(Tin, heading)
% helper
%Reads tables written for Tg rats
heading = string(heading(:))';                 % row of strings, preserve order
have   = string(Tin.Properties.VariableNames);
 
missingCols = setdiff(heading, have, 'stable');   % heading cols not in this file
missingCols = missingCols(:);    
 
for c = 1:numel(missingCols)
    Tin.(char(missingCols(c))) = nan(height(Tin), 1);
end
 
% Hard error if, after filling, a heading column still isn't present
stillMissing = setdiff(heading, string(Tin.Properties.VariableNames), 'stable');
if ~isempty(stillMissing)
    error('HELP_SelectColumnsToSchema:headingUnmet', ...
        'Could not satisfy heading columns: %s', strjoin(cellstr(stillMissing), ', '));
end
 
Tout = Tin(:, cellstr(heading));         