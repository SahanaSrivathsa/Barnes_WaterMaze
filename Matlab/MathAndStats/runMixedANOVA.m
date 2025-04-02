function anovaResults = runMixedANOVA(tbl, measureVars, factorName)
% runMixedANOVA performs a mixed-design repeated measures ANOVA.
%
% Within-subject factor table with Day and fit a repeated measures model with Age
% as a between-subject factor.
% Then run ranova with 'WithinModel' set to Day to test the main effect of Day
% and the Age*Day interaction.
%
% INPUTS:
%   tbl         - Table with one row per subject. Must include a column
%                 'Age' (the group) and repeated measures variables (e.g., 'Day1'-'Day4').
%   measureVars - Cell array of strings naming the repeated measures (e.g., 'Day1', 'Day2', 'Day3', etc.)
%   factorName  - (Optional) String for the within-subject factor (default is 'Day').
%
% OUTPUT:
%   anovaResults - Table of ANOVA results including main effects and interaction.
%
%SS 2024

if nargin < 3
    factorName = 'Day';
end

numLevels = numel(measureVars);
withinDesign = table((1:numLevels)', 'VariableNames', {factorName});
formula = sprintf('%s-%s ~ Age', measureVars{1}, measureVars{end});
rm = fitrm(tbl, formula, 'WithinDesign', withinDesign);
anovaResults = ranova(rm, 'WithinModel', factorName);
end
