function postHocResults = runTukeyPostHocMixed(tbl, measureVars, factorVar)
% runTukeyPostHocMixed runs Tukey's post-hoc pairwise comparisons on the 
% within-subject factor (default 'Day') within each Age group for a 
% mixed-design repeated measures ANOVA.
%
% Create the within-subject factor table and fit the repeated measures model
% with Age as a between-subject factor. Use multcompare to run Tukey comparisons 
% among levels of the within-subject factor.
%
% INPUTS:
%   tbl         - Table with one row per subject. Must include a column 'Age'
%                 (the group) and repeated measures variables (e.g., day1, day2, day3, etc.)
%   measureVars - Cell array of strings naming the repeated measures.
%                 e.g., {'Day1', 'Day2', 'Day3', ...}
%   factorVar   - (Optional) String naming the within-subject factor. Default is 'Day'.
%
% OUTPUT:
%   postHocResults - The results of the Tukey post hoc tests.
%
% SS 2025

if nargin < 3
    factorVar = 'Day';
end

numLevels = numel(measureVars);
withinDesign = table((1:numLevels)', 'VariableNames', {factorVar});
formula = sprintf('%s-%s ~ Age', measureVars{1}, measureVars{end});
rm = fitrm(tbl, formula, 'WithinDesign', withinDesign);
postHocResults = multcompare(rm, factorVar, 'By', 'Age', 'ComparisonType', 'tukey-kramer');
end
