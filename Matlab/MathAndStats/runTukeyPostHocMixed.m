function postHocResults = runTukeyPostHocMixed(tbl, measureVars)
% runTukeyPostHocMixed runs Tukeys post-hoc pairwise comparisons on the 
% Day factor within each Age group for a mixed-design repeated measures ANOVA.
% Originally designed for WaterMaze and W maze data
%
%
% Create the within-subject factor table and fit the repmeasures model with Age as a between-subject factor.
%Use multcompareto run Tukey comparisons among Day levels for each Age group.
%
% INPUTS:
%   tbl         - Table with one row per subject.Group here is 'Age' 
%               and repeated measures variables- here is day
%   measureVars - Cell array of strings naming the repeated measures,
%                   here day1, day2, day 3 etc
% OUTPUT:
%   postHocResults - The results of the Tukey post hoc tests.
% SS 2025
 
numLevels = numel(measureVars);
withinDesign = table((1:numLevels)', 'VariableNames', {'Day'});
formula = sprintf('%s-%s ~ Age', measureVars{1}, measureVars{end});
rm = fitrm(tbl, formula, 'WithinDesign', withinDesign);
postHocResults = multcompare(rm, 'Day', 'By', 'Age', 'ComparisonType', 'tukey-kramer');
end
