function anovaResults = runMixedANOVA(tbl, measureVars)
% runMixedANOVA performs a mixed-design repeated measures ANOVA.
%
% Within-subject factor table with Day and fit a repeated measures model with Age
% as a between-subject factor.
% Then run ranova with 'WithinModel' set to 'Day' to test the main effect of Day
%and the Age*Day interaction.
%
% INPUTS:
%   tbl         - Table with one row per subject. Must include a column
%                   'Age'- here the group
%                 and repeated measures variables (here 'Day1'-'Day4').
%   measureVars - Cell array of strings naming the repeated measures,
%                   here day1, day2, day 3 etc
%
% OUTPUT:
%   anovaResults - Table of ANOVA results including main effects and interaction.
%

%SS 2024

numLevels = numel(measureVars);
withinDesign = table((1:numLevels)', 'VariableNames', {'Day'});
formula = sprintf('%s-%s ~ Age', measureVars{1}, measureVars{end});
rm = fitrm(tbl, formula, 'WithinDesign', withinDesign);
anovaResults = ranova(rm, 'WithinModel', 'Day');
end
