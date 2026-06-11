function PLOT_MWM_AgeSubplots(T, varName, varargin)
% PLOT_MWM_AgeSubplots(T, 'Platform_CIPL', 'SaveDir', proc_dir)
%
% One subplot (tile) per age-group in T.AgeGroup. Within each tile the four
% Sex x Genotype groups (T.Group4) are drawn as learning curves across Day,
% mean +/- SEM with rat as the unit of replication.
%
% COLOUR SCHEME:
%   F WT  -- light red    (sex: female -> red,  genotype: WT  -> light)
%   F APP -- dark red     (sex: female -> red,  genotype: APP -> dark)
%   M WT  -- light blue   (sex: male   -> blue, genotype: WT  -> light)
%   M APP -- dark blue    (sex: male   -> blue, genotype: APP -> dark)
%
% Override with the 'Colors' name-value arg.
% 2026

% ---- defaults ----
p.SaveDir = '';
p.Groups  = ["F WT","F APP","M WT","M APP"];
p.Colors  = [1.00 0.60 0.60;    % F WT  light red / pink
             0.75 0.00 0.00;    % F APP dark red
             0.50 0.70 1.00;    % M WT  light blue
             0.00 0.20 0.70];   % M APP dark blue / navy
p.YLabel  = varName;
p.NCols   = [];
for k = 1:2:numel(varargin)              % LOOP: tiny fixed arg list
    p.(varargin{k}) = varargin{k+1};
end

% ---- validate required columns ----
req  = ["Animal","Day","AgeGroup","Group4", string(varName)];
have = string(T.Properties.VariableNames);
miss = req(~ismember(req, have));
if ~isempty(miss)
    error('PLOT_MWM_AgeSubplots:missingCols', ...
        'T missing required column(s): %s', strjoin(cellstr(miss), ', '));
end
if ~iscategorical(T.AgeGroup)
    error('PLOT_MWM_AgeSubplots:ageType','T.AgeGroup must be categorical (use AssignAgeBin).');
end
T.Group4 = categorical(string(T.Group4), p.Groups);   % enforce order; drop unknowns

% ---- per-rat per-day mean, then across-rat mean / SEM ----
ratDay = groupsummary(T, {'AgeGroup','Group4','Animal','Day'}, 'mean', varName);
ratCol = ['mean_' varName];
grp    = groupsummary(ratDay, {'AgeGroup','Group4','Day'}, {'mean','std','nnz'}, ratCol);
grp.SEM = grp.(['std_' ratCol]) ./ sqrt(grp.(['nnz_' ratCol]));

% ---- age groups that actually have data, in ordinal order ----
ageCats = categories(removecats(T.AgeGroup));
ageCats = ageCats(ismember(ageCats, unique(cellstr(ratDay.AgeGroup))));
nAge    = numel(ageCats);
if nAge == 0
    error('PLOT_MWM_AgeSubplots:noData','No populated age groups for %s.', varName);
end
if isempty(p.NCols), nCols = ceil(sqrt(nAge)); else, nCols = p.NCols; end
nRows = ceil(nAge / nCols);

yTop = max(grp.(['mean_' ratCol]) + grp.SEM, [], 'omitnan');
if isempty(yTop) || ~isfinite(yTop), yTop = 1; end

fig = figure('Name', ['MWM age subplots: ' varName], 'Position', [80 60 1200 800]);
tl  = tiledlayout(nRows, nCols, 'TileSpacing','compact', 'Padding','compact');

hLines = gobjects(numel(p.Groups), 1);
for a = 1:nAge                                  % LOOP: one tile per age group (plotting)
    nexttile; hold on
    isAge = grp.AgeGroup == ageCats{a};
    for g = 1:numel(p.Groups)                   % LOOP: 4 group lines (plotting)
        m = grp(isAge & grp.Group4 == p.Groups(g), :);
        if isempty(m), continue; end
        m = sortrows(m, 'Day');
        h = errorbar(m.Day, m.(['mean_' ratCol]), m.SEM, '-o', ...
            'Color', p.Colors(g,:), 'LineWidth', 2.5, ...
            'MarkerFaceColor', p.Colors(g,:), 'DisplayName', p.Groups(g));
        if ~isgraphics(hLines(g)) || ~isvalid(hLines(g)), hLines(g) = h; end
    end
    title(ageCats{a}, 'FontWeight','bold');
    xlim([0.9 max(T.Day)+0.1]); xticks(1:max(T.Day));
    if isfinite(yTop) && yTop > 0, ylim([0 yTop]); end
    if exist('pubify_figure_axis_robust','file'), pubify_figure_axis_robust(12,12); end
    hold off
end

xlabel(tl, 'Day', 'FontSize', 14, 'FontWeight','bold');
ylabel(tl, p.YLabel, 'FontSize', 14, 'FontWeight','bold', 'Interpreter','tex');
title(tl, varName, 'FontSize', 16, 'FontWeight','bold', 'Interpreter','none');
valid = arrayfun(@(x) isgraphics(x) && isvalid(x), hLines);
legend(hLines(valid), cellstr(p.Groups(valid)), 'Orientation','horizontal', ...
       'Location','southoutside', 'FontSize', 12);

if ~isempty(p.SaveDir)
    if ~exist(p.SaveDir,'dir'), mkdir(p.SaveDir); end
    saveas(fig, fullfile(p.SaveDir, ['MWM_AgeSubplots_' varName]), 'png');
end
close(fig)
end