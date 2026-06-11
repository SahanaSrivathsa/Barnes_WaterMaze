function PLOT_MWM_CoefAgeSummary(ageResults, ageLabels, ageCenters, varName, varargin)
%
% Coefficient-vs-age summary plot from a set of per-age LME fits. 
% For each non-intercept fixed-effect term, plots Estimate +/- 95% CI 
% vs the age-bin centre

% ---- defaults ----
p.SaveDir      = '';
p.ExcludeTerms = "(Intercept)";
for k = 1:2:numel(varargin), p.(varargin{k}) = varargin{k+1}; end

ageLabels  = string(ageLabels(:));
ageCenters = ageCenters(:);

% 
firstR = [];
for a = 1:numel(ageResults)                         
    if ~isempty(ageResults{a}), firstR = ageResults{a}; break; end
end
if isempty(firstR)
    warning('PLOT_MWM_CoefAgeSummary:noResults', 'No fitted models for %s.', varName);
    return
end

allTerms = string(firstR.fixedEffects.Term);
terms    = allTerms(~ismember(allTerms, p.ExcludeTerms));
nTerm    = numel(terms);

%
nAge = numel(ageResults);
est  = nan(nTerm, nAge);
lo   = nan(nTerm, nAge);
hi   = nan(nTerm, nAge);
for a = 1:nAge  % LOOP: per-age extraction
    R = ageResults{a};
    if isempty(R), continue; end
    feT = R.fixedEffects;
    feNames = string(feT.Term);
    for t = 1:nTerm   
        idx = feNames == terms(t);
        if any(idx)
            est(t,a) = feT.Estimate(idx);
            lo(t,a)  = feT.Lower(idx);
            hi(t,a)  = feT.Upper(idx);
        end
    end
end


nCols = ceil(sqrt(nTerm));
nRows = ceil(nTerm / nCols);

fig = figure('Name', ['Coef vs Age: ' varName], 'Position', [100 80 1200 800]);
tl  = tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');

for t = 1:nTerm                                      % LOOP: one panel per term (plotting)
    nexttile; hold on
    yline(0, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);   % effect = 0 reference

    valid = ~isnan(est(t,:)) & ~isnan(ageCenters(:)');
    if any(valid)
        x   = ageCenters(valid);
        y   = est(t,valid);
        neg = y - lo(t,valid);
        pos = hi(t,valid) - y;
        errorbar(x, y, neg, pos, 'o-', 'LineWidth', 1.8, ...
            'MarkerFaceColor', 'k', 'Color', 'k', 'CapSize', 6);
    end
    title(terms(t), 'Interpreter', 'none', 'FontSize', 11);
    xlim([min(ageCenters,[],'omitnan')-1, max(ageCenters,[],'omitnan')+1]);
    if exist('pubify_figure_axis_robust','file'), pubify_figure_axis_robust(10,10); end
    hold off
end

xlabel(tl, 'Age (months)',                'FontSize', 14, 'FontWeight', 'bold');
ylabel(tl, 'Coefficient estimate (95% CI)','FontSize', 14, 'FontWeight', 'bold');
title(tl,  ['Per-age coefficients: ' varName], 'FontSize', 16, ...
           'FontWeight', 'bold', 'Interpreter', 'none');

if ~isempty(p.SaveDir)
    if ~exist(p.SaveDir, 'dir'), mkdir(p.SaveDir); end
    saveas(fig, fullfile(p.SaveDir, ['MWM_CoefAge_' varName]), 'png');
end
close(fig);
end