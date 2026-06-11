function [hBar, ctr] = HELP_DrawGroupedBars(ax, vals, colors, seriesNames, xlabels, faceAlpha, jit)
% [hBar, ctr] = HELP_DrawGroupedBars(ax, vals, colors, seriesNames, xlabels, faceAlpha, jit)
%
% Draw grouped bars (across-rat mean), SEM error bars, and jittered per-rat
% dots into axes AX. Shared by the single-figure and panel bar plotters.
%
% INPUT
%   ax          : target axes
%   vals        : nX x nSeries cell, each a vector of per-rat values
%   colors      : nSeries x 3
%   seriesNames : string array (legend / DisplayName per series)
%   xlabels     : x tick labels (length nX)
%   faceAlpha   : bar face alpha (e.g. 0.75)
%   jit         : scatter x-jitter width (e.g. 0.12)
% OUTPUT
%   hBar : bar handles (1 x nSeries)
%   ctr  : nX x nSeries bar-centre x positions (for sigstar / annotation)
%
% SS 2026

[nX, nS] = size(vals);
M = nan(nX,nS); SEM = nan(nX,nS);
for x = 1:nX                                     % LOOP: cells x series (small)
    for s = 1:nS
        v = vals{x,s}; v = v(~isnan(v));
        if ~isempty(v), M(x,s) = mean(v); SEM(x,s) = std(v)/sqrt(numel(v)); end
    end
end

hBar = bar(ax, 1:nX, M, 'grouped', 'BarWidth', 0.9);
for s = 1:nS
    hBar(s).FaceColor   = colors(s,:);
    hBar(s).FaceAlpha   = faceAlpha;
    hBar(s).DisplayName = char(seriesNames(s));
end
drawnow;
ctr = nan(nX,nS);
for s = 1:nS, ctr(:,s) = hBar(s).XEndPoints'; end

hold(ax, 'on');
for s = 1:nS
    errorbar(ax, ctr(:,s), M(:,s), SEM(:,s), 'k', 'LineStyle','none', ...
        'LineWidth', 1.2, 'HandleVisibility','off');
end
for x = 1:nX
    for s = 1:nS
        v = vals{x,s}; if isempty(v), continue; end
        xj = ctr(x,s) + (rand(size(v))-0.5)*jit;
        scatter(ax, xj, v, 16, colors(s,:), 'filled', ...
            'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.6, 'HandleVisibility','off');
    end
end
xticks(ax, 1:nX); xticklabels(ax, xlabels);
end