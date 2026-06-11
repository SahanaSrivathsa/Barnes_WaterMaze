function [hMean, pairByLevel] = HELP_DrawDay3v4Axes(ax, rat, ratCol, grpCol, levels, days, colors, indivAlpha)
% [hMean, pairByLevel] = HELP_DrawDay3v4Axes(ax, rat, ratCol, grpCol, levels, days, colors, indivAlpha)
%
% Draw within-rat day d1->d2 trajectories into axes AX: one faint line+circles
% per rat in its level colour, plus a bold mean line per level. Shared by the
% single-figure and age-subplot day3-vs-day4 plotters.
%


d1 = days(1); d2 = days(2);
nL = numel(levels);
hMean = gobjects(nL,1);
pairByLevel = cell(nL,1);
gstr = string(rat.(grpCol));
hold(ax, 'on');
for li = 1:nL                                    % LOOP: per level
    col  = colors(li,:);
    rl   = rat(gstr == levels(li), :);
    rats = unique(string(rl.Animal));
    P    = nan(numel(rats), 2);
    for r = 1:numel(rats)                        % LOOP: per rat (faint line)
        v1 = rl.(ratCol)(string(rl.Animal)==rats(r) & rl.Day==d1);
        v2 = rl.(ratCol)(string(rl.Animal)==rats(r) & rl.Day==d2);
        if isempty(v1) || isempty(v2), continue; end
        P(r,:) = [v1(1) v2(1)];
        plot(ax, [d1 d2], P(r,:), '-o', 'Color', [col indivAlpha], ...
            'MarkerFaceColor', col, 'MarkerEdgeColor','none', ...
            'MarkerSize', 4, 'LineWidth', 1, 'HandleVisibility','off');
    end
    pairByLevel{li} = P;
    m = mean(P, 1, 'omitnan');
    hMean(li) = plot(ax, [d1 d2], m, '-o', 'Color', col, 'LineWidth', 3, ...
        'MarkerFaceColor', col, 'MarkerSize', 8, 'DisplayName', char(levels(li)));
end
xlim(ax, [d1-0.3 d2+0.3]); xticks(ax, [d1 d2]);
end
 