% From workspace
S=slopeResults.window_10_60s;
% Colors
cYoung = [0.2196, 0.5569, 0.2353];
cOld   = [0.4157, 0.1059, 0.6039];


young = double(S.youngSlopes(:));
old   = double(S.oldSlopes(:));

% drop NaNs just to be safe
young = young(~isnan(young));
old   = old(~isnan(old));

% stats
m  = [mean(young,'omitnan'), mean(old,'omitnan')];
se = [std(young,0,'omitnan')/sqrt(numel(young)), ...
      std(old,  0,'omitnan')/sqrt(numel(old))];

% Welch two-sample t-test 
[~, p_value] = ttest2(young, old, 'Vartype','unequal');  

%Plot: bar + SEM + jittered scatter
jitSigma = 0.04;   % jitter spread
capSize  = 10;     % errorbar cap
barAlpha = 0.60;

figure('Color','w','Units','pixels','Position',[100 100 640 480]);
ax = axes; hold(ax,'on');

b = bar(1:2, m, 'FaceColor','flat');            % XEndPoints used to center error bars
b.CData     = [cYoung; cOld];
b.FaceAlpha = barAlpha;


% jittered circular markers
x1 = 1 + jitSigma*randn(size(young));
x2 = 2 + jitSigma*randn(size(old));
scatter(x1, young, 36, cYoung, 'o', 'filled', 'MarkerEdgeColor','k');   % :contentReference[oaicite:3]{index=3}
scatter(x2, old,   36, cOld,   'o', 'filled', 'MarkerEdgeColor','k');
errorbar(b.XEndPoints, m, se, 'k', 'LineStyle','none', ...
    'CapSize',capSize, 'LineWidth',3);                                  

set(ax,'XLim',[0.5 2.5], 'XTick',1:2, 'XTickLabel',{'Young','Old'});
ylabel(ax,'Slope');
title(ax,'Slope of Speed (10-60s)');


% p-value label
yTop = max(m + se);
dy   = max(0.05*yTop, 0.2);
ylim(ax, [min(0, min([young;old])*1.05), yTop + 1.2*dy]);
text(mean(b.XEndPoints), yTop + 0.2*dy, sprintf('p = %.4f', p_value), ...
    'HorizontalAlignment','center','Parent',ax);
pubify_figure_axis_robust(18,18)
