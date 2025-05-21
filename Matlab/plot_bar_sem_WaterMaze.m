function plot_bar_sem_WaterMaze(uniqueDays, meanYoung, semYoung, meanOld, semOld, dataYoung, dataOld, tukeyResults)
% plot_bar_sem creates a grouped bar plot with error bars, jittered individual points,
% and significance markers (using sigstar) based on pairwise comparisons.
%
% Inputs:
%   uniqueDays  - Vector of day
%   meanYoung   - Mean values for the young group per day.
%   semYoung    - SEM values for the young group per day.
%   meanOld     - Mean values for the old group per day.
%   semOld      - SEM values for the old group per day.
%   dataYoung   - Cell array of individual data points for young (one cell per day).
%   dataOld     - Cell array of individual data points for old (one cell per day).
%   tukeyResults- Table of Tukey post hoc results (with columns: Age, Day_1, Day_2, pValue, etc.).

% This function plots:
%   - A grouped bar plot with error bars.
%   - Jittered scatter points overlaid on the bars.
%   - Significance markers: 
%       - Between-age comparisons on the same day in gray. - Uses a t test
%       - Within-group day-to-day comparisons -based on tukey post hoc

%SS 2025
%% PARAMETERS THAT MIGHT CHANGE
% colors
 oldColor=[0.2196,0.5569,0.2353]; %green
 youngColor=[0.4157,0.1059,0.6039];% purple
%Extract_varargin; %Change color defaults
%% 
nDays = numel(uniqueDays);
nDays = numel(uniqueDays);
meanData = [meanYoung(:), meanOld(:)];
barWidth = 0.8;
hBar = bar(uniqueDays, meanData, 'grouped', 'BarWidth', barWidth);

hBar(1).FaceColor = youngColor;
hBar(2).FaceColor = oldColor;
hBar(1).FaceAlpha = 0.35;
hBar(2).FaceAlpha = 0.35;

drawnow; % Update bar properties
barCenters = zeros(nDays,2);
for i = 1:2
    barCenters(:,i) = hBar(i).XEndPoints;
end

hold on;

% Plot individual data points with jitter
jitterAmount = 0.1;
for d = 1:nDays
    % Young points
    x_jitter = barCenters(d,1) + (rand(size(dataYoung{d}))-0.5)*jitterAmount;
    scatter(x_jitter, dataYoung{d}, 12, 'MarkerFaceColor', youngColor,...
        'MarkerEdgeColor', 'none', 'MarkerFaceAlpha',0.85);
    % Old points
    x_jitter = barCenters(d,2) + (rand(size(dataOld{d}))-0.5)*jitterAmount;
    scatter(x_jitter, dataOld{d}, 12, 'MarkerFaceColor', oldColor,...
        'MarkerEdgeColor', 'none', 'MarkerFaceAlpha',0.85);
end

% Add error bars
errorbar(barCenters(:,1), meanYoung, semYoung, 'k', 'linestyle', 'none', 'LineWidth', 1.5);
errorbar(barCenters(:,2), meanOld, semOld, 'k', 'linestyle', 'none', 'LineWidth', 1.5);

% Prepare cell arrays for sigstar comparisons
sigPairs = {};
sigP = []; 
sigColors = {};

%% Between-age comparisons-each day
alphaFW = 0.05;            % family-wise alpha
m       = nDays;           % number of day-wise age comparisons

for d = 1:nDays
    [~, p] = ttest2(dataYoung{d}, dataOld{d});
    
    % --- Bonferroni adjust ---
    pAdj = p * m;                       % Bonferroni: p × (#tests)
    pAdj = min(pAdj, 1);                % never exceed 1
    
    if pAdj < alphaFW
        sigPairs{end+1}  = [barCenters(d,1), barCenters(d,2)];
        sigP(end+1)      = pAdj;        % store the *corrected* p
        sigColors{end+1} = [0.55 0.55 0.55];   % gray
    end
end
%% Within-group (day-day) comparisons from Tukey results
% for i = 1:nDays-1
%     j = i+1;  % Only consecutive days
%     % Filter Tukey results for day pair in young and in old separately:
%     idx_y = strcmp(tukeyResults.Age, 'young') & (tukeyResults.Day_1 == uniqueDays(i)) & (tukeyResults.Day_2 == uniqueDays(j));
%     idx_o = strcmp(tukeyResults.Age, 'old')   & (tukeyResults.Day_1 == uniqueDays(i)) & (tukeyResults.Day_2 == uniqueDays(j));
%     if any(idx_y)
%         p_y = min(tukeyResults.pValue(idx_y));
%     else
%         p_y = [];
%     end
%     if any(idx_o)
%         p_o = min(tukeyResults.pValue(idx_o));
%     else
%         p_o = [];
%     end
% 
%     sig_y = ~isempty(p_y) && (p_y < 0.05);
%     sig_o = ~isempty(p_o) && (p_o < 0.05);
% 
%     if sig_y && sig_o
%         % Both groups significant: plot one combined bar spanning from the
%         % average of young and old centers on day i to that on day j.
%         p_comb = min(p_y, p_o);
%         x1 = mean(barCenters(i,:)); % average x for day i (young & old)
%         x2 = mean(barCenters(j,:)); % average x for day j
%         markerColor = [0 0 0];       % black for both significant
%         sigPairs{end+1} = [x1, x2];
%         sigP(end+1) = p_comb;
%         sigColors{end+1} = markerColor;
%     else
%         % Only one group is significant, add individual bar:
%         if sig_y
%             p_comb = p_y;
%             markerColor = youngColor;
%             sigPairs{end+1} = [barCenters(i,1), barCenters(j,1)];
%             sigP(end+1) = p_comb;
%             sigColors{end+1} = markerColor;
%         end
%         if sig_o
%             p_comb = p_o;
%             markerColor = oldColor;
%             sigPairs{end+1} = [barCenters(i,2), barCenters(j,2)];
%             sigP(end+1) = p_comb;
%             sigColors{end+1} = markerColor;
%         end
%     end
% end

% Call sigstar once with the prepared groups and numeric p-values.
hSig = sigstar(sigPairs, sigP);
% Set the color for each significance marker (bar and stars) once.
% for k = 1:length(hSig)
%     set(hSig(k), 'Color', sigColors{k});
%     set(hSig(k,2), 'Color', sigColors{k});
% end

end
