%% Tg_Strategy_ExtraPlots_SS.m
%  Extra strategy figures + saved comparison stats, all off Build_Tg_StrategyTable:
%   1) Day-N WT vs APP for each strategy and each strategy group, both age
%      groupings - Welch's t-test for comparisons.  Loops over days 1-4.
%   2) Strategy-group usage (WT | APP) for each day, both age groupings.
%   3) Day d1 vs Day d2 within-rat, per strategy group.
%      Consecutive pairs [1 2],[2 3],[3 4] + full span [1 4].
%
%  Stats saved:
%      Stats_DayN_GenoByAge.csv      (Plot 1, all days)
%      Stats_Day3v4_AgeSubplots.csv  (Plot 3 AgeSubplots, all day pairs)
%      Stats_Day3v4_byAge.csv        (Plot 3 WithinRat,   all day pairs)
%
%  USES: Build_Tg_StrategyTable,
%        PLOT_Tg_Day4GenoByAge, PLOT_Tg_Day4Groups_Panel,
%        PLOT_Tg_GroupUsage_Panel,
%        PLOT_Tg_Day3v4_AgeSubplots, PLOT_Tg_Day3v4_WithinRat
%
%  SS 2026

%% INPUT PARAMS
clear all;
base_dir   = 'D:\NARP_Data\RTrack_NARPMale';
strat_xlsx = fullfile(base_dir, 'Tg_MWM_results_06-07-2026');
fig_dir    = fullfile(base_dir, 'Figures', 'ExtraPlots');

ageEdges    = [3.5 4.5 5.5 6.5 7.5 8.5 9.5 10.5 11.5 14.5 17.5 22.5];
ageLabels   = ["4mo","5mo","6mo","7mo","8mo","9mo","10mo","11mo","12-13mo","15-16mo","20-21mo"];
youngMidCut = 9;
midOldCut   = 15;

genoCols   = [0.50 0.70 1.00; 0.00 0.20 0.70];
grpCols    = [0.3961 0.2627 0.1294; 1.0000 0.7020 0.4000; 0 0.5020 0];
greyShades = [0.72 0.72 0.72; 0.45 0.45 0.45; 0.15 0.15 0.15];

% ---- days to iterate -------------------------------------------------------
allDays    = 1:4;                          % Plot 1 + 2: one panel per day
dayPairs   = [1 2; 2 3; 3 4; 1 4];        % Plot 3: consecutive + full span
% ----------------------------------------------------------------------------

dir_p1    = fullfile(fig_dir, 'Day_GenoByAge');
dir_panel = fullfile(fig_dir, 'Panels');
dir_p3    = fullfile(fig_dir, 'Day_v_Day');
dir_stats = fullfile(fig_dir, 'Stats');
for d = {dir_p1, dir_panel, dir_p3, dir_stats}
    if ~exist(d{1},'dir'), mkdir(d{1}); end
end

%% Load + reshape
data1 = readtable(strat_xlsx);
T = Build_Tg_StrategyTable(data1, ...
    'AgeEdges', ageEdges, 'AgeLabels', ageLabels, ...
    'YoungMidCut', youngMidCut, 'MidOldCut', midOldCut);
T = T(~ismissing(T.Group4) & ~ismissing(T.AgeMonth), :);

strategies    = {'thigmotaxis','circling','randomPath','scanning', ...
                 'chaining','directedSearch','correctedPath','directPath'};
grpVars       = {'PlatformIndependent','Procedural','Allocentric'};
grpYLab       = {'P(Platform-Independent)','P(Procedural)','P(Allocentric)'};
coarseLevels  = ["Young","Mid","Old"];

%% Plot 1 - WT vs APP, each day
statsDay = table();
for dy = allDays                                   % LOOP: one iteration per day
    % (a) individual figure per strategy
    for s = 1:numel(strategies)
        for bc = {'AgeGroup','AgeMonth'}
            S = PLOT_Tg_Day4GenoByAge(T, strategies{s}, bc{1}, ...
                'Day', dy, 'Colors', genoCols, 'YLabel', 'P(strategy)', ...
                'SaveDir', fullfile(dir_p1, sprintf('Day%d', dy)));
            statsDay = [statsDay; S]; %#ok<AGROW>
        end
    end
    % (b) 1x3 strategy-group panel (coarse + fine)
    for bc = {'AgeGroup','AgeMonth'}
        S = PLOT_Tg_Day4Groups_Panel(T, bc{1}, 'Day', dy, ...
            'Groups', grpVars, 'Colors', genoCols, ...
            'SaveDir', fullfile(dir_panel, sprintf('Day%d', dy)));
        statsDay = [statsDay; S]; %#ok<AGROW>
    end
end
writetable(statsDay, fullfile(dir_stats, 'Stats_DayN_GenoByAge.csv'));

%% Plot 2 - Strategy-group usage (WT | APP), each day
for dy = allDays                                   % LOOP: one iteration per day
    for bc = {'AgeGroup','AgeMonth'}
        PLOT_Tg_GroupUsage_Panel(T, bc{1}, 'Day', dy, 'Colors', grpCols, ...
            'SaveDir', fullfile(dir_panel, sprintf('Day%d', dy)));
    end
end

%% Plot 3 - Day d1 vs Day d2, per strategy group, all day pairs
statsSub  = table();
statsGrey = table();
for dp = 1:size(dayPairs,1)                        % LOOP: one iteration per day pair
    d1 = dayPairs(dp,1);  d2 = dayPairs(dp,2);
    subDir = fullfile(dir_p3, sprintf('Day%dvDay%d', d1, d2));
    for gv = 1:numel(grpVars)
        S = PLOT_Tg_Day3v4_AgeSubplots(T, grpVars{gv}, ...
            'AgeLevels', coarseLevels, 'Colors', genoCols, 'Days', [d1 d2], ...
            'YLabel', grpYLab{gv}, 'SaveDir', subDir);
        statsSub = [statsSub; S]; %#ok<AGROW>

        S = PLOT_Tg_Day3v4_WithinRat(T, grpVars{gv}, ...
            'GroupBy', 'AgeGroup', 'Levels', coarseLevels, 'Colors', greyShades, ...
            'TitleTag', 'byAge', 'Days', [d1 d2], 'YLabel', grpYLab{gv}, ...
            'SaveDir', subDir);
        statsGrey = [statsGrey; S]; %#ok<AGROW>
    end
end
writetable(statsSub,  fullfile(dir_stats, 'Stats_DayPairs_AgeSubplots.csv'));
writetable(statsGrey, fullfile(dir_stats, 'Stats_DayPairs_byAge.csv'));

fprintf('\nDone. Figures under %s\n       Stats  under %s\n', fig_dir, dir_stats);