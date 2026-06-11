%% Tg_Strategy_ExtraPlots_SS.m
%  Three extra strategy figures, all off Build_Tg_StrategyTable:
%   1) Day-4 WT vs APP per strategy, across coarse AND fine age bins
%   2) WT-only strategy-group usage (collapsed over days) across age bins
%   3) Day3 vs Day4 within-rat lines for Procedural, one fig per coarse age group
%
%  USES:
%    Build_Tg_StrategyTable, PLOT_Tg_Day4GenoByAge,
%    PLOT_Tg_GroupUsageByAge, PLOT_Tg_Day3v4_WithinRat
%
%  SS 2026


%% INPUT PARAMS
clear all;
base_dir   = 'D:\NARP_Data\RTrack_NARPMale';
strat_xlsx = fullfile(base_dir, 'Tg_MWM_results_06-07-2026');   % add .xlsx if readtable errors
fig_dir    = fullfile(base_dir, 'Figures', 'ExtraPlots');

ageEdges   = [3.5 4.5 5.5 6.5 7.5 8.5 9.5 10.5 11.5 14.5 17.5 22.5];
ageLabels  = ["4mo","5mo","6mo","7mo","8mo","9mo","10mo","11mo","12-13mo","15-16mo","20-21mo"];
youngMidCut = 9;
midOldCut   = 15;

% Colours
genoCols = [0.50 0.70 1.00;    % WT  light blue
            0.00 0.20 0.70];   % APP dark blue
grpCols  = [0.3961 0.2627 0.1294;   % brown  - platform-independent
            1.0000 0.7020 0.4000;   % orange - procedural
            0      0.5020 0     ];   % green  - allocentric
greyShades = [0.72 0.72 0.72;       % Young light grey
              0.45 0.45 0.45;       % Mid   mid grey
              0.15 0.15 0.15];      % Old   dark grey

% Output subfolders
dir_p1    = fullfile(fig_dir, 'Day4_GenoByAge');
dir_p2    = fullfile(fig_dir, 'GroupUsageByAge');
dir_p3    = fullfile(fig_dir, 'Day3v4');
dir_stats = fullfile(fig_dir, 'Stats');
for d = {dir_p1, dir_p2, dir_p3, dir_stats}
    if ~exist(d{1},'dir'), mkdir(d{1}); end
end

%% Load + reshape
data1 = readtable(strat_xlsx);
T = Build_Tg_StrategyTable(data1, ...
    'AgeEdges', ageEdges, 'AgeLabels', ageLabels, ...
    'YoungMidCut', youngMidCut, 'MidOldCut', midOldCut);
T = T(~ismissing(T.Group4) & ~ismissing(T.AgeMonth), :);

strategies = {'thigmotaxis','circling','randomPath','scanning', ...
              'chaining','directedSearch','correctedPath','directPath'};
grpVars    = {'PlatformIndependent','Procedural','Allocentric'};
grpYLab    = {'P(Platform-Independent)','P(Procedural)','P(Allocentric)'};

%% Plot 1 - Day-4 WT vs APP for each strategy AND each strategy group
plot1Vars = [strategies, grpVars];
plot1YLab = [repmat({'P(strategy)'},1,numel(strategies)), repmat({'P(strategy group)'},1,numel(grpVars))];
statsDay4 = table();
for v = 1:numel(plot1Vars)                       % LOOP: per variable
    for bc = {'AgeGroup','AgeMonth'}
        S = PLOT_Tg_Day4GenoByAge(T, plot1Vars{v}, bc{1}, ...
            'Day', 4, 'Colors', genoCols, 'YLabel', plot1YLab{v}, 'SaveDir', dir_p1);
        if isempty(statsDay4), statsDay4 = S; else, statsDay4 = [statsDay4; S]; end %#ok<AGROW>
    end
end
writetable(statsDay4, fullfile(dir_stats, 'Stats_Day4_GenoByAge.csv'));

%% Plot 2 - strategy-group usage, Day 4, WT and APP, coarse + fine (4 figures)
for geno = ["WT","APP"]
    for bc = {'AgeGroup','AgeMonth'}
        PLOT_Tg_GroupUsageByAge(T, bc{1}, 'Genotype', geno, 'Day', 4, ...
            'Colors', grpCols, 'SaveDir', dir_p2);
    end
end
for bc = {'AgeGroup','AgeMonth'}
    PLOT_Tg_GroupUsage_Panel(T, bc{1}, 'Day', 4, 'Colors', grpCols, 'SaveDir', dir_panel);
end
 
%% Plot 3 - Day3 vs Day4 within-rat, per strategy group
coarseLevels = ["Young","Mid","Old"];
statsSub  = table();
statsGrey = table();
for gv = 1:numel(grpVars)                         % LOOP: per strategy group
    % (a) Young/Mid/Old as 3 subplots, WT vs APP
    S = PLOT_Tg_Day3v4_AgeSubplots(T, grpVars{gv}, ...
        'AgeLevels', coarseLevels, 'Colors', genoCols, 'Days', [3 4], ...
        'YLabel', grpYLab{gv}, 'SaveDir', dir_p3);
    if isempty(statsSub), statsSub = S; else, statsSub = [statsSub; S]; end %#ok<AGROW>

    % (b) age groups overlaid in grey shades (genotype pooled)
    S = PLOT_Tg_Day3v4_WithinRat(T, grpVars{gv}, ...
        'GroupBy', 'AgeGroup', 'Levels', coarseLevels, 'Colors', greyShades, ...
        'TitleTag', 'byAge', 'Days', [3 4], 'YLabel', grpYLab{gv}, 'SaveDir', dir_p3);
    if isempty(statsGrey), statsGrey = S; else, statsGrey = [statsGrey; S]; end %#ok<AGROW>
end
writetable(statsSub,  fullfile(dir_stats, 'Stats_Day3v4_AgeSubplots.csv'));
writetable(statsGrey, fullfile(dir_stats, 'Stats_Day3v4_byAge.csv'));

fprintf('\nDone. Figures under %s ; stats under %s\n', fig_dir, dir_stats);