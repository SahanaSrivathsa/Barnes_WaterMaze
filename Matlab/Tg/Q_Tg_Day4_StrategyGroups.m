function Q_Tg_Day4_StrategyGroups()
% Q_Tg_Day4_StrategyGroups  How WT vs APP strategy use differs by age on Day 4.
%
% Example call:
%   Q_Tg_Day4_StrategyGroups()
%
% Q1. Loads the per-rat-per-day strategy table (Build_Tg_StrategyTable output),
% then:
%   (1) HEADLINE 4-way mixed ANOVA  AgeGroup x StrategyGroup x Genotype x Day
%       (between: Genotype, AgeGroup; within: Day, StrategyGroup). NOTE the
%       three strategy-group probabilities sum to ~1 (compositional), so the
%       StrategyGroup factor is structurally non-independent -- read the 4-way
%       as descriptive context only.
%   (2) VALID FOLLOW-UP: three AgeGroup x Genotype x Day mixed ANOVAs, one per
%       strategy group (avoids the sum-to-1 problem).
%   (3) Planned pairwise WT vs APP on Day 4 within each AgeGroup x StrategyGroup
%       (9 cells), Bonferroni m=9, drawn via PLOT_Tg_Day4Groups_Panel with only
%       corrected p-values feeding sigstar.
%
% USES: STATS_MixedANOVA, PLOT_Tg_Day4Groups_Panel.
% SS 2026

%% Input params
in_csv  = 'D:\NARP_Data\RTrack_NARPMale\StrategyProcessed\Tg_Strategy_RatDayMeans.csv';   
out_dir = fullfile(fileparts(in_csv), 'Stats_Q1_Day4');
genoCols = [0.50 0.70 1.00; 0.00 0.20 0.70];                      % WT light / APP dark
grpVars  = {'PlatformIndependent','Procedural','Allocentric'};
grpLab   = {'Platform-Independent','Procedural','Allocentric'};
if ~exist(out_dir,'dir'), mkdir(out_dir); end
 
%% Load + set factor orders
T = readtable(in_csv);
T.Genotype = categorical(string(T.Genotype), ["WT","APP"]);
T.AgeGroup = categorical(string(T.AgeGroup), ["Young","Mid","Old"]);
fineOrder  = ["4mo","5mo","6mo","7mo","8mo","9mo","10mo","11mo","12-13mo","15-16mo","20-21mo"];
T.AgeMonth = categorical(string(T.AgeMonth), fineOrder);
 
%% (1) Day-4 VALUES table (Genotype x AgeGroup x StrategyGroup: mean, SEM, n)
L  = stack(T, grpVars, 'NewDataVariableName','Value', 'IndexVariableName','StrategyGroup');
d4 = L(L.Day==4, :);
V  = groupsummary(d4, {'StrategyGroup','AgeGroup','Genotype'}, {'mean','std'}, 'Value');
V.SEM = V.std_Value ./ sqrt(V.GroupCount);
V = renamevars(V, {'mean_Value','GroupCount'}, {'mean','n'});
V = removevars(V, 'std_Value');
writetable(V, fullfile(out_dir,'Day4_Values_GenoAgeStrat.csv'));
 
%% (2) Per-strategy-group Genotype x AgeGroup two-way ANOVA (Type III, between)
A = table();
for gi = 1:numel(grpVars)                         % LOOP: 3 strategy groups
    Tg = T(T.Day==4, :);
    Rg = STATS_TwoWayANOVA(Tg.(grpVars{gi}), Tg.Genotype, Tg.AgeGroup, ...
            'Names', {'Genotype','AgeGroup'});
    Rg.effects.StrategyGroup = repmat(string(grpLab{gi}), height(Rg.effects), 1);
    A = [A; Rg.effects]; %#ok<AGROW>
end
A = movevars(A, 'StrategyGroup', 'Before', 'Effect');
writetable(A, fullfile(out_dir,'ANOVA_Day4_perGroup_GenoxAge.csv'));
 
%% (3) Planned Day-4 pairwise (Bonferroni) + figures: coarse AND fine age
S_coarse = PLOT_Tg_Day4Groups_Panel(T, 'AgeGroup', 'Day', 4, ...
        'Groups', grpVars, 'GroupLabels', grpLab, 'Colors', genoCols, ...
        'Test', 'ttest2', 'Correction', 'bonferroni', 'SaveDir', out_dir);
writetable(S_coarse, fullfile(out_dir,'Pairwise_Day4_WTvsAPP_coarse_Bonferroni.csv'));
 
% Fine bins are sparse (some cells n<2 are skipped); correction is across all
% valid fine tests (~30), so almost everything drops out -- this is expected.
S_fine = PLOT_Tg_Day4Groups_Panel(T, 'AgeMonth', 'Day', 4, ...
        'Groups', grpVars, 'GroupLabels', grpLab, 'Colors', genoCols, ...
        'Test', 'ttest2', 'Correction', 'bonferroni', 'SaveDir', out_dir);
writetable(S_fine, fullfile(out_dir,'Pairwise_Day4_WTvsAPP_fine_Bonferroni.csv'));
 
%% Print summary
fprintf('\n--- Q_Tg_Day4_StrategyGroups complete ---\n');
fprintf('Day-4 values    : Day4_Values_GenoAgeStrat.csv\n');
fprintf('per-group ANOVAs: ANOVA_Day4_perGroup_GenoxAge.csv\n');
fprintf('Day-4 pairwise  : Pairwise_Day4_WTvsAPP_coarse_Bonferroni.csv + ..._fine_Bonferroni.csv\n');
fprintf('Output dir      : %s\n', out_dir);
end