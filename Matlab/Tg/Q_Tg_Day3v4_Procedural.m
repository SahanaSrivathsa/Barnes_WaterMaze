function Q_Tg_Day3v4_Procedural()
% Q_Tg_Day3v4_Procedural  Is the Day3->Day4 procedural change opposite in WT vs APP?
%
% Example call:
%   Q_Tg_Day3v4_Procedural()
%
% Q2 (two stats, Old rats, Procedural strategy group):
%   STAT 1 - "how it changes": within-genotype Day3->Day4 change. Produced by
%            PLOT_Tg_Day3v4_AgeSubplots (paired t per genotype x age group), which
%            also draws the Day3-vs-4 figure for Young/Mid/Old.
%   STAT 2 - "how the pattern differs": Genotype x Day interaction from a 2x2
%            mixed ANOVA (between Genotype, within Day in {3,4}). A significant
%            interaction with opposite-sign simple effects is the formal test of
%            "opposite". Run for Old (headline) and, for completeness, all ages.
% The interaction equals a two-sample test on each rat's (Day4-Day3) change
% score; expected Old result F ~ 20.0, p ~ 1e-4 (WT falls, APP rises).
%
% USES: PLOT_Tg_Day3v4_AgeSubplots, STATS_MixedANOVA.
% SS 2026

%% Input params
in_csv  = 'D:\NARP_Data\RTrack_NARPMale\StrategyProcessed\Tg_Strategy_RatDayMeans.csv';
out_dir = fullfile(fileparts(in_csv), 'Stats_Q2_Day3v4');
genoCols = [0.50 0.70 1.00; 0.00 0.20 0.70];
ages     = ["Young","Mid","Old"];
days     = [3 4];
if ~exist(out_dir,'dir'), mkdir(out_dir); end
 
%% Load + set factor orders
T = readtable(in_csv);
T.Genotype = categorical(string(T.Genotype), ["WT","APP"]);
T.AgeGroup = categorical(string(T.AgeGroup), ages);
 
%% STAT 1 - within-genotype Day3->Day4 change (paired) + figure
S1 = PLOT_Tg_Day3v4_AgeSubplots(T, 'Procedural', ...
        'AgeCol','AgeGroup', 'AgeLevels', ages, 'Colors', genoCols, ...
        'Days', days, 'Test', 'ttest', 'YLabel', 'P(Procedural)', 'SaveDir', out_dir);
writetable(S1, fullfile(out_dir,'WithinGenotype_Day3v4_paired.csv'));
 
%% STAT 2 - Day4-Day3 change (Delta): directionality + WT-vs-APP interaction
% The Genotype x Day interaction of a 2-day design equals a two-sample test on
% the per-rat change score Delta = P(Day4) - P(Day3); PLOT_Tg_DeltaByAge runs
% that (verified: Old F = t^2 ~ 19.8, p ~ 1e-4) plus the one-sample Delta-vs-0
% directionality test, both Bonferroni-corrected, and draws the Delta bar+dots.
S2 = PLOT_Tg_DeltaByAge(T, 'Procedural', ...
        'Days', days, 'AgeCol','AgeGroup', 'AgeLevels', ages, ...
        'Colors', genoCols, 'Test', 'ttest2', 'Correction', 'bonferroni', ...
        'SaveDir', out_dir);
writetable(S2, fullfile(out_dir,'Delta_Day3v4_directionality_and_interaction.csv'));
 
%% Print summary
fprintf('\n--- Q_Tg_Day3v4_Procedural complete ---\n');
fprintf('Within-genotype paired : WithinGenotype_Day3v4_paired.csv\n');
fprintf('Delta + interaction    : Delta_Day3v4_directionality_and_interaction.csv\n');
fprintf('Output dir             : %s\n', out_dir);
end
 