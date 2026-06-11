%% Tg_Strategy_AgeGroupAnalysis_SS.m
%  Runs PLOT_MWM_AgeSubplots + PLOT_MWM_CoefAgeSummary for STRATEGY data:
% - 8  strategies
% - 3 strategy groups (Platform-Independent/Procedural/Allocentric)
% - entropy (all strategies entropy, 3-group )
%
%  Age-subplots 
%    Coarse-Young / Mid / Old
%    FineMonths-per  month bin (uses AgeAssignBin )
%  USES:
%    Build_Tg_StrategyTable, AssignAgeBin, STATS_MWM_LME,
%    PLOT_MWM_AgeSubplots, PLOT_MWM_CoefAgeSummary
%
%  SS 2026

%% INPUT PARAMS
clear all;
base_dir   = 'D:\NARP_Data\RTrack_NARPMale';
strat_xlsx = fullfile(base_dir, 'Tg_MWM_results_06-07-2026');   % add .xlsx if readtable errors
proc_dir   = fullfile(base_dir, 'StrategyProcessed');
fig_dir    = fullfile(base_dir, 'Figures', 'AgeSummary');
 
% Month bins (same as Analyze_Tg)
ageEdges   = [3.5 4.5 5.5 6.5 7.5 8.5 9.5 10.5 11.5 14.5 17.5 22.5];
ageLabels  = ["4mo","5mo","6mo","7mo","8mo","9mo","10mo","11mo","12-13mo","15-16mo","20-21mo"];
ageCenters = (ageEdges(1:end-1) + ageEdges(2:end)) / 2;
 
% Month cutoffs
youngMidCut = 9;
midOldCut   = 15;
 
% Output subfolders (all under base_dir)
dir_coarse = fullfile(fig_dir, 'Coarse');
dir_fine   = fullfile(fig_dir, 'FineMonths');
dir_coef   = fullfile(fig_dir, 'CoefAge');
for d = {proc_dir, dir_coarse, dir_fine, dir_coef}
    if ~exist(d{1},'dir'), mkdir(d{1}); end
end
 
%% Load and reshape
data1 = readtable(strat_xlsx);
T = Build_Tg_StrategyTable(data1, ...
    'AgeEdges', ageEdges, 'AgeLabels', ageLabels, ...
    'YoungMidCut', youngMidCut, 'MidOldCut', midOldCut);
 
% Drop rows with no metadata
T = T(~ismissing(T.Group4) & ~ismissing(T.AgeMonth), :);
 
% Same rows, different tiling factor for the two subplot styles
Tcoarse = T;                      % AgeGroup = Young/Mid/Old
Tfine   = T;                      % AgeGroup = month bin
Tfine.AgeGroup = T.AgeMonth;
 
%% Groups
dvDefs = { ...
    'thigmotaxis',          'P(strategy)'; ...
    'circling',             'P(strategy)'; ...
    'randomPath',           'P(strategy)'; ...
    'scanning',             'P(strategy)'; ...
    'chaining',             'P(strategy)'; ...
    'directedSearch',       'P(strategy)'; ...
    'correctedPath',        'P(strategy)'; ...
    'directPath',           'P(strategy)'; ...
    'PlatformIndependent',  'P(strategy group)'; ...
    'Procedural',           'P(strategy group)'; ...
    'Allocentric',          'P(strategy group)'; ...
    'Entropy_All',          'Entropy (bits)'; ...
    'Entropy_Group',        'Entropy (bits)' };
nDV = size(dvDefs,1);
 
%% Per-rat per-day means
dvNames = dvDefs(:,1)';
RatDay  = groupsummary(T, {'Animal','Sex','Genotype','AgeMonth','AgeGroup','Group4','Day'}, ...
    'mean', dvNames);
for v = 1:nDV    % rename mean_X -> X
    old = ['mean_' dvNames{v}];
    [tf, idx] = ismember(old, RatDay.Properties.VariableNames);
    if tf, RatDay.Properties.VariableNames{idx} = dvNames{v}; end
end
writetable(RatDay, fullfile(proc_dir, 'Tg_Strategy_RatDayMeans.csv'));
 
% Month bins with data
ages_present = categories(removecats(RatDay.AgeMonth));
[~, ix] = ismember(string(ages_present), ageLabels);
centers_pres = ageCenters(ix);
 
%% LOOP THROUGH DVs
for v = 1:nDV
    vn   = dvDefs{v,1};
    ylab = dvDefs{v,2};
 
    % Age-subplots: coarse (Young/Mid/Old) and fine (month bins)
    PLOT_MWM_AgeSubplots(Tcoarse, vn, 'SaveDir', dir_coarse, 'YLabel', ylab);
    PLOT_MWM_AgeSubplots(Tfine,   vn, 'SaveDir', dir_fine,   'YLabel', ylab);
 
    % Per-month LME - coefficient-vs-age summary
    ageResults = cell(numel(ages_present), 1);
    for a = 1:numel(ages_present)    % LOOP: one fitlme per month bin
        ageLab = ages_present{a};
        Tage   = RatDay(RatDay.AgeMonth == ageLab, :);
        if height(Tage) < 4, continue; end
        try
            R = STATS_MWM_LME(Tage, vn, 'Factors', "Genotype");
            ageResults{a} = R;
            writetable(R.fixedEffects, fullfile(proc_dir, ...
                sprintf('LME_Strategy_%s_%s.csv', vn, strrep(ageLab,' ','_'))));
            fprintf('\n[%s] %-20s @ %-8s | n_obs=%3d n_rats=%2d | resid=%s', ...
                R.model_type, vn, ageLab, R.n_obs, R.n_rats, R.norm_flag);
        catch ME
            fprintf('\n[skip] %s @ %s: %s', vn, ageLab, ME.message);
        end
    end
    PLOT_MWM_CoefAgeSummary(ageResults, ages_present, centers_pres, vn, ...
        'SaveDir', dir_coef);
end
fprintf('\nDone. Figures under %s\n', fig_dir);