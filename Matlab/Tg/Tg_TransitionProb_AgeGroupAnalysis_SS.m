%% Tg_TransitionProb_AgeGroupAnalysis_SS.m
%  Reproduce the MWM age-summary figures (PLOT_MWM_AgeSubplots +
%  PLOT_MWM_CoefAgeSummary) for TRANSITION-MATRIX metrics on the Tg data:
%    - TransEntropy : Shannon entropy of the within-day transition matrix
%    - SpectralGap  : 1 - |lambda_2| of the within-day transition matrix
%
%  Age-subplots are produced TWICE per metric:
%    Coarse/      one tile per Young / Mid / Old
%    FineMonths/  one tile per populated month bin
%  Coefficient-vs-age (CoefAge/) uses the fine month bins, one fitlme per bin.
%
%  USES
%    Build_Tg_TransitionMetrics, matrixEntropy, AssignAgeBin,
%    STATS_MWM_LME, PLOT_MWM_AgeSubplots, PLOT_MWM_CoefAgeSummary
%
%  SS 2026

%% INPUT PARAMS
clear all;
base_dir   = 'D:\NARP_Data\RTrack_NARPMale';
strat_xlsx = fullfile(base_dir, 'Tg_MWM_results_06-07-2026');   % add .xlsx if readtable errors
proc_dir   = fullfile(base_dir, 'StrategyProcessed');
fig_dir    = fullfile(base_dir, 'Figures', 'AgeSummary_Transition');

% Month bins (same as Analyze_Tg)
ageEdges   = [3.5 4.5 5.5 6.5 7.5 8.5 9.5 10.5 11.5 12.5 13.5 14.5 17.5 22.5];
ageLabels  = ["4mo","5mo","6mo","7mo","8mo","9mo","10mo","11mo","12mo","13mo","14mo","15-16mo","20-21mo"];
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

%% Load + compute per-rat-per-day transition metrics
data1 = readtable(strat_xlsx);
T = Build_Tg_TransitionMetrics(data1, ...
    'AgeEdges', ageEdges, 'AgeLabels', ageLabels, ...
    'YoungMidCut', youngMidCut, 'MidOldCut', midOldCut);
% 'StateSource','column','StateColumn','strategy'  % <- to use RTrack's discrete call
% 'MakeStochastic', true                            % <- proper row-stochastic spectral gap

T = T(~ismissing(T.Group4) & ~ismissing(T.AgeMonth), :);
writetable(T, fullfile(proc_dir, 'Tg_Transition_RatDayMetrics.csv'));

%% DV definitions: {column, ylabel}
dvDefs = { ...
    'TransEntropy', 'Transition entropy (bits)'; ...
    'SpectralGap',  'Spectral gap (1-|\lambda_2|)' };
nDV = size(dvDefs,1);

% Month bins present -> x positions for the coef plot
ages_present = categories(removecats(T.AgeMonth));
[~, ix]      = ismember(string(ages_present), ageLabels);
centers_pres = ageCenters(ix);

%% Loop both metrics through both binnings + coef summary
for v = 1:nDV
    vn   = dvDefs{v,1};
    ylab = dvDefs{v,2};

    % drop rat-days where this metric is undefined (sparse days)
    Tv = T(~isnan(T.(vn)), :);

    % same rows, different tiling factor
    Tcoarse = Tv;
    Tfine   = Tv;
    Tfine.AgeGroup = Tv.AgeMonth;

    % Age-subplots: coarse (Young/Mid/Old) and fine (month bins)
    PLOT_MWM_AgeSubplots(Tcoarse, vn, 'SaveDir', dir_coarse, 'YLabel', ylab);
    PLOT_MWM_AgeSubplots(Tfine,   vn, 'SaveDir', dir_fine,   'YLabel', ylab);

    % Per-month LME - coefficient-vs-age summary
    ageResults = cell(numel(ages_present), 1);
    for a = 1:numel(ages_present)    % LOOP: one fitlme per month bin
        ageLab = ages_present{a};
        Tage   = Tv(Tv.AgeMonth == ageLab, :);
        if height(Tage) < 4, continue; end
        try
            R = STATS_MWM_LME(Tage, vn, 'Factors', "Genotype");
            ageResults{a} = R;
            writetable(R.fixedEffects, fullfile(proc_dir, ...
                sprintf('LME_Transition_%s_%s.csv', vn, strrep(ageLab,' ','_'))));
            fprintf('\n[%s] %-14s @ %-8s | n_obs=%3d n_rats=%2d | resid=%s', ...
                R.model_type, vn, ageLab, R.n_obs, R.n_rats, R.norm_flag);
        catch ME
            fprintf('\n[skip] %s @ %s: %s', vn, ageLab, ME.message);
        end
    end
    PLOT_MWM_CoefAgeSummary(ageResults, ages_present, centers_pres, vn, ...
        'SaveDir', dir_coef);
end
fprintf('\nDone. Figures under %s\n', fig_dir);