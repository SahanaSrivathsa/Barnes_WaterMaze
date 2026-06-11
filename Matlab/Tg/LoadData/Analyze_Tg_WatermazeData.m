%% Analyze_Tg_WatermazeData.m
%  spatial learning and probe trial analysis, per age gorup
% Sex x Genotype is the within-age comparison per age group
% SHOUDL RUN  Load_Tg_WatermazeData.m FIRST

%% SPATIAL TRIALS

%% Load
fname = fullfile(proc_dir, 'Tg_AllMorrisWaterMazeData_Spatial.csv');
T     = readtable(fname);
T     = sortrows(T, ["Animal","Trial"]);
T.Day = floor((T.Trial - 1) / 6) + 1; %6 trials a day      


if ~isnumeric(T.Age), T.Age = str2double(string(T.Age)); end

%Group by month
ageEdges   = [3.5 4.5 5.5 6.5 7.5 8.5 9.5 10.5 11.5 14.5 17.5 22.5];
ageLabels  = ["4mo","5mo","6mo","7mo","8mo","9mo","10mo","11mo","12-13mo","15-16mo","20-21mo"];
ageCenters = (ageEdges(1:end-1) + ageEdges(2:end)) / 2;   % months -> x-axis for coef summary
T.AgeGroup = AssignAgeBin(T.Age, ageEdges, ageLabels);
% Sex x genotype 4 groups
genIn  = string(T.Genotype);
genStr = strings(height(T), 1);
genStr(contains(genIn, "APP", 'IgnoreCase', true)) = "APP";
genStr(genIn == "WT") = "WT";
T.Group4 = categorical(string(T.Sex) + " " + genStr, ["F WT","F APP","M WT","M APP"]);


hasMeta = T.Age > 0 & strlength(string(T.Sex)) > 0 & strlength(genStr) > 0;
varList = {'Platform_CIPL', 'Duration', 'Distance', 'PathEfficiency'};

% Per-rat per-day mean
RatDay = groupsummary(T(hasMeta,:), ...
    {'Animal','Sex','Genotype','Age','AgeGroup','Group4','Day'}, 'mean', varList);
for v = 1:numel(varList)                  % rename mean_X -> X
    old = ['mean_' varList{v}];
    [~, idx] = ismember(old, RatDay.Properties.VariableNames);
    if idx > 0, RatDay.Properties.VariableNames{idx} = varList{v}; end
end
writetable(RatDay, fullfile(proc_dir, 'Tg_MWM_Spatial_RatDayMeans.csv'));
 
%% Per-variable: age-subplot plot + per-age LME + coefficient-vs-age summary
ages_present_sp = categories(removecats(RatDay.AgeGroup));
[~, ix]          = ismember(string(ages_present_sp), ageLabels);
centers_sp       = ageCenters(ix);
%% Age-subplot plots 

for v = 1:numel(varList)
    vn = varList{v};

    PLOT_MWM_AgeSubplots(T, vn, 'SaveDir', proc_dir);

    % stats- separate per age
    ageResults = cell(numel(ages_present_sp), 1);
    for a = 1:numel(ages_present_sp) 
        ageLab = ages_present_sp{a};
        Tage   = RatDay(RatDay.AgeGroup == ageLab, :);
        if height(Tage) < 4, continue; end

        try
            R = STATS_MWM_LME(Tage, vn);
            ageResults{a} = R;
            outName = sprintf('LME_Spatial_%s_%s.csv', vn, strrep(ageLab,' ','_'));
            writetable(R.fixedEffects, fullfile(proc_dir, outName));
            fprintf('\n[%s] %-12s @ %-8s | n_obs=%3d n_rats=%2d | resid=%s', ...
                R.model_type, vn, ageLab, R.n_obs, R.n_rats, R.norm_flag);
        catch ME
            fprintf('\n[skip] %s @ %s: %s', vn, ageLab, ME.message);
        end
    end
    % OVerall fig of coefficients
       PLOT_MWM_CoefAgeSummary(ageResults, ages_present_sp, centers_sp, vn, ...
        'SaveDir', proc_dir);
end
fprintf('\n');

%% PROBE TRIALS

fname  = fullfile(proc_dir, 'Tg_AllMorrisWaterMazeData_Probe.csv');
TProbe = readtable(fname);
if ~isnumeric(TProbe.Age), TProbe.Age = str2double(string(TProbe.Age)); end

% Same age bining as spatial
TProbe.AgeGroup = AssignAgeBin(TProbe.Age, ageEdges, ageLabels);
genIn  = string(TProbe.Genotype);
genStr = strings(height(TProbe), 1);
genStr(contains(genIn, "APP", 'IgnoreCase', true)) = "APP";
genStr(genIn == "WT") = "WT";
TProbe.Group4 = categorical(string(TProbe.Sex) + " " + genStr, ["F WT","F APP","M WT","M APP"]);

varListProbe = {'Duration','Distance','Q1_Time','Platform_Entries', ...
                'Platform_CIPL','Platform_PathEfficiencyToEntry'};

ages_present_pb = categories(removecats(TProbe.AgeGroup));
[~, ix]          = ismember(string(ages_present_pb), ageLabels);
centers_pb       = ageCenters(ix);
for v = 1:numel(varListProbe)
    vn = varListProbe{v};
    if ~ismember(vn, TProbe.Properties.VariableNames), continue; end
    if all(isnan(TProbe.(vn))), continue; end     
    
     ageResults = cell(numel(ages_present_pb), 1);
    for a = 1:numel(ages_present_pb) 
        ageLab = ages_present_pb{a};
        Tage   = TProbe(TProbe.AgeGroup == ageLab, :);
        if height(Tage) < 4, continue; end

        try
            R = STATS_MWM_LME(Tage, vn);
             ageResults{a} = R;
            outName = sprintf('LME_Probe_%s_%s.csv', vn, strrep(ageLab,' ','_'));
            writetable(R.fixedEffects, fullfile(proc_dir, outName));
            fprintf('\n[%s] %-30s @ %-8s | n_obs=%3d | resid=%s', ...
                R.model_type, vn, ageLab, R.n_obs, R.norm_flag);
        catch ME
            fprintf('\n[skip] %s @ %s: %s', vn, ageLab, ME.message);
        end
    end

    PLOT_MWM_CoefAgeSummary(ageResults, ages_present_pb, centers_pb, vn, ...
        'SaveDir', proc_dir);
end
fprintf('\n');