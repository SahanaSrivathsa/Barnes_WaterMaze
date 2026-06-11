% Code for Morris Water Maze
% This script runs all analysis for the Morris WaterMaze Code and Plotting
% Assumes spatial files like 'Coh2-NARP-Spatial.xlsx' with Spatial, Probe,
% and Cued sheets. Age, Sex, and Genotype come from TgF344_Aging_AllRats.xlsx
% matched on Barnes ID.

%% GET DATA
clear all;
tg_age_file = 'D:\NARP_Data\RTrack_NARPMale\TgF344_Aging_AllRats.xlsx';
base_dir  = 'D:\NARP_Data\RTrack_NARPMale\Spatial_Sheets';
proc_dir= fullfile(base_dir,'Processed');
if ~exist(proc_dir, "dir")
    mkdir(proc_dir);
end

tg_all_rats = readtable(tg_age_file, 'VariableNamingRule', 'preserve');

% Build Barnes ID lookup from tg_age sheet
barnesIds = str2double(string(tg_all_rats.("Barnes ID")));
ageRaw = tg_all_rats.("Age at Behavior Start");
if isnumeric(ageRaw)
    ages = double(ageRaw(:));
else
    ages = str2double(regexp(string(ageRaw), '\d+\.?\d*', 'match', 'once'));
end
sexes = string(tg_all_rats.Sex);
genotypes = string(tg_all_rats.("APP/WT"));

% Find only unique barnesIds -from when rats were duplicates
[~, ~, grp] = unique(barnesIds);
nPerId = accumarray(grp, 1);
isDup = nPerId(grp) > 1;
doubleRats = tg_all_rats(isDup, :);
if ~isempty(doubleRats)
    doubleRats.BarnesID = barnesIds(isDup);
end
ratMeta = table( ...
    barnesIds(nPerId(grp) == 1), ...
    ages(nPerId(grp) == 1), ...
    sexes(nPerId(grp) == 1), ...
    genotypes(nPerId(grp) == 1), ...
    'VariableNames', {'Animal','Age','Sex','Genotype'});

writetable(doubleRats, fullfile(proc_dir, 'Tg_DuplicateBarnesIDs.csv'));
if ~isempty(doubleRats)
    fprintf('\nWarning: %d aging-sheet row(s) have duplicate Barnes IDs.\n', height(doubleRats));
    fprintf('Saved details to %s\n', fullfile(proc_dir, 'Tg_DuplicateBarnesIDs.csv'));
end
%Initialize
files = dir(fullfile(base_dir,'*.xlsx'));
AllSpatial = table();
AllProbe   = table();
AllVisual  = table();
missingRats = [];
% sheeet names
spatialAliases = {'Spatial'};
probeAliases   = {'Probe', 'Spatial Probe'};
cuedAliases    = {'Cued',  'Visual'};

% headings in each sheet
spHeading = {'Test','Animal','Trial','Duration','Distance','PathEfficiency','Platform_CIPL'};
pbHeading = {'Animal','Trial','Duration','Distance','Q1_Time','Platform_Entries', ...
            'Platform_CIPL','Platform_PathEfficiencyToEntry'};
vsHeading = {'Animal','Trial','Duration','Distance','PathEfficiency', ...
            'Platform_CIPL','Platform_PathEfficiencyToEntry'};
% Loop through and fill tables (doing all three at once )

for i = 1:numel(files)
    fname = fullfile(base_dir, files(i).name);

    matchCoh = regexp(files(i).name,'Coh(\d+)','tokens','once');
    if ~isempty(matchCoh)
        cohortNum = str2double(matchCoh{1});
    else
        cohortNum = NaN;
    end

    spSheet = GetSheetName(fname, spatialAliases);
spT     = readtable(fname,'Sheet',spSheet,'VariableNamingRule','preserve');
spT     = NormalizeMWMHeaders(spT);
[spT, spMiss] = SelectColumnsToHeading(spT, spHeading);

pbSheet = GetSheetName(fname, probeAliases);
pbT     = readtable(fname,'Sheet',pbSheet,'VariableNamingRule','preserve');
pbT     = NormalizeMWMHeaders(pbT);
[pbT, pbMiss] = SelectColumnsToHeading(pbT, pbHeading);

vsSheet = GetSheetName(fname, cuedAliases);
vsT = readtable(fname,'Sheet',vsSheet,'VariableNamingRule','preserve');
vsT   = NormalizeMWMHeaders(vsT);
[vsT, vsMiss] = SelectColumnsToHeading(vsT, vsHeading);

allMiss = [spMiss; pbMiss; vsMiss];
if ~isempty(allMiss)
    fprintf('\n[%s] missing columns (NaN-filled): %s', ...
        files(i).name, strjoin(cellstr(unique(allMiss,'stable')), ', '));
end

    allMiss = [spMiss; pbMiss; vsMiss];
    if ~isempty(allMiss)
        fprintf('\n[%s] expected columns missing (NaN-filled): %s', ...
            files(i).name, strjoin(cellstr(unique(allMiss,'stable')), ', '));
    end
    % Save meta data for all tables at the same time(easier)
    trialTables = {spT, pbT, vsT};
    for s = 1:numel(trialTables)
        T = trialTables{s};
        ids = str2double(string(T.Animal));
        [found, loc] = ismember(ids, ratMeta.Animal);
        T.Age = nan(height(T), 1);
        T.Sex = strings(height(T), 1);
        T.Genotype = strings(height(T), 1);
        T.Age(found) = ratMeta.Age(loc(found));
        T.Sex(found) = ratMeta.Sex(loc(found));
        T.Genotype(found) = ratMeta.Genotype(loc(found));
        missingRats = [missingRats; ids(~found)]; %#ok<AGROW>
        T.Cohort = repmat(cohortNum, height(T), 1);
        trialTables{s} = T;
    end

    AllSpatial = [AllSpatial; trialTables{1}];
    AllProbe   = [AllProbe;   trialTables{2}];
    AllVisual  = [AllVisual;  trialTables{3}];

    fprintf('\nMorris WaterMaze data added for cohort: %d', cohortNum);
end

missingRats = unique(missingRats, 'stable');
if ~isempty(missingRats)
    missingTbl = table(missingRats, 'VariableNames', {'BarnesID'});
    writetable(missingTbl, fullfile(proc_dir, 'Tg_MissingFromAgingSheet.csv'));
    fprintf('\nWarning: %d Barnes ID(s) in spatial files were not found in the aging sheet.\n', numel(missingRats));
end

writetable(AllSpatial, fullfile(proc_dir,'Tg_AllMorrisWaterMazeData_Spatial.csv'));
writetable(AllProbe,   fullfile(proc_dir,'Tg_AllMorrisWaterMazeData_Probe.csv'));
writetable(AllVisual,  fullfile(proc_dir,'Tg_AllMorrisWaterMazeData_Visual.csv'));


