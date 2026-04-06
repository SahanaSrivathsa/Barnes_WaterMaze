base_dir='/Users/miasponseller/Desktop/Lab/Rtrack/CAS';
processed_dir=fullfile(base_dir,'StrategyProcessed');
fig_dir=fullfile(base_dir,'Figures');
% Make the dirs if they do not exist
if ~exist(processed_dir,"dir")
    mkdir(processed_dir);
end
if ~exist(fig_dir,"dir")
    mkdir(fig_dir);
end

% Load data from Excel files
strat_sheet = readtable(fullfile(base_dir,'CAS_MWM_results_05-24-2025.xlsx'));  % From RTrack, contains Track_ID, Strategy, Age
all_rats_spatial = readtable(fullfile(base_dir,'CAS_AllRats_Spatial.csv'));  % From Python Analysis contains Test_No, Trial_No, Cohort

%% Analysis Parameters that are constant - mostly colours for groups
age_grps= ["6mo", "15mo", "23mo"];
perf_grps=["normal", "poor", "good"];

clrMap  = {[0.2196,0.5569,0.2353],[0.1216, 0.4667, 0.7059], [0.4157,0.1059,0.6039]}; % green, blue, purple

% Colors for strat groups
strat_colors = containers.Map;

strat_colors('NonGoal') = [102, 194, 165]; % teal
strat_colors('Procedural') = [252, 141, 98]; % orange
strat_colors('Allocentric') = [141, 160, 203]; % periwinkle?

% Colors for performance groups
perf_colors = containers.Map;

perf_colors('good')   = [1.000, 0.843, 0.000];  % gold
perf_colors('normal') = [0.439, 0.502, 0.565];  % gray
perf_colors('poor')   = [0.863, 0.078, 0.235];  % crimson

% Strategy groups
strat_groups = {
    {'thigmotaxis', 'circling', 'random_path'}, 'NonGoal';
    {'scanning', 'chaining'}, 'Procedural';
    {'directed_search', 'corected_search', 'direct_path', 'perseverance'}, 'Allocentric'};

% Define the strategy column names (adjust names if needed)
strategyNames = {'thigmotaxis','circling','random_path','scanning',...
    'chaining','directed_search','corrected_search','direct_path','perseverance'};
nStrategies = numel(strategyNames);

strategy_titles = {'Thigmotaxis','Circling','Random Path','Scanning',...
    'Chaining','Directed Search','Corrected Search','Direct Path','Perseverance'};

% CHECK FOR CONSISTENCY in datasets (Trial no and full values)

%--- For strat_sheet: Check that each unique x_TargetID has 24 unique trials ---
grp1 = varfun(@(x) numel(unique(x)), strat_sheet, 'GroupingVariables', 'x_TargetID', 'InputVariables', 'x_Trial');
rats_to_remove = grp1.x_TargetID(grp1.Fun_x_Trial ~= 24);

if ~isempty(rats_to_remove)
    fprintf('Removing the following x_TargetID from strat_sheet (incomplete trials):\n');
    disp(rats_to_remove);
    strat_sheet = strat_sheet(~ismember(strat_sheet.x_TargetID, rats_to_remove),:);
else
    fprintf('All x_TargetID in strat_sheet have 24 unique trials.\n');
end

%--- For all_rats_spatial: Chck that each Animal has 24 unique trials ---
grp2 = varfun(@(x) numel(unique(x)), all_rats_spatial, 'GroupingVariables', 'Animal', 'InputVariables', 'Trial');
animals_incomplete = grp2.Animal(grp2.Fun_Trial ~= 24);

if ~isempty(animals_incomplete)
    fprintf('Removing the following Animals from all_rats_spatial (incomplete trials):\n');
    disp(animals_incomplete);   
    all_rats_spatial = all_rats_spatial(~ismember(all_rats_spatial.Animal, animals_incomplete), :);
else
    fprintf('All Animals in all_rats_spatial have 24 unique trials.\n');
end

%--- For all_rats_spatial: Remove Animals with any NaN in Platform_CIPL ---
% grp_nan = varfun(@(x) any(isnan(x)), all_rats_spatial, 'GroupingVariables', 'Animal', 'InputVariables', 'Platform_CIPL');
% animals_with_nan = grp_nan.Animal(grp_nan.Fun_Platform_CIPL);
% 
% if ~isempty(animals_with_nan)
%     fprintf('Removing the following Animals from all_rats_spatial (NaN in Platform_CIPL):\n');
%     disp(animals_with_nan);
%     all_rats_spatial = all_rats_spatial(~ismember(all_rats_spatial.Animal, animals_with_nan), :);
% else
%     fprintf('No Animals with NaN in Platform_CIPL in all_rats_spatial.\n');
% end

%--- Synchronize strat_sheet and all_rats_spatial: Keep only animals that exist in both datasets ---
commonAnimals = intersect(unique(strat_sheet.x_TargetID), unique(string(all_rats_spatial.Animal)));

strat_sheet = strat_sheet(ismember(strat_sheet.x_TargetID, commonAnimals), :);
all_rats_spatial = all_rats_spatial(ismember(string(all_rats_spatial.Animal), commonAnimals), :);

%% Plot the strategy differences 

% Get rats of each age group 
young_rats=unique(all_rats_spatial.Animal(all_rats_spatial.Age=="6mo"));
middle_rats=unique(all_rats_spatial.Animal(all_rats_spatial.Age=="15mo"));
old_rats=unique(all_rats_spatial.Animal(all_rats_spatial.Age=="23mo"));

% Get intersection of each group for each age group and plot strategy types
for ii = 1:length(age_grps)
    age_str = age_grps(ii); % current age group

    % animals in this age group
    rats = unique(all_rats_spatial.Animal(all_rats_spatial.Age == age_str));
    is_age = all_rats_spatial.Age == age_str;

    % animals for each performance group within this age group
    good_rats = unique(all_rats_spatial.Animal(is_age & all_rats_spatial.Performance == "good"));
    norm_rats = unique(all_rats_spatial.Animal(is_age & all_rats_spatial.Performance == "normal"));
    poor_rats = unique(all_rats_spatial.Animal(is_age & all_rats_spatial.Performance == "poor"));

   % for day=1:numel(unique(data2.Day))
   day=unique(days(ii))
   good_str=strat_sheet.name(strat_sheet.x_TargetID==good_rats && strat_sheet.x_Day==day);

    figure 
    
end

% strat groups
groupLabels = strategyGroups(:,2);       % {'NonGoal','Procedural','Allocentric'}
groupNames  = {'Non‑Goal Oriented','Procedural','Allocentric'}; % for display
n_groups     = numel(groupLabels);

rats = unique(spatial_sheet.x_TargetID);
n_rats = numel(rats);
unique_days = unique(spatial_sheet.Day);
n_days = numel(unique_days);

for a = 1:length(age_grps)
    age_str = age_grps(a);

    for p = 1:length(perf_grps)
        perf_str = perf_str(p);

        % select rats from this age and performance group
        idx = (all_rats_spatial.Age == age_str) & (all_rats_spatial.Performance == perf_str);
        rats = unique(all_rats_spatial.Animal(idx));
        n_rats = numel(rats);


