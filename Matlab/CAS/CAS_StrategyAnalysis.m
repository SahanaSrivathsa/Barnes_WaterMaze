base_dir='C:\DATA\WaterMaze\CAS';
processed_dir=fullfile(base_dir,'StrategyProcessed2');
fig_dir=fullfile(base_dir,'Figures');
% Make the dirs if they do not exist
if ~exist(processed_dir,"dir")
    mkdir(processed_dir);
end
if ~exist(fig_dir,"dir")
    mkdir(fig_dir);
end

%Load data from Excel files
data1 = readtable(fullfile(base_dir,'CAS_MWM_results_05-24-2025.xlsx'));  % From RTrack, contains Track_ID, Strategy, Age
data2 = readtable(fullfile(base_dir,'CAS_AllRats_Spatial.csv'));  % From Matlab Analysis contains Test_No, Cohort, Platform_CIPL

%% Analysis Parameters that are constant - mostly colours for groups
age_grps= ["6mo","15mo","23mo"];
perf_grps=["normal","poor","good"];

%clrMap  = {[0.2196,0.5569,0.2353],[0.4157,0.1059,0.6039]}; % green, purple

% Define the strategy column names (adjust names if needed)
strategyNames = {'thigmotaxis','circling','random_path','scanning',...
    'chaining','directed_search','corrected_search','direct_path','perseverance'};
nStrategies = numel(strategyNames);

strategy_titles={'Thigmotaxis','Circling','Random Path','Scanning',...
    'Chaining','Directed Search','Corrected Search','Direct Path','Perseverance'};

% CHECK FOR CONSISTENCY in datasets (Trial no and full values)

%--- For data1: Check that each unique x_TargetID has 24 unique trials ---
grp1 = varfun(@(x) numel(unique(x)), data1, 'GroupingVariables', 'x_TargetID', 'InputVariables', 'x_Trial');
% The new variable name created by varfun is 'Fun_x_Trial'
rats_to_remove = grp1.x_TargetID(grp1.Fun_x_Trial ~= 24);

if ~isempty(rats_to_remove)
    fprintf('Removing the following x_TargetID from data1 (incomplete trials):\n');
    disp(rats_to_remove);
    data1 = data1(~ismember(data1.x_TargetID, rats_to_remove), :);
    data2=data2(~ismember(data2.Animal,rats_to_remove),:);
else
    fprintf('All x_TargetID in data1 have 24 unique trials.\n');
end

%% Plot the strategy differences 

% Get rats of each age group 
young_rats=unique(data2.Animal(data2.Age=="6mo"));
middle_rats=unique(data2.Animal(data2.Age=="15mo"));
old_rats=unique(data2.Animal(data2.Age=="23mo"));

% Get intersection of each group for each age group and plot strategy types

for ii=1:length(age_grps)
    rats=unique(data2.Animal(data2.Age==age_grps(ii)));
    perf
    % Get the rats that intersecct fo reach performance type and young
    good_rats=unique(rats(strcmp(data2.Performance,'good')&&strcmp(data2.Age,age_str)));
    norm_rats=unique(data2.Animal(data2.Performance=='normal'));
    poor_rats=unique(data2.Animal(data2.Performance=='poor'));
end