%% ------------- MUST RUN SECTION FOR ALL FURTHER ANALYSIS ------------- %%
base_dir='/Users/miasponseller/Desktop/Lab/Rtrack/Tg/MatlabFiles';
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
data1 = readtable(fullfile('/Users/miasponseller/Desktop/Lab/Rtrack/Tg/Tg_Results_NoCoh1.xlsx'));  % From RTrack, contains Track_ID, Strategy, Age
data2 = readtable(fullfile('/Users/miasponseller/Desktop/Lab/Rtrack/Tg/MatlabFiles/Processed/AllMorrisWaterMazeData_Spatial.csv'));  % From Matlab Analysis contains Test_No, Cohort, Platform_CIPL

% -----------------------------------------------------------------------
% NEW: Age binning and group label construction
% -----------------------------------------------------------------------
% Bin continuous age into Young / Mid / Old
% Cutoffs: Young <= 9, Mid = 10-15, Old = 16+
% (Modify cutoffs here if needed)
ageCutYoungMid = 9;   % upper bound for Young  (age <= 9  → Young)
ageCutMidOld   = 15;  % upper bound for Mid    (age <= 15 → Mid)

% --- Clean the Age column: ensure it is a numeric double vector,
%     regardless of whether Excel delivered it as cell, string, or double.
if iscell(data1.Age)
    % e.g. {'4','5.5','4 mo', ...}  — strip any non-numeric characters then convert
    data1.Age = cellfun(@(x) str2double(regexp(num2str(x),'[\d.]+','match','once')), ...
                        data1.Age);
elseif isstring(data1.Age) || ischar(data1.Age)
    % string array or char — same approach
    data1.Age = arrayfun(@(x) str2double(regexp(x,'[\d.]+','match','once')), ...
                         string(data1.Age));
end
% data1.Age is now a numeric double column vector.

data1.AgeGroup = repmat("Old", height(data1), 1);
data1.AgeGroup(data1.Age <= ageCutYoungMid) = "Young";
data1.AgeGroup(data1.Age > ageCutYoungMid & data1.Age <= ageCutMidOld) = "Mid";

% Ensure Sex and APP columns are string arrays for consistency
data1.Sex = string(data1.Sex);
data1.APP  = string(data1.APP);

% Build a combined Group label: AgeGroup_Sex_APP  e.g. "Young_M_WT"
% (This is used in sections that want all-combinations colour mapping)
data1.Group = data1.AgeGroup + "_" + data1.Sex + "_" + data1.APP;

% -----------------------------------------------------------------------
% NEW: Derive the full list of groups actually present in the data
%      and assign a colour to each one.
% -----------------------------------------------------------------------
grpList = unique(data1.Group);          % e.g. ["Mid_F_APP/+","Mid_F_WT", ...]
nGrps   = numel(grpList);

% Build a colour palette with nGrps distinct colours.
% We use a perceptually-uniform base set; if there are more groups than
% base colours the colormap will cycle through the full hsv wheel.
baseCols = [
    0.2196, 0.5569, 0.2353;   % green
    0.4157, 0.1059, 0.6039;   % purple
    0.8392, 0.1529, 0.1569;   % red
    0.1216, 0.4706, 0.7059;   % blue
    0.9020, 0.6235, 0.0000;   % amber
    0.0000, 0.6275, 0.6275;   % teal
    0.9176, 0.3176, 0.1412;   % vermilion
    0.3373, 0.7059, 0.9137;   % sky blue
    0.8000, 0.4745, 0.6549;   % pink
    0.4980, 0.4980, 0.4980;   % grey
    0.6275, 0.3216, 0.1765;   % brown
    0.5608, 0.6902, 0.1961;   % lime
];
if nGrps <= size(baseCols,1)
    clrMap = num2cell(baseCols(1:nGrps,:), 2);
else
    % Fall back to evenly-spaced hsv colours
    c = hsv(nGrps);
    clrMap = num2cell(c, 2);
end

% Helper: get colour for a given group label string
getColor = @(g) clrMap{strcmp(grpList, g)};

% -----------------------------------------------------------------------
% comboList is built AFTER data filtering — see block just before %% 1)
% -----------------------------------------------------------------------

% -----------------------------------------------------------------------
% Legacy two-group references kept for any section that still uses them
% (will be overridden inside sections that use flexible grouping)
% -----------------------------------------------------------------------

% Analysis Parameters
strategyNames = {'thigmotaxis','circling','randomPath','scanning',...
    'chaining','directedSearch','correctedPath','directPath','perseverance'};
nStrategies = numel(strategyNames);

strategy_titles={'Thigmotaxis','Circling','Random Path','Scanning',...
    'Chaining','Directed Search','Corrected Search','Direct Path','Perseverance'};

% CHECK FOR CONSISTENCY in datasets (Trial no and full values)

%--- For data1: Check that each unique x_TargetID has 24 unique trials ---
grp1 = varfun(@(x) numel(unique(x)), data1, 'GroupingVariables', 'x_TargetID', 'InputVariables', 'x_Trial');
rats_to_remove = grp1.x_TargetID(grp1.Fun_x_Trial ~= 24);

if ~isempty(rats_to_remove)
    fprintf('Removing the following x_TargetID from data1 (incomplete trials):\n');
    disp(rats_to_remove);
    data1 = data1(~ismember(data1.x_TargetID, rats_to_remove), :);
else
    fprintf('All x_TargetID in data1 have 24 unique trials.\n');
end

%--- For data2: Check that each Animal has 24 unique trials ---
grp2 = varfun(@(x) numel(unique(x)), data2, 'GroupingVariables', 'Animal', 'InputVariables', 'Trial');
animals_incomplete = grp2.Animal(grp2.Fun_Trial ~= 24);

if ~isempty(animals_incomplete)
    fprintf('Removing the following Animals from data2 (incomplete trials):\n');
    disp(animals_incomplete);
    data2 = data2(~ismember(data2.Animal, animals_incomplete), :);
else
    fprintf('All Animals in data2 have 24 unique trials.\n');
end

%--- For data2: Remove Animals with any NaN in Platform_CIPL ---
grp_nan = varfun(@(x) any(isnan(x)), data2, 'GroupingVariables', 'Animal', 'InputVariables', 'Platform_CIPL');
animals_with_nan = grp_nan.Animal(grp_nan.Fun_Platform_CIPL);

if ~isempty(animals_with_nan)
    fprintf('Removing the following Animals from data2 (NaN in Platform_CIPL):\n');
    disp(animals_with_nan);
    data2 = data2(~ismember(data2.Animal, animals_with_nan), :);
else
    fprintf('No Animals with NaN in Platform_CIPL in data2.\n');
end

%--- Synchronize data1 and data2: Keep only animals that exist in both datasets ---
commonAnimals = intersect(unique(data1.x_TargetID), unique(string(data2.Animal)));

data1 = data1(ismember(data1.x_TargetID, commonAnimals), :);
data2 = data2(ismember(string(data2.Animal), commonAnimals), :);

% Get Platform Scores for CIPL
platformScores = nan(height(data1), 1);
 
% Pre-convert data2.Cohort to uppercase string array for comparison
data2CohortStr = upper(string(data2.Cohort));
 
for i = 1:height(data1)
    tID = data1.Track_ID{i};
 
    % Extract cohort label and test/trial numbers from Track_ID
    % Handles: Coh10M_test5  or  Coh10_test5
    tokens = regexp(tID, 'Coh(\d+[A-Za-z]*)_test(\d+)', 'tokens', 'once');
    if isempty(tokens)
        tokens = regexp(tID, 'Coh(\d+)_test(\d+)', 'tokens', 'once');
    end
    if isempty(tokens), continue; end
 
    cohortStr = upper(string(tokens{1}));  % e.g. "10M"
    animalNum = str2double(data1.x_TargetID{i});
    trialNum  = data1.x_Trial(i);
 
    % Match on Cohort string + Animal number + Trial number
    idx2 = find(data2CohortStr == cohortStr & ...
                data2.Animal   == animalNum & ...
                data2.Trial    == trialNum, 1);
 
    if ~isempty(idx2)
        platformScores(i) = data2.Platform_CIPL(idx2);
    end
end

validIdx = ~isnan(platformScores);
data1 = data1(validIdx, :);
platformScores = platformScores(validIdx);

validAnimals = str2double(data1.x_TargetID);
validTrials = data1.x_Trial;

keepIdx = false(height(data2), 1);
for i = 1:length(validAnimals)
    match = data2.Animal == validAnimals(i) & data2.Trial == validTrials(i);
    if any(match)
        keepIdx = keepIdx | match;
    end
end

data2 = data2(keepIdx, :);

% -----------------------------------------------------------------------
% Build comboList AFTER all filtering so masks align with filtered data1
% and platformScores. One entry per unique Sex x Age x APP combination.
% -----------------------------------------------------------------------
[uniqueRatsAll, ia] = unique(data1.x_TargetID);
comboSex = data1.Sex(ia);
comboAge = data1.Age(ia);
comboAPP = data1.APP(ia);

makeComboLabel = @(sex, age, app) sprintf('%s_%g_%s', sex, age, ...
    strrep(strrep(char(app), '/', ''), ' ', ''));

ratComboLabel = arrayfun(@(i) makeComboLabel(comboSex(i), comboAge(i), comboAPP(i)), ...
    (1:numel(uniqueRatsAll))', 'UniformOutput', false);
ratComboLabel = string(ratComboLabel);

[comboLabels, ~, comboIdx] = unique(ratComboLabel);
nCombos = numel(comboLabels);
comboN  = accumarray(comboIdx, 1);

if nCombos <= size(baseCols,1)
    comboCols = num2cell(baseCols(1:nCombos,:), 2);
else
    c = hsv(nCombos); comboCols = num2cell(c, 2);
end

comboList = struct();
for ci = 1:nCombos
    comboList(ci).label   = char(comboLabels(ci));
    comboList(ci).n       = comboN(ci);
    comboList(ci).color   = comboCols{ci};
    % Build mask into the already-filtered data1
    comboList(ci).mask    = ismember(data1.x_TargetID, ...
        uniqueRatsAll(comboIdx == ci));
    comboList(ci).ratMask = (comboIdx == ci);
end

comboTitle = @(ci, extra) sprintf('%s  (n=%d)%s', ...
    comboList(ci).label, comboList(ci).n, ['   ' extra]);

%% 1) CIPL-Strategy Plots - All trials
% -----------------------------------------------------------------------
% NEW: Set groupBy to choose which factor to colour by in this section.
%   Options: 'AgeGroup'  → Young / Mid / Old
%            'Sex'       → M / F  (or whatever values are in data1.Sex)
%            'APP'       → WT / APP/+
%            'Group'     → all combinations (AgeGroup_Sex_APP)
%            'Age'       → individual numeric age values (one colour each)
% -----------------------------------------------------------------------
groupBy_S1 = 'AgeGroup';   % <-- CHANGE THIS PER ANALYSIS

polyOrder = 2;

for s = 1:nStrategies
    stratProb = data1.(strategyNames{s});

    % Get unique levels for the chosen grouping factor
    if strcmp(groupBy_S1, 'Age')
        levels = unique(data1.Age);
        getLvl = @(row) data1.Age(row);
        lvlStr = @(l) sprintf('%.1f mo', l);
        matchFn = @(l) data1.Age == l;
    else
        levels = unique(data1.(groupBy_S1));
        matchFn = @(l) strcmp(data1.(groupBy_S1), l);
        lvlStr  = @(l) char(l);
    end

    % Build a local colour map for the levels in this section
    nLvls_s1 = numel(levels);
    if nLvls_s1 <= size(baseCols,1)
        localClr = num2cell(baseCols(1:nLvls_s1,:), 2);
    else
        c = hsv(nLvls_s1); localClr = num2cell(c,2);
    end

    f = figure; hold on;

    for li = 1:nLvls_s1
        lvl = levels(li);
        if strcmp(groupBy_S1,'Age')
            idx = data1.Age == lvl;
        else
            idx = strcmp(data1.(groupBy_S1), lvl);
        end
        x_vals = platformScores(idx);
        y_vals = stratProb(idx);
        % Filter non-positive for display
        pos = x_vals > 0;
        scatter(x_vals(pos), y_vals(pos), 20, localClr{li}, 'filled',...
            'MarkerFaceAlpha', 0.5, 'DisplayName', lvlStr(lvl));
    end

    title(sprintf('%s: CIPL', strategy_titles{s}));
    xlabel('CIPL Score (m.s)');
    ylabel('Probability of Strategy Use');
    legend('show');
    xlim([0 60]);
    ylim([0 1]);
    set(gca, 'FontSize', 14, 'FontWeight', 'bold', ...
         'LineWidth', 1.5, 'Box', 'off', 'TickDir', 'out');
    hold off;

    saveas(f, fullfile(fig_dir, 'CIPL_Strategy', sprintf('CIPL_%s', strategyNames{s})), 'png');
end

% ---- Per Sex×Age×APP combo figures for Section 1 ---- %
for s = 1:nStrategies
    stratProb = data1.(strategyNames{s});
    for ci = 1:nCombos
        idx     = comboList(ci).mask;
        x_vals  = platformScores(idx);
        y_vals  = stratProb(idx);
        pos     = x_vals > 0;

        f = figure; hold on;
        scatter(x_vals(pos), y_vals(pos), 30, comboList(ci).color, 'filled', ...
            'MarkerFaceAlpha', 0.6);

        title(sprintf('%s: CIPL  —  %s  (n=%d)', ...
            strategy_titles{s}, comboList(ci).label, comboList(ci).n));
        xlabel('CIPL Score (m.s)');
        ylabel('Probability of Strategy Use');
        xlim([0 60]); ylim([0 1]);
        set(gca, 'FontSize', 14, 'FontWeight', 'bold', ...
         'LineWidth', 1.5, 'Box', 'off', 'TickDir', 'out');
        hold off;

        saveas(f, fullfile(fig_dir, 'CIPL_Strategy', ...
            sprintf('CIPL_%s_%s.png', strategyNames{s}, comboList(ci).label)));
        close;
    end
end

%% 2) CIPL vs group Probabilities
% -----------------------------------------------------------------------
% NEW: Set groupBy to choose colouring factor for this section.
%   Same options as Section 1.
% -----------------------------------------------------------------------
groupBy_S2 = 'AgeGroup';   % <-- CHANGE THIS PER ANALYSIS

polyOrder = 2;
strategyGroups = {
    {'thigmotaxis', 'circling', 'randomPath'}, 'Platform-Independent';
    {'scanning', 'chaining'}, 'Procedural';
    {'directedSearch', 'correctedPath', 'directPath','perseverance'}, 'Allocentric'
    };
groupNames = strategyGroups(:, 2);

for s = 1:numel(groupNames)
    currentStrategies = strategyGroups{s, 1};
    groupProb=data1.(currentStrategies{1});
    for ii=2:numel(currentStrategies)
        groupProb=groupProb+data1.(currentStrategies{ii});
    end
    stratProb = groupProb;

    % Build local levels & colours for chosen factor
    if strcmp(groupBy_S2,'Age')
        levels = unique(data1.Age);
        lvlStr = @(l) sprintf('%.1f mo', l);
    else
        levels = unique(data1.(groupBy_S2));
        lvlStr = @(l) char(l);
    end
    nLvls_s2 = numel(levels);
    if nLvls_s2 <= size(baseCols,1)
        localClr = num2cell(baseCols(1:nLvls_s2,:), 2);
    else
        c = hsv(nLvls_s2); localClr = num2cell(c,2);
    end

    f=figure; hold on;

    for li = 1:nLvls_s2
        lvl = levels(li);
        if strcmp(groupBy_S2,'Age')
            idx = data1.Age == lvl;
        else
            idx = strcmp(data1.(groupBy_S2), lvl);
        end
        scatter(platformScores(idx), stratProb(idx), 20, localClr{li}, 'filled',...
            'MarkerFaceAlpha',0.5, 'DisplayName', lvlStr(lvl));
    end

    title(sprintf('%s', groupNames{s}),'FontSize',18,'FontWeight','bold');
    xlabel('CIPL Score (m.s)','FontSize',14,'FontWeight','bold');
    ylabel('Probability','FontSize',14,'FontWeight','bold');
    legend('show');
    xlim([0 60]);
    ylim([0 1]);
    set(gca, 'FontSize', 14, 'FontWeight', 'bold', ...
         'LineWidth', 1.5, 'Box', 'off', 'TickDir', 'out');
    hold off;

    saveas(f, fullfile(fig_dir, 'CIPL_Strategy', sprintf('CIPL_%s',groupNames{s})),'png');
end

% ---- Per Sex×Age×APP combo figures for Section 2 ---- %
% for s = 1:numel(groupNames)
%     currentStrategies = strategyGroups{s, 1};
%     groupProb = data1.(currentStrategies{1});
%     for ii = 2:numel(currentStrategies)
%         groupProb = groupProb + data1.(currentStrategies{ii});
%     end
% 
%     for ci = 1:nCombos
%         idx    = comboList(ci).mask;
%         x_vals = platformScores(idx);
%         y_vals = groupProb(idx);
% 
%         f = figure; hold on;
%         scatter(x_vals, y_vals, 30, comboList(ci).color, 'filled', ...
%             'MarkerFaceAlpha', 0.6);
% 
%         title(sprintf('%s  —  %s  (n=%d)', ...
%             groupNames{s}, comboList(ci).label, comboList(ci).n), ...
%             'FontSize', 16, 'FontWeight', 'bold');
%         xlabel('CIPL Score (m.s)', 'FontSize', 14, 'FontWeight', 'bold');
%         ylabel('Probability',      'FontSize', 14, 'FontWeight', 'bold');
%         xlim([0 60]); ylim([0 1]);
%         set(gca, 'FontSize', 14, 'FontWeight', 'bold', ...
%          'LineWidth', 1.5, 'Box', 'off', 'TickDir', 'out');
%         hold off;
% 
%         saveas(f, fullfile(fig_dir, 'CIPL_Strategy', ...
%             sprintf('CIPL_%s_%s.png', groupNames{s}, comboList(ci).label)));
%         close;
%     end
% end

% ---- Per individual numeric age figures for Section 2 ---- %
% One figure per strategy group: each unique age value is its own
% coloured scatter series so you can compare e.g. 5mo vs 5.5mo vs 6mo.
uniqueAges_S2 = unique(data1.Age);
nAges_S2      = numel(uniqueAges_S2);
if nAges_S2 <= size(baseCols,1)
    ageCols_S2 = num2cell(baseCols(1:nAges_S2,:), 2);
else
    c = hsv(nAges_S2); ageCols_S2 = num2cell(c, 2);
end
 
for s = 1:numel(groupNames)
    currentStrategies_age = strategyGroups{s, 1};
    groupProb_age = data1.(currentStrategies_age{1});
    for ii = 2:numel(currentStrategies_age)
        groupProb_age = groupProb_age + data1.(currentStrategies_age{ii});
    end
 
    f = figure; hold on;
    for ai = 1:nAges_S2
        age_val = uniqueAges_S2(ai);
        idx_age = data1.Age == age_val;
        nRats_age = numel(unique(data1.x_TargetID(idx_age)));
        scatter(platformScores(idx_age), groupProb_age(idx_age), 20, ...
            ageCols_S2{ai}, 'filled', 'MarkerFaceAlpha', 0.5, ...
            'DisplayName', sprintf('%.4g mo (n=%d)', age_val, nRats_age));
    end
 
    title(sprintf('%s: by Age', groupNames{s}), 'FontSize', 18, 'FontWeight', 'bold');
    xlabel('CIPL Score (m.s)', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Probability',      'FontSize', 14, 'FontWeight', 'bold');
    legend('show', 'Location', 'bestoutside');
    xlim([0 60]); ylim([0 1]);
    set(gca, 'FontSize', 14, 'FontWeight', 'bold', ...
        'LineWidth', 1.5, 'Box', 'off', 'TickDir', 'out');
    hold off;
 
    saveas(f, fullfile(fig_dir, 'CIPL_Strategy', ...
        sprintf('CIPL_%s_byAge.png', groupNames{s})));
    close;
end

% ---- Per individual numeric age figures for Section 2 ---- %
% One figure per strategy group: each unique age value is its own
% coloured scatter series so you can compare e.g. 5mo vs 5.5mo vs 6mo.
uniqueAges_S2 = unique(data1.Age);
nAges_S2      = numel(uniqueAges_S2);
if nAges_S2 <= size(baseCols,1)
    ageCols_S2 = num2cell(baseCols(1:nAges_S2,:), 2);
else
    c = hsv(nAges_S2); ageCols_S2 = num2cell(c, 2);
end
 
for s = 1:numel(groupNames)
    currentStrategies_age = strategyGroups{s, 1};
    groupProb_age = data1.(currentStrategies_age{1});
    for ii = 2:numel(currentStrategies_age)
        groupProb_age = groupProb_age + data1.(currentStrategies_age{ii});
    end
 
    for ai = 1:nAges_S2
        age_val   = uniqueAges_S2(ai);
        idx_age   = data1.Age == age_val;
        nRats_age = numel(unique(data1.x_TargetID(idx_age)));
 
        f = figure; hold on;
        scatter(platformScores(idx_age), groupProb_age(idx_age), 30, ...
            ageCols_S2{ai}, 'filled', 'MarkerFaceAlpha', 0.6);
 
        title(sprintf('%s  —  %.4g mo  (n=%d)', groupNames{s}, age_val, nRats_age), ...
            'FontSize', 16, 'FontWeight', 'bold');
        xlabel('CIPL Score (m.s)', 'FontSize', 14, 'FontWeight', 'bold');
        ylabel('Probability',      'FontSize', 14, 'FontWeight', 'bold');
        xlim([0 60]); ylim([0 1]);
        set(gca, 'FontSize', 14, 'FontWeight', 'bold', ...
            'LineWidth', 1.5, 'Box', 'off', 'TickDir', 'out');
        hold off;
 
        % filename: replace '.' with 'p' so e.g. 5.5mo → 5p5mo
        ageLabel = strrep(sprintf('%.4g', age_val), '.', 'p');
        saveas(f, fullfile(fig_dir, 'CIPL_Strategy', ...
            sprintf('CIPL_%s_%smo.png', groupNames{s}, ageLabel)));
        close;
    end
end


%% SEX + GENOTYPE 
% ---- Per Sex figures for Section 2 ---- %
% One figure per strategy group per sex (M / F), all ages combined.
uniqueSexes_S2 = unique(data1.Sex);
nSexes_S2      = numel(uniqueSexes_S2);
if nSexes_S2 <= size(baseCols,1)
    sexCols_S2 = num2cell(baseCols(1:nSexes_S2,:), 2);
else
    c = hsv(nSexes_S2); sexCols_S2 = num2cell(c, 2);
end
 
for s = 1:numel(groupNames)
    currentStrategies_sex = strategyGroups{s, 1};
    groupProb_sex = data1.(currentStrategies_sex{1});
    for ii = 2:numel(currentStrategies_sex)
        groupProb_sex = groupProb_sex + data1.(currentStrategies_sex{ii});
    end
 
    for si = 1:nSexes_S2
        sex_val   = uniqueSexes_S2(si);
        idx_sex   = strcmp(data1.Sex, sex_val);
        nRats_sex = numel(unique(data1.x_TargetID(idx_sex)));
 
        f = figure; hold on;
        scatter(platformScores(idx_sex), groupProb_sex(idx_sex), 30, ...
            sexCols_S2{si}, 'filled', 'MarkerFaceAlpha', 0.6);
 
        title(sprintf('%s  —  %s  (n=%d)', groupNames{s}, sex_val, nRats_sex), ...
            'FontSize', 16, 'FontWeight', 'bold');
        xlabel('CIPL Score (m.s)', 'FontSize', 14, 'FontWeight', 'bold');
        ylabel('Probability',      'FontSize', 14, 'FontWeight', 'bold');
        xlim([0 60]); ylim([0 1]);
        set(gca, 'FontSize', 14, 'FontWeight', 'bold', ...
            'LineWidth', 1.5, 'Box', 'off', 'TickDir', 'out');
        hold off;
 
        saveas(f, fullfile(fig_dir, 'CIPL_Strategy', ...
            sprintf('CIPL_%s_%s.png', groupNames{s}, sex_val)));
        close;
    end
end
 
% ---- Per Genotype figures for Section 2 ---- %
% One figure per strategy group per genotype (WT / APP/+), all ages combined.
uniqueAPP_S2 = unique(data1.APP);
nAPP_S2      = numel(uniqueAPP_S2);
if nAPP_S2 <= size(baseCols,1)
    appCols_S2 = num2cell(baseCols(1:nAPP_S2,:), 2);
else
    c = hsv(nAPP_S2); appCols_S2 = num2cell(c, 2);
end
 
for s = 1:numel(groupNames)
    currentStrategies_app = strategyGroups{s, 1};
    groupProb_app = data1.(currentStrategies_app{1});
    for ii = 2:numel(currentStrategies_app)
        groupProb_app = groupProb_app + data1.(currentStrategies_app{ii});
    end
 
    for gi = 1:nAPP_S2
        app_val   = uniqueAPP_S2(gi);
        idx_app   = strcmp(data1.APP, app_val);
        nRats_app = numel(unique(data1.x_TargetID(idx_app)));
        % safe filename: replace / and spaces
        appLabel  = strrep(strrep(char(app_val), '/', ''), ' ', '');
 
        f = figure; hold on;
        scatter(platformScores(idx_app), groupProb_app(idx_app), 30, ...
            appCols_S2{gi}, 'filled', 'MarkerFaceAlpha', 0.6);
 
        title(sprintf('%s  —  %s  (n=%d)', groupNames{s}, char(app_val), nRats_app), ...
            'FontSize', 16, 'FontWeight', 'bold');
        xlabel('CIPL Score (m.s)', 'FontSize', 14, 'FontWeight', 'bold');
        ylabel('Probability',      'FontSize', 14, 'FontWeight', 'bold');
        xlim([0 60]); ylim([0 1]);
        set(gca, 'FontSize', 14, 'FontWeight', 'bold', ...
            'LineWidth', 1.5, 'Box', 'off', 'TickDir', 'out');
        hold off;
 
        saveas(f, fullfile(fig_dir, 'CIPL_Strategy', ...
            sprintf('CIPL_%s_%s.png', groupNames{s}, appLabel)));
        close;
    end
end
 
% ---- Per Sex×Genotype figures for Section 2 ---- %
% One figure per strategy group per Sex-Genotype combo (M-WT, M-APP+, F-WT, F-APP+).
% Build unique Sex×Genotype pairs actually present in the data.
[uniqueRats_sg, ia_sg] = unique(data1.x_TargetID);
sgSex = data1.Sex(ia_sg);
sgAPP = data1.APP(ia_sg);
sgLabel = sgSex + "-" + sgAPP;  % e.g. "M-WT", "F-APP/+"
[uniqueSG, ~, sgIdx] = unique(sgLabel);
nSG = numel(uniqueSG);
if nSG <= size(baseCols,1)
    sgCols = num2cell(baseCols(1:nSG,:), 2);
else
    c = hsv(nSG); sgCols = num2cell(c, 2);
end
 
for s = 1:numel(groupNames)
    currentStrategies_sg = strategyGroups{s, 1};
    groupProb_sg = data1.(currentStrategies_sg{1});
    for ii = 2:numel(currentStrategies_sg)
        groupProb_sg = groupProb_sg + data1.(currentStrategies_sg{ii});
    end
 
    for sgi = 1:nSG
        sg_val    = uniqueSG(sgi);
        % mask into data1 rows where Sex-APP matches this combo
        idx_sg    = ismember(data1.x_TargetID, uniqueRats_sg(sgIdx == sgi));
        nRats_sg  = sum(sgIdx == sgi);
        % safe filename label: replace / and spaces
        sgFileLabel = strrep(strrep(char(sg_val), '/', ''), ' ', '');
 
        f = figure; hold on;
        scatter(platformScores(idx_sg), groupProb_sg(idx_sg), 30, ...
            sgCols{sgi}, 'filled', 'MarkerFaceAlpha', 0.6);
 
        title(sprintf('%s  —  %s  (n=%d)', groupNames{s}, char(sg_val), nRats_sg), ...
            'FontSize', 16, 'FontWeight', 'bold');
        xlabel('CIPL Score (m.s)', 'FontSize', 14, 'FontWeight', 'bold');
        ylabel('Probability',      'FontSize', 14, 'FontWeight', 'bold');
        xlim([0 60]); ylim([0 1]);
        set(gca, 'FontSize', 14, 'FontWeight', 'bold', ...
            'LineWidth', 1.5, 'Box', 'off', 'TickDir', 'out');
        hold off;
 
        saveas(f, fullfile(fig_dir, 'CIPL_Strategy', ...
            sprintf('CIPL_%s_%s.png', groupNames{s}, sgFileLabel)));
        close;
    end
end

%% Sex and Genotype - one plot

% SEX
for s = 1:numel(groupNames)
    currentStrategies_sex = strategyGroups{s, 1};
    groupProb_sex = data1.(currentStrategies_sex{1});
    for ii = 2:numel(currentStrategies_sex)
        groupProb_sex = groupProb_sex + data1.(currentStrategies_sex{ii});
    end

    f = figure; hold on;

    for si = 1:nSexes_S2
        sex_val   = uniqueSexes_S2(si);
        idx_sex   = strcmp(data1.Sex, sex_val);
        nRats_sex = numel(unique(data1.x_TargetID(idx_sex)));

        scatter(platformScores(idx_sex), groupProb_sex(idx_sex), 30, ...
            sexCols_S2{si}, 'filled', 'MarkerFaceAlpha', 0.6);
    end

    title(sprintf('%s — Sex Combined', groupNames{s}), ...
        'FontSize', 16, 'FontWeight', 'bold');
    xlabel('CIPL Score (m.s)', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Probability',      'FontSize', 14, 'FontWeight', 'bold');
    xlim([0 60]); ylim([0 1]);

    set(gca, 'FontSize', 14, 'FontWeight', 'bold', ...
        'LineWidth', 1.5, 'Box', 'off', 'TickDir', 'out');

    legend(uniqueSexes_S2, 'Location', 'best');
    hold off;

    saveas(f, fullfile(fig_dir, 'CIPL_Strategy', ...
        sprintf('CIPL_%s_SexCombined.png', groupNames{s})));
    close;
end

% GENOTYPE
for s = 1:numel(groupNames)
    currentStrategies_app = strategyGroups{s, 1};
    groupProb_app = data1.(currentStrategies_app{1});
    for ii = 2:numel(currentStrategies_app)
        groupProb_app = groupProb_app + data1.(currentStrategies_app{ii});
    end

    f = figure; hold on;

    for gi = 1:nAPP_S2
        app_val   = uniqueAPP_S2(gi);
        idx_app   = strcmp(data1.APP, app_val);
        nRats_app = numel(unique(data1.x_TargetID(idx_app)));

        scatter(platformScores(idx_app), groupProb_app(idx_app), 30, ...
            appCols_S2{gi}, 'filled', 'MarkerFaceAlpha', 0.6);
    end

    title(sprintf('%s — Genotype Combined', groupNames{s}), ...
        'FontSize', 16, 'FontWeight', 'bold');
    xlabel('CIPL Score (m.s)', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Probability',      'FontSize', 14, 'FontWeight', 'bold');
    xlim([0 60]); ylim([0 1]);

    set(gca, 'FontSize', 14, 'FontWeight', 'bold', ...
        'LineWidth', 1.5, 'Box', 'off', 'TickDir', 'out');

    legend(uniqueAPP_S2, 'Location', 'best');
    hold off;

    saveas(f, fullfile(fig_dir, 'CIPL_Strategy', ...
        sprintf('CIPL_%s_APPCombined.png', groupNames{s})));
    close;
end

% SEX & GENOTYPE
for s = 1:numel(groupNames)
    currentStrategies_sg = strategyGroups{s, 1};
    groupProb_sg = data1.(currentStrategies_sg{1});
    for ii = 2:numel(currentStrategies_sg)
        groupProb_sg = groupProb_sg + data1.(currentStrategies_sg{ii});
    end

    f = figure; hold on;

    for sgi = 1:nSG
        sg_val = uniqueSG(sgi);

        idx_sg = ismember(data1.x_TargetID, uniqueRats_sg(sgIdx == sgi));
        nRats_sg = sum(sgIdx == sgi);

        scatter(platformScores(idx_sg), groupProb_sg(idx_sg), 30, ...
            sgCols{sgi}, 'filled', 'MarkerFaceAlpha', 0.6);
    end

    title(sprintf('%s — Sex × Genotype Combined', groupNames{s}), ...
        'FontSize', 16, 'FontWeight', 'bold');
    xlabel('CIPL Score (m.s)', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Probability',      'FontSize', 14, 'FontWeight', 'bold');
    xlim([0 60]); ylim([0 1]);

    set(gca, 'FontSize', 14, 'FontWeight', 'bold', ...
        'LineWidth', 1.5, 'Box', 'off', 'TickDir', 'out');

    legend(uniqueSG, 'Location', 'best');
    hold off;

    saveas(f, fullfile(fig_dir, 'CIPL_Strategy', ...
        sprintf('CIPL_%s_SexGenotypeCombined.png', groupNames{s})));
    close;
end

%% 3) Strategy-Day Plots (individual and group)
% -----------------------------------------------------------------------
% NEW: Set groupBy to choose colouring/grouping factor for ANOVA and plots.
%   Options: 'AgeGroup', 'Sex', 'APP', 'Group'
%   Note: 'Age' (numeric) is not suitable here as it would give one
%         "group" per animal — use 'AgeGroup' instead.
% -----------------------------------------------------------------------
groupBy_S3 = 'AgeGroup';   % <-- CHANGE THIS PER ANALYSIS

apaTbl = table( ...
    strings(0,1) , strings(0,1) , ...
    zeros (0,1)  , ...
    zeros (0,1)  , zeros (0,1) , ...
    zeros (0,1)  , zeros (0,1)  , ...
    zeros (0,1)  , ...
    'VariableNames', ...
    {'Strategy','Effect','df','SS','MS','F','p','eta2'});

% Derive levels and colours for chosen factor
levels_S3 = unique(data1.(groupBy_S3));
nLvls_S3  = numel(levels_S3);
if nLvls_S3 <= size(baseCols,1)
    clrMap_S3 = num2cell(baseCols(1:nLvls_S3,:), 2);
else
    c = hsv(nLvls_S3); clrMap_S3 = num2cell(c,2);
end

for s = 1:nStrategies

    %-------------- Individual Rat Plot --------------%
    stratProb = data1.(strategyNames{s});
    uniqueDays = unique(data1.Day);
    uniqueRats = unique(data1.x_TargetID);

    f1 = figure;
    hold on;

    mean_strat = {};
    validRatIdx = 0;

    % Track first-plot flag per level for legend
    firstPlot = true(1, nLvls_S3);

    for r = 1:numel(uniqueRats)
        ratID = uniqueRats{r};
        isCurrentRat = strcmp(data1.x_TargetID, ratID);

        % Determine this rat's level for the chosen grouping factor
        ratLvl = data1.(groupBy_S3){find(isCurrentRat,1)};
        lvlIdx = find(strcmp(levels_S3, ratLvl), 1);
        color  = clrMap_S3{lvlIdx};

        meanUse = nan(1, numel(uniqueDays));
        for d = 1:numel(uniqueDays)
            idx = isCurrentRat & (data1.Day == uniqueDays(d));
            if any(idx)
                meanUse(d) = mean(stratProb(idx));
            end
        end

        if any(isnan(meanUse))
            fprintf('Missing days for %s in strategy %s\n', string(ratID), strategyNames{s});
            continue;
        end

        validRatIdx = validRatIdx + 1;
        mean_strat{validRatIdx, 1} = ratID;
        mean_strat{validRatIdx, 2} = ratLvl;   % store the group label
        mean_strat(validRatIdx, 3:6) = num2cell(meanUse);

        hSwarm = swarmchart(uniqueDays, meanUse, 20, 'filled', ...
            'MarkerFaceColor', color, 'MarkerFaceAlpha', 0.5);
        hSwarm.Annotation.LegendInformation.IconDisplayStyle = 'off';

        if firstPlot(lvlIdx)
            plot(uniqueDays, meanUse, '-', 'Color', color, 'LineWidth', 1, ...
                'DisplayName', char(levels_S3(lvlIdx)));
            firstPlot(lvlIdx) = false;
        else
            plot(uniqueDays, meanUse, '-', 'Color', color, 'LineWidth', 1, ...
                'HandleVisibility','off');
        end
    end

    title(strategy_titles{s});
    xlabel('Day');
    ylabel('Probability of Strategy Use');
    xticks(uniqueDays);
    legend('Location', 'Northeast');
    ylim([0 0.5]);
    xlim([0.75 4.25]);
    pubify_figure_axis_robust(14,14);
    hold off;
    saveas(f1, fullfile(fig_dir,'StrategyUse', sprintf('IndividualRats_%s', strategyNames{s})), 'png');
    close;

    % Convert cell array to table for ANOVA
    % Column 2 is now whatever groupBy_S3 is labelled
    mean_strat_table = cell2table(mean_strat, 'VariableNames', ...
        {'RatID', 'Age', 'Day1', 'Day2', 'Day3', 'Day4'});
    % NOTE: the 'Age' column here actually contains the groupBy_S3 values.
    % The runMixedANOVA helper uses 'Age' as the between-subjects factor
    % column name — rename if your helper uses a different name.

    anovaResults = runMixedANOVA(mean_strat_table, {'Day1','Day2','Day3','Day4'});
    postHocResults = runTukeyPostHocMixed(mean_strat_table, {'Day1','Day2','Day3','Day4'});

    writetable(anovaResults, fullfile(processed_dir, sprintf('Anova_%s.csv', strategyNames{s})));
    writetable(postHocResults, fullfile(processed_dir, sprintf('PostHoc_Tukey_%s.csv', strategyNames{s})));

    rowNames = string(anovaResults.Properties.RowNames);

    if ~strcmpi(strategyNames{s},'perseverance')
        wanted   = ["Age","(Intercept):Day","Age:Day"];
        tidyLab  = ["Group","Day","Group×Day"];   % updated label to reflect flexible grouping

        keepIdx  = find(ismember(rowNames, wanted));
        if isempty(keepIdx),  warning('%s has no mixed-ANOVA rows',strategyNames{s});  end
        L        = numel(keepIdx);

        SS   = anovaResults.SumSq(keepIdx);
        df   = anovaResults.DF   (keepIdx);
        MS   = anovaResults.MeanSq(keepIdx);
        Fval = anovaResults.F    (keepIdx);
        pVal = anovaResults.pValue(keepIdx);

        eta2 = nan(L,1);
        for j = 1:L
            k    = keepIdx(j);
            eRow = find(startsWith(rowNames(k+1:end),"Error"),1,'first') + k;
            if ~isempty(eRow)
                eta2(j) = SS(j) / ( SS(j) + anovaResults.SumSq(eRow) );
            end
        end

        tmp = table( ...
            repmat(string(strategy_titles{s}), L,1)  , ...
            tidyLab(ismember(wanted,rowNames(keepIdx))).' , ...
            df   , SS , MS , Fval , pVal , eta2 , ...
            'VariableNames',{'Strategy','Effect','df','SS','MS','F','p','eta2'});

        apaTbl = [apaTbl ; tmp];
    end

    %-------------- Overall Bar Plot - Mean Per Rat Per Day --------------%
    dataMeanGrps = cell(nLvls_S3, numel(uniqueDays));
    for li = 1:nLvls_S3
        lvl = levels_S3(li);
        for d = 1:numel(uniqueDays)
            colName = sprintf('Day%d', d);
            dataMeanGrps{li,d} = mean_strat_table{strcmp(mean_strat_table.Age, lvl), colName};
        end
    end

    % For the bar plot helper we still pass Young/Old style — here we pass
    % the first two levels if there are only two, otherwise we only plot
    % the first two and print a warning. For >2 groups, see note below.
    if nLvls_S3 == 2
        meanGrp1 = cellfun(@mean, dataMeanGrps(1,:));
        semGrp1  = cellfun(@(x) std(x)/sqrt(numel(x)), dataMeanGrps(1,:));
        meanGrp2 = cellfun(@mean, dataMeanGrps(2,:));
        semGrp2  = cellfun(@(x) std(x)/sqrt(numel(x)), dataMeanGrps(2,:));

        f3=figure;
        plot_bar_sem_WaterMaze(uniqueDays, meanGrp1, semGrp1, meanGrp2, semGrp2, ...
            dataMeanGrps(1,:)', dataMeanGrps(2,:)', postHocResults);
        ylim([0 0.55]);
        pubify_figure_axis_robust(16,16);
        hold off;
        exportgraphics(f3, fullfile(fig_dir, 'StrategyUse',...
            sprintf('1_RatMeanPerDay_%s.png', strategyNames{s})),'Resolution', 450 );
        close
    else
        % NEW: For >2 groups, build a simple multi-line plot instead.
        % plot_bar_sem_WaterMaze only supports 2 groups; extend it if needed.
        f3 = figure; hold on;
        for li = 1:nLvls_S3
            meanVals = cellfun(@mean, dataMeanGrps(li,:));
            semVals  = cellfun(@(x) std(x)/sqrt(numel(x)), dataMeanGrps(li,:));
            errorbar(uniqueDays, meanVals, semVals, '-o', ...
                'Color', clrMap_S3{li}, 'LineWidth', 2, ...
                'DisplayName', char(levels_S3(li)));
        end
        xlabel('Day'); ylabel('Probability of Strategy Use');
        title(strategy_titles{s});
        legend('Location','Northeast');
        ylim([0 0.55]);
        xticks(uniqueDays);
        pubify_figure_axis_robust(16,16);
        hold off;
        exportgraphics(f3, fullfile(fig_dir, 'StrategyUse',...
            sprintf('1_RatMeanPerDay_%s.png', strategyNames{s})),'Resolution', 450 );
        close
    end
end

writetable(apaTbl, fullfile(processed_dir,'ANOVA_APA_AllStrategies.csv'));
disp(apaTbl)
saveTablePNG_APA(apaTbl, ...
    fullfile(fig_dir,'ANOVA_APA_AllStrategies.png'), ...
    'Title','Repeated Measures ANOVA for Strategies');

% ---- Per Sex×Age×APP combo figures for Section 3 ---- %
% One figure per strategy per combo: mean strategy use per day as a line plot.
uniqueDays_c3 = unique(data1.Day);
for s = 1:nStrategies
    stratProb = data1.(strategyNames{s});
    for ci = 1:nCombos
        idx     = comboList(ci).mask;
        nRats_c = comboList(ci).n;

        % Compute per-rat daily means
        rats_c  = unique(data1.x_TargetID(idx));
        ratDayMeans = nan(numel(rats_c), numel(uniqueDays_c3));
        for r = 1:numel(rats_c)
            rIdx = strcmp(data1.x_TargetID, rats_c{r});
            for d = 1:numel(uniqueDays_c3)
                dIdx = rIdx & (data1.Day == uniqueDays_c3(d));
                if any(dIdx)
                    ratDayMeans(r,d) = mean(stratProb(dIdx));
                end
            end
        end

        meanVals = mean(ratDayMeans, 1, 'omitnan');
        if nRats_c > 1
            semVals = std(ratDayMeans, 0, 1, 'omitnan') ./ sqrt(nRats_c);
        else
            semVals = zeros(size(meanVals));
        end

        f = figure; hold on;
        % Individual rat lines (thin, transparent)
        for r = 1:nRats_c
            plot(uniqueDays_c3, ratDayMeans(r,:), '-o', ...
                'Color', [comboList(ci).color 0.25], ...
                'LineWidth', 1, 'MarkerSize', 3, ...
                'MarkerFaceColor', comboList(ci).color, ...
                'MarkerEdgeColor', 'none', ...
                'HandleVisibility', 'off');
        end
        % Group mean ± SEM
        errorbar(uniqueDays_c3, meanVals, semVals, '-o', ...
            'Color', comboList(ci).color, 'LineWidth', 2.5, ...
            'MarkerFaceColor', comboList(ci).color, ...
            'MarkerEdgeColor', 'none', 'MarkerSize', 7, ...
            'DisplayName', comboList(ci).label);

        title(sprintf('%s  —  %s  (n=%d)', ...
            strategy_titles{s}, comboList(ci).label, nRats_c));
        xlabel('Day'); ylabel('Probability of Strategy Use');
        xticks(uniqueDays_c3);
        xlim([0.75 4.25]); ylim([0 0.55]);
        pubify_figure_axis_robust(14,14);
        hold off;

        saveas(f, fullfile(fig_dir, 'StrategyUse', ...
            sprintf('1_RatMeanPerDay_%s_%s.png', strategyNames{s}, comboList(ci).label)));
        close;
    end
end

%% 4) Strategy by Group - Non-spatial/Procedural/Allocentric
% -----------------------------------------------------------------------
% NEW: Set groupBy to choose colouring/grouping factor for this section.
% -----------------------------------------------------------------------
groupBy_S4 = 'AgeGroup';   % <-- CHANGE THIS PER ANALYSIS

strategyGroups = {
    {'thigmotaxis', 'circling', 'randomPath'}, 'Non-Goal Oriented';
    {'scanning', 'chaining'}, 'Procedural';
    {'directedSearch', 'correctedPath', 'directPath','perseverance'}, 'Allocentric'
    };
groupNames = strategyGroups(:, 2);
nGroups = numel(groupNames);

% Derive levels and colours for chosen factor
levels_S4 = unique(data1.(groupBy_S4));
nLvls_S4  = numel(levels_S4);
if nLvls_S4 <= size(baseCols,1)
    clrMap_S4 = num2cell(baseCols(1:nLvls_S4,:), 2);
else
    c = hsv(nLvls_S4); clrMap_S4 = num2cell(c,2);
end

groupStrat=cell(1,3);
uniqueDays = unique(data1.Day);

for g = 1:nGroups
    f1 = figure;
    hold on;

    currentStrategies = strategyGroups{g, 1};
    groupProb=data1.(currentStrategies{1});
    for ii=2:numel(currentStrategies)
        groupProb=groupProb+data1.(currentStrategies{ii});
    end
    groupStrat{g}=groupProb;

    mean_strat = {};
    validRatIdx = 0;
    firstPlot = true(1, nLvls_S4);

    uniqueRats = unique(data1.x_TargetID);
    for r = 1:numel(uniqueRats)
        ratID = uniqueRats{r};
        isCurrentRat = strcmp(data1.x_TargetID, ratID);

        ratLvl = data1.(groupBy_S4){find(isCurrentRat,1)};
        lvlIdx = find(strcmp(levels_S4, ratLvl), 1);
        color  = clrMap_S4{lvlIdx};

        meanUse = nan(1, numel(uniqueDays));
        for d = 1:numel(uniqueDays)
            dayIdx = isCurrentRat & (data1.Day == uniqueDays(d));
            if any(dayIdx)
                meanUse(d) = mean(groupProb(dayIdx));
            end
        end

        if any(isnan(meanUse))
            fprintf('Missing days for %s in group %s\n', string(ratID), groupNames{g});
            continue;
        end

        validRatIdx = validRatIdx + 1;
        mean_strat{validRatIdx, 1} = ratID;
        mean_strat{validRatIdx, 2} = ratLvl;
        mean_strat(validRatIdx, 3:6) = num2cell(meanUse);

        hSwarm = swarmchart(uniqueDays, meanUse, 20, 'filled', ...
            'MarkerFaceColor', color, 'MarkerFaceAlpha', 0.5);
        hSwarm.Annotation.LegendInformation.IconDisplayStyle = 'off';

        if firstPlot(lvlIdx)
            plot(uniqueDays, meanUse, '-', 'Color', color, 'LineWidth', 1, ...
                'DisplayName', char(levels_S4(lvlIdx)));
            firstPlot(lvlIdx) = false;
        else
            plot(uniqueDays, meanUse, '-', 'Color', color, 'LineWidth', 1, ...
                'HandleVisibility', 'off');
        end
    end

    title(groupNames{g});
    xlabel('Day');
    ylabel('Probability of Strategy Use');
    xticks(uniqueDays);
    legend('Location', 'northeastoutside');
    ylim([0 1]);
    xlim([0.75 4.25]);
    pubify_figure_axis_robust(14, 14);
    hold off;

    saveas(f1, fullfile(fig_dir, 'StrategyUse',sprintf('IndividualRats_%s', groupNames{g})), 'png');
    close;

    mean_strat_table = cell2table(mean_strat, 'VariableNames', ...
        {'RatID', 'Age', 'Day1', 'Day2', 'Day3', 'Day4'});

    anovaResults = runMixedANOVA(mean_strat_table, {'Day1', 'Day2', 'Day3', 'Day4'});
    postHocResults = runTukeyPostHocMixed(mean_strat_table, {'Day1', 'Day2', 'Day3', 'Day4'});

    writetable(anovaResults, fullfile(processed_dir, sprintf('Anova_%s.csv', groupNames{g})));
    writetable(postHocResults, fullfile(processed_dir, sprintf('PostHoc_Tukey_%s.csv', groupNames{g})));

    %-------------- Figure for Mean Per Rat Per Day --------------%
    dataMeanGrps = cell(nLvls_S4, numel(uniqueDays));
    for li = 1:nLvls_S4
        lvl = levels_S4(li);
        for d = 1:numel(uniqueDays)
            colName = sprintf('Day%d', d);
            dataMeanGrps{li,d} = mean_strat_table{strcmp(mean_strat_table.Age, lvl), colName};
        end
    end

    if nLvls_S4 == 2
        meanGrp1 = cellfun(@mean, dataMeanGrps(1,:));
        semGrp1  = cellfun(@(x) std(x)/sqrt(numel(x)), dataMeanGrps(1,:));
        meanGrp2 = cellfun(@mean, dataMeanGrps(2,:));
        semGrp2  = cellfun(@(x) std(x)/sqrt(numel(x)), dataMeanGrps(2,:));

        f3 = figure;
        plot_bar_sem_WaterMaze(uniqueDays, meanGrp1, semGrp1, meanGrp2, semGrp2, ...
            dataMeanGrps(1,:)', dataMeanGrps(2,:)', postHocResults);
        title(sprintf('%s: Rat Mean Per Day', groupNames{g}), 'FontSize', 16);
        xlabel('Day');
        ylabel('Probability of Strategy Use');
        ylim([0 1.25]);
        pubify_figure_axis_robust(14, 14);
        hold off;
        saveas(f3, fullfile(fig_dir, 'StrategyUse',sprintf('RatMeanPerDay_%s', groupNames{g})), 'png');
        close;
    else
        f3 = figure; hold on;
        for li = 1:nLvls_S4
            meanVals = cellfun(@mean, dataMeanGrps(li,:));
            semVals  = cellfun(@(x) std(x)/sqrt(numel(x)), dataMeanGrps(li,:));
            errorbar(uniqueDays, meanVals, semVals, '-o', ...
                'Color', clrMap_S4{li}, 'LineWidth', 2, ...
                'DisplayName', char(levels_S4(li)));
        end
        xlabel('Day'); ylabel('Probability of Strategy Use');
        title(sprintf('%s: Rat Mean Per Day', groupNames{g}), 'FontSize', 16);
        legend('Location','northeastoutside');
        ylim([0 1.25]);
        xticks(uniqueDays);
        pubify_figure_axis_robust(14, 14);
        hold off;
        saveas(f3, fullfile(fig_dir, 'StrategyUse',sprintf('RatMeanPerDay_%s', groupNames{g})), 'png');
        close;
    end
end

% ---- Per Sex×Age×APP combo figures for Section 4 ---- %
% One figure per strategy group per combo.
uniqueDays_c4 = unique(data1.Day);
strategyGroups_c4 = {
    {'thigmotaxis', 'circling', 'randomPath'}, 'Non-Goal Oriented';
    {'scanning', 'chaining'}, 'Procedural';
    {'directedSearch', 'correctedPath', 'directPath','perseverance'}, 'Allocentric'
    };
for g = 1:size(strategyGroups_c4,1)
    currentStrategies_c4 = strategyGroups_c4{g,1};
    groupProb_c4 = data1.(currentStrategies_c4{1});
    for ii = 2:numel(currentStrategies_c4)
        groupProb_c4 = groupProb_c4 + data1.(currentStrategies_c4{ii});
    end
    gName_c4 = strategyGroups_c4{g,2};

    for ci = 1:nCombos
        idx     = comboList(ci).mask;
        nRats_c = comboList(ci).n;
        rats_c  = unique(data1.x_TargetID(idx));

        ratDayMeans = nan(numel(rats_c), numel(uniqueDays_c4));
        for r = 1:numel(rats_c)
            rIdx = strcmp(data1.x_TargetID, rats_c{r});
            for d = 1:numel(uniqueDays_c4)
                dIdx = rIdx & (data1.Day == uniqueDays_c4(d));
                if any(dIdx)
                    ratDayMeans(r,d) = mean(groupProb_c4(dIdx));
                end
            end
        end

        meanVals = mean(ratDayMeans, 1, 'omitnan');
        if nRats_c > 1
            semVals = std(ratDayMeans, 0, 1, 'omitnan') ./ sqrt(nRats_c);
        else
            semVals = zeros(size(meanVals));
        end

        f = figure; hold on;
        for r = 1:nRats_c
            plot(uniqueDays_c4, ratDayMeans(r,:), '-o', ...
                'Color', [comboList(ci).color 0.25], 'LineWidth', 1, ...
                'MarkerSize', 3, 'MarkerFaceColor', comboList(ci).color, ...
                'MarkerEdgeColor', 'none', 'HandleVisibility', 'off');
        end
        errorbar(uniqueDays_c4, meanVals, semVals, '-o', ...
            'Color', comboList(ci).color, 'LineWidth', 2.5, ...
            'MarkerFaceColor', comboList(ci).color, ...
            'MarkerEdgeColor', 'none', 'MarkerSize', 7);

        title(sprintf('%s  —  %s  (n=%d)', gName_c4, comboList(ci).label, nRats_c), ...
            'FontSize', 14);
        xlabel('Day'); ylabel('Probability of Strategy Use');
        xticks(uniqueDays_c4);
        xlim([0.75 4.25]); ylim([0 1.25]);
        pubify_figure_axis_robust(14,14);
        hold off;

        saveas(f, fullfile(fig_dir, 'StrategyUse', ...
            sprintf('RatMeanPerDay_%s_%s.png', gName_c4, comboList(ci).label)));
        close;
    end
end

%% 5) Entropy Calculation
% -----------------------------------------------------------------------
% NEW: Set groupBy to choose colouring/grouping factor for this section.
% -----------------------------------------------------------------------
groupBy_S5 = 'AgeGroup';   % <-- CHANGE THIS PER ANALYSIS

% Derive levels and colours for chosen factor
levels_S5 = unique(data1.(groupBy_S5));
nLvls_S5  = numel(levels_S5);
if nLvls_S5 <= size(baseCols,1)
    clrMap_S5 = num2cell(baseCols(1:nLvls_S5,:), 2);
else
    c = hsv(nLvls_S5); clrMap_S5 = num2cell(c,2);
end

%------------------ Compute Entropy for Each Trial ------------------%
pMat = data1{:, 11:19};
entropy_vals = -sum(pMat .* log2(pMat + eps), 2);
data1.entropy = entropy_vals;

uniqueRats = unique(data1.x_TargetID);
uniqueDays = unique(data1.Day);
nRats = numel(uniqueRats);

meanDayEntropy = nan(nRats, numel(uniqueDays));
AgeCell = cell(nRats,1);   % stores the groupBy_S5 label per rat

for r = 1:nRats
    ratID = uniqueRats{r};
    idxRat = strcmp(data1.x_TargetID, ratID);
    AgeCell{r} = data1.(groupBy_S5){find(idxRat,1)};
    for d = 1:numel(uniqueDays)
        idxDay = idxRat & (data1.Day == uniqueDays(d));
        meanDayEntropy(r,d) = mean(data1.entropy(idxDay));
    end
end

entropyTable = table(uniqueRats, AgeCell, meanDayEntropy(:,1), meanDayEntropy(:,2), ...
    meanDayEntropy(:,3), meanDayEntropy(:,4), ...
    'VariableNames', {'RatID', 'Age', 'Day1', 'Day2', 'Day3', 'Day4'});

anovaResults = runMixedANOVA(entropyTable, {'Day1','Day2','Day3','Day4'});
postHocResults = runTukeyPostHocMixed(entropyTable, {'Day1','Day2','Day3','Day4'});
writetable(anovaResults, fullfile(processed_dir, 'Entropy_ANOVA.csv'));
writetable(postHocResults, fullfile(processed_dir, 'Entropy_PostHoc_Tukey.csv'));

%------------------Entropy for Each Trial for Each Rat ------------------%
fe1=figure;
hold on;
uniqueRats = unique(data1.x_TargetID);
jitterAmount = 0.1;

for r = 1:length(uniqueRats)
    ratID = uniqueRats{r};
    idx = strcmp(data1.x_TargetID, ratID);

    trialNumbers = data1.x_Trial(idx);
    ratLvl = data1.(groupBy_S5){find(idx,1)};
    lvlIdx = find(strcmp(levels_S5, ratLvl),1);
    thisColor = clrMap_S5{lvlIdx};

    % Jitter direction offset by level index so groups spread apart
    jitterOffset = (lvlIdx - (nLvls_S5+1)/2) * 0.15;
    jitteredTrials = trialNumbers + jitterOffset + ...
        (rand(size(trialNumbers))-0.5)*jitterAmount;

    entropy_r = data1.entropy(idx);
    scatter(jitteredTrials, entropy_r, 36, 'MarkerFaceColor', thisColor, ...
        'MarkerEdgeColor', 'none');
end
xlabel('Trial');
ylabel('Entropy');
title('Entropy: All Trials and All Rats');
pubify_figure_axis_robust(14,14);
hold off;
saveas(fe1, fullfile(fig_dir, 'Entropy','EntropyAllTrials'), 'png');
close;

%------------------ Mean Entropy per Rat per Day ------------------%
meanEntropy = zeros(numel(uniqueDays), nLvls_S5);
semEntropy  = zeros(numel(uniqueDays), nLvls_S5);
dataEntropyPerLvl = cell(numel(uniqueDays), nLvls_S5);

for d = 1:numel(uniqueDays)
    for li = 1:nLvls_S5
        lvl = levels_S5(li);
        idx = strcmp(entropyTable.Age, lvl);
        vals = meanDayEntropy(idx, d);
        meanEntropy(d,li) = mean(vals,'omitnan');
        semEntropy(d,li)  = std(vals,'omitnan') / sqrt(sum(~isnan(vals)));
        dataEntropyPerLvl{d,li} = vals;
    end
end

%---------------- PLOT PER RAT PER DAY MEAN ENTROPY ---------------------------%
if nLvls_S5 == 2
    fe2=figure; hold on;
    plot_bar_sem_WaterMaze(uniqueDays, meanEntropy(:,1), semEntropy(:,1), ...
        meanEntropy(:,2), semEntropy(:,2), ...
        dataEntropyPerLvl(:,1), dataEntropyPerLvl(:,2), postHocResults);
    title('Mean Entropy');
    xlabel('Day'); ylabel('Entropy');
    pubify_figure_axis_robust(14,14);
    xticks(uniqueDays);
    hold off;
    saveas(fe2, fullfile(fig_dir, 'Entropy', 'Entropy_MeanRat'), 'png');
    close;
else
    fe2 = figure; hold on;
    for li = 1:nLvls_S5
        errorbar(uniqueDays, meanEntropy(:,li), semEntropy(:,li), '-o', ...
            'Color', clrMap_S5{li}, 'LineWidth', 2, ...
            'DisplayName', char(levels_S5(li)));
    end
    title('Mean Entropy'); xlabel('Day'); ylabel('Entropy');
    legend('Location','Northeast');
    xticks(uniqueDays);
    pubify_figure_axis_robust(14,14);
    hold off;
    saveas(fe2, fullfile(fig_dir, 'Entropy', 'Entropy_MeanRat'), 'png');
    close;
end

% Per-rat trajectory plot with group mean lines
fe3 = figure; clf; hold on;
alphaIndiv = 0.2;
lwIndiv    = 1.5;
lwMean     = 6;

for r = 1:nRats
    lvlIdx   = find(strcmp(levels_S5, AgeCell{r}), 1);
    thisCol  = [clrMap_S5{lvlIdx}, alphaIndiv];
    plot(uniqueDays, meanDayEntropy(r,:), '-o', ...
        'Color', thisCol, 'LineWidth', lwIndiv, ...
        'MarkerSize', 4, 'MarkerFaceColor', thisCol(1:3), ...
        'MarkerEdgeColor', 'none');
end

for li = 1:nLvls_S5
    ageMask = strcmp(AgeCell, levels_S5(li));
    grpMean = mean(meanDayEntropy(ageMask,:), 1, 'omitnan');
    plot(uniqueDays, grpMean, '-o', ...
        'Color', clrMap_S5{li}, 'LineWidth', lwMean, ...
        'MarkerSize', 6, 'MarkerFaceColor', clrMap_S5{li}, ...
        'DisplayName', char(levels_S5(li)));
end

xlabel('Day'); ylabel('Entropy');
title('Entropy Across All 8 Strategy Types');
xticks(uniqueDays);
legend('Location','Northeast');
pubify_figure_axis_robust(14,14);
hold off;
saveas(fe3, fullfile(fig_dir, 'Entropy', 'Entropy_All8_MeanRat'), 'png');
close(fe3);

% ---- Per Sex×Age×APP combo figures for Section 5 ---- %
% One figure per combo: entropy mean per day.
uniqueDays_c5 = unique(data1.Day);
% Use meanDayEntropy and AgeCell already built above, but index by combo mask.
% We need per-rat entropy aligned to comboList, so recompute using comboList masks.
for ci = 1:nCombos
    idx     = comboList(ci).mask;
    nRats_c = comboList(ci).n;
    rats_c  = unique(data1.x_TargetID(idx));

    ratEntr = nan(numel(rats_c), numel(uniqueDays_c5));
    for r = 1:numel(rats_c)
        rIdx = strcmp(data1.x_TargetID, rats_c{r});
        for d = 1:numel(uniqueDays_c5)
            dIdx = rIdx & (data1.Day == uniqueDays_c5(d));
            if any(dIdx)
                ratEntr(r,d) = mean(data1.entropy(dIdx));
            end
        end
    end

    meanE = mean(ratEntr, 1, 'omitnan');
    if nRats_c > 1
        semE = std(ratEntr, 0, 1, 'omitnan') ./ sqrt(nRats_c);
    else
        semE = zeros(size(meanE));
    end

    fe_c = figure; hold on;
    for r = 1:nRats_c
        plot(uniqueDays_c5, ratEntr(r,:), '-o', ...
            'Color', [comboList(ci).color 0.25], 'LineWidth', 1, ...
            'MarkerSize', 3, 'MarkerFaceColor', comboList(ci).color, ...
            'MarkerEdgeColor', 'none', 'HandleVisibility', 'off');
    end
    errorbar(uniqueDays_c5, meanE, semE, '-o', ...
        'Color', comboList(ci).color, 'LineWidth', 2.5, ...
        'MarkerFaceColor', comboList(ci).color, ...
        'MarkerEdgeColor', 'none', 'MarkerSize', 7);

    title(sprintf('Entropy  —  %s  (n=%d)', comboList(ci).label, nRats_c));
    xlabel('Day'); ylabel('Entropy');
    xticks(uniqueDays_c5);
    pubify_figure_axis_robust(14,14);
    hold off;

    saveas(fe_c, fullfile(fig_dir, 'Entropy', ...
        sprintf('Entropy_MeanRat_%s.png', comboList(ci).label)));
    close;
end

%% 6) Group Strategy Based Entropy
% -----------------------------------------------------------------------
% NEW: Set groupBy to choose colouring/grouping factor for this section.
% -----------------------------------------------------------------------
groupBy_S6 = 'AgeGroup';   % <-- CHANGE THIS PER ANALYSIS

levels_S6 = unique(data1.(groupBy_S6));
nLvls_S6  = numel(levels_S6);
if nLvls_S6 <= size(baseCols,1)
    clrMap_S6 = num2cell(baseCols(1:nLvls_S6,:), 2);
else
    c = hsv(nLvls_S6); clrMap_S6 = num2cell(c,2);
end

%------------------ Compute Group Strategy Entropy ------------------%
pMat = [groupStrat{1}';groupStrat{2}';groupStrat{3}']';
entropy_vals = -sum(pMat .* log2(pMat + eps), 2);
data1.entropy = entropy_vals;

uniqueRats = unique(data1.x_TargetID);
uniqueDays = unique(data1.Day);
nRats = numel(uniqueRats);

meanDayEntropy = nan(nRats, numel(uniqueDays));
AgeCell = cell(nRats,1);

for r = 1:nRats
    ratID = uniqueRats{r};
    idxRat = strcmp(data1.x_TargetID, ratID);
    AgeCell{r} = data1.(groupBy_S6){find(idxRat,1)};
    for d = 1:numel(uniqueDays)
        idxDay = idxRat & (data1.Day == uniqueDays(d));
        meanDayEntropy(r,d) = mean(data1.entropy(idxDay));
    end
end

entropyTable = table(uniqueRats, AgeCell, meanDayEntropy(:,1), meanDayEntropy(:,2), ...
    meanDayEntropy(:,3), meanDayEntropy(:,4), ...
    'VariableNames', {'RatID', 'Age', 'Day1', 'Day2', 'Day3', 'Day4'});

anovaResults = runMixedANOVA(entropyTable, {'Day1','Day2','Day3','Day4'});
disp(anovaResults)
postHocResults = runTukeyPostHocMixed(entropyTable, {'Day1','Day2','Day3','Day4'});
writetable(anovaResults, fullfile(processed_dir, 'GroupStrat_Entropy_ANOVA.csv'));
writetable(postHocResults, fullfile(processed_dir, 'GroupStrat_Entropy_PostHoc_Tukey.csv'));

%------------------Entropy for Each Trial for Each Rat ------------------%
fe1=figure; hold on;
uniqueRats = unique(data1.x_TargetID);
jitterAmount = 0.1;

for r = 1:length(uniqueRats)
    ratID = uniqueRats{r};
    idx = strcmp(data1.x_TargetID, ratID);
    trialNumbers = data1.x_Trial(idx);
    ratLvl = data1.(groupBy_S6){find(idx,1)};
    lvlIdx = find(strcmp(levels_S6, ratLvl),1);
    thisColor = clrMap_S6{lvlIdx};
    jitterOffset = (lvlIdx - (nLvls_S6+1)/2) * 0.15;
    jitteredTrials = trialNumbers + jitterOffset + ...
        (rand(size(trialNumbers))-0.5)*jitterAmount;
    entropy_r = data1.entropy(idx);
    scatter(jitteredTrials, entropy_r, 36, 'MarkerFaceColor', thisColor, ...
        'MarkerEdgeColor', 'none');
end
xlabel('Trial'); ylabel('Entropy');
title('Entropy: All Trials - Grouped By Strategy Type');
pubify_figure_axis_robust(14,14);
hold off;
saveas(fe1, fullfile(fig_dir, 'Entropy', 'GroupStrat_EntropyAllTrials'), 'png');
close;

%------------------ Mean Entropy per Rat per Day ------------------%
meanEntropy = zeros(numel(uniqueDays), nLvls_S6);
semEntropy  = zeros(numel(uniqueDays), nLvls_S6);
dataEntropyPerLvl = cell(numel(uniqueDays), nLvls_S6);

for d = 1:numel(uniqueDays)
    for li = 1:nLvls_S6
        lvl = levels_S6(li);
        idx = strcmp(entropyTable.Age, lvl);
        vals = meanDayEntropy(idx, d);
        meanEntropy(d,li) = mean(vals,'omitnan');
        semEntropy(d,li)  = std(vals,'omitnan') / sqrt(sum(~isnan(vals)));
        dataEntropyPerLvl{d,li} = vals;
    end
end

if nLvls_S6 == 2
    fe2=figure; hold on;
    plot_bar_sem_WaterMaze(uniqueDays, meanEntropy(:,1), semEntropy(:,1), ...
        meanEntropy(:,2), semEntropy(:,2), ...
        dataEntropyPerLvl(:,1), dataEntropyPerLvl(:,2), postHocResults);
    title('Entropy (Rats Day Means) - Grouped By Strategy Type');
    xlabel('Day'); ylabel('Entropy');
    pubify_figure_axis_robust(14,14);
    xticks(uniqueDays);
    hold off;
    saveas(fe2, fullfile(fig_dir,'Entropy', 'GroupStrat_Entropy_MeanRat'), 'png');
    close;
else
    fe2 = figure; hold on;
    for li = 1:nLvls_S6
        errorbar(uniqueDays, meanEntropy(:,li), semEntropy(:,li), '-o', ...
            'Color', clrMap_S6{li}, 'LineWidth', 2, ...
            'DisplayName', char(levels_S6(li)));
    end
    title('Entropy (Rats Day Means) - Grouped By Strategy Type');
    xlabel('Day'); ylabel('Entropy');
    legend('Location','Northeast');
    xticks(uniqueDays);
    pubify_figure_axis_robust(14,14);
    hold off;
    saveas(fe2, fullfile(fig_dir,'Entropy', 'GroupStrat_Entropy_MeanRat'), 'png');
    close;
end

% Per-rat line plot
fe3 = figure; hold on;
nRats = size(meanDayEntropy,1);
for r = 1:nRats
    lvlIdx   = find(strcmp(levels_S6, AgeCell{r}),1);
    thisColor = clrMap_S6{lvlIdx};
    plot(uniqueDays, meanDayEntropy(r,:), '-o', 'Color', thisColor, 'LineWidth', 1.5);
end
xlabel('Day'); ylabel('Entropy');
title('Entropy Per Rat');
pubify_figure_axis_robust(14,14);
xticks(uniqueDays);
hold off;
saveas(fe3, fullfile(fig_dir, 'Entropy', 'GroupStrat_Entropy_RatChanges'), 'png');
close;

fe4=figure; hold on;
dataEntropyYoung = cell(numel(uniqueDays),1);
dataEntropyOld   = cell(numel(uniqueDays),1);
for d = 1:numel(uniqueDays)
    dataEntropyYoung{d}=entropy_vals(uniqueDays(d) & strcmp(data1.Age,'young'));
    dataEntropyOld{d}=entropy_vals(data1.Day == uniqueDays(d) & strcmp(data1.Age,'old'));
end
plot_bar_sem_WaterMaze(uniqueDays, meanEntropy(:,1), semEntropy(:,1), ...
    meanEntropy(:,2), semEntropy(:,2), dataEntropyYoung, dataEntropyOld, postHocResults);
title(' Entropy (All Trials) - Grouped By Strategy Type');
xlabel('Day'); ylabel('Entropy');
pubify_figure_axis_robust(14,14);
xticks(uniqueDays)
hold off;
saveas(fe4, fullfile(fig_dir, 'Entropy', 'GroupStrat_Entropy_Mean_AllTrials'), 'png');
close;

% ---- Per Sex×Age×APP combo figures for Section 6 ---- %
% Uses data1.entropy which was just overwritten with group-strategy entropy.
uniqueDays_c6 = unique(data1.Day);
for ci = 1:nCombos
    idx     = comboList(ci).mask;
    nRats_c = comboList(ci).n;
    rats_c  = unique(data1.x_TargetID(idx));

    ratEntr6 = nan(numel(rats_c), numel(uniqueDays_c6));
    for r = 1:numel(rats_c)
        rIdx = strcmp(data1.x_TargetID, rats_c{r});
        for d = 1:numel(uniqueDays_c6)
            dIdx = rIdx & (data1.Day == uniqueDays_c6(d));
            if any(dIdx)
                ratEntr6(r,d) = mean(data1.entropy(dIdx));
            end
        end
    end

    meanE6 = mean(ratEntr6, 1, 'omitnan');
    if nRats_c > 1
        semE6 = std(ratEntr6, 0, 1, 'omitnan') ./ sqrt(nRats_c);
    else
        semE6 = zeros(size(meanE6));
    end

    fe_c6 = figure; hold on;
    for r = 1:nRats_c
        plot(uniqueDays_c6, ratEntr6(r,:), '-o', ...
            'Color', [comboList(ci).color 0.25], 'LineWidth', 1, ...
            'MarkerSize', 3, 'MarkerFaceColor', comboList(ci).color, ...
            'MarkerEdgeColor', 'none', 'HandleVisibility', 'off');
    end
    errorbar(uniqueDays_c6, meanE6, semE6, '-o', ...
        'Color', comboList(ci).color, 'LineWidth', 2.5, ...
        'MarkerFaceColor', comboList(ci).color, ...
        'MarkerEdgeColor', 'none', 'MarkerSize', 7);

    title(sprintf('Group Strategy Entropy  —  %s  (n=%d)', comboList(ci).label, nRats_c));
    xlabel('Day'); ylabel('Entropy');
    xticks(uniqueDays_c6);
    pubify_figure_axis_robust(14,14);
    hold off;

    saveas(fe_c6, fullfile(fig_dir, 'Entropy', ...
        sprintf('GroupStrat_Entropy_MeanRat_%s.png', comboList(ci).label)));
    close;
end

%% Change platform groups -test
data1.pd(data1.pd==4)=3;
data1.pd(data1.pd==2)=1;

%% 7) Individual Strategy x Platform Group
% -----------------------------------------------------------------------
% NEW: Set groupBy to choose colouring/grouping factor for this section.
% -----------------------------------------------------------------------
groupBy_S7 = 'AgeGroup';   % <-- CHANGE THIS PER ANALYSIS

levels_S7 = unique(data1.(groupBy_S7));
nLvls_S7  = numel(levels_S7);

% Generate shades from light→dark for each level (replaces youngShades/oldShades)
% Each row of platformShades{li} is one shade for that level
platformShades = cell(1, nLvls_S7);
for li = 1:nLvls_S7
    if li <= size(baseCols,1)
        baseC = baseCols(li,:);
    else
        baseC = hsv(nLvls_S7); baseC = baseC(li,:);
    end
    blendLvls = linspace(0.6, 0, 4)';
    platformShades{li} = blendLvls .* 1 + (1-blendLvls) .* baseC;
end

% Signature colours per level for sigstar
sigColors_lvl = cell(1, nLvls_S7);
for li = 1:nLvls_S7
    if li <= size(baseCols,1)
        sigColors_lvl{li} = baseCols(li,:);
    else
        c = hsv(nLvls_S7); sigColors_lvl{li} = c(li,:);
    end
end
clrBothSig = [0.3 0.3 0.3];

uniqueDays=unique(data1.Day);
uniqueRats=unique(data1.x_TargetID);
uniquePlatforms=unique(data1.pd);
uniquePlatforms = uniquePlatforms(~isnan(uniquePlatforms));

for s = 1:numel(strategyNames)
    stratName  = strategyNames{s};
    stratTitle = strategy_titles{s};

    nPlat = numel(uniquePlatforms);
    measureVars = arrayfun(@(x) sprintf('P%d', x), uniquePlatforms, 'UniformOutput', false);

    master_data = {};

    for d = 1:numel(uniqueDays)
        dVal = uniqueDays(d);
        dayData = data1(data1.Day == dVal, :);

        for r = 1:numel(uniqueRats)
            ratID = uniqueRats{r};
            isCurrentRat = strcmp(dayData.x_TargetID, ratID);
            if ~any(isCurrentRat), continue; end

            % NEW: use groupBy_S7 label instead of hardcoded Age
            ratLvl = dayData.(groupBy_S7){find(isCurrentRat,1)};
            if iscell(ratLvl), ratLvl = ratLvl{1}; end

            meanUse = nan(1, nPlat);
            for p = 1:nPlat
                pVal = uniquePlatforms(p);
                idx = isCurrentRat & (dayData.pd == pVal);
                if any(idx)
                    meanUse(p) = mean(dayData.(stratName)(idx));
                end
            end

            row = {ratID, ratLvl, dVal};
            for p = 1:nPlat, row = [row, {meanUse(p)}]; end
            master_data(end+1,:) = row;
        end
    end

    measureVars = arrayfun(@(x) sprintf('P%d', x), uniquePlatforms, 'UniformOutput', false);
    measureVars = measureVars(:)';
    varNames = [{'RatID','Age','Day'}, measureVars];
    master_table = cell2table(master_data, 'VariableNames', varNames);

    all_anova = table(); all_tukey = table();
    for d = 1:numel(uniqueDays)
        dVal = uniqueDays(d);
        day_table = master_table(master_table.Day == dVal, :);
        anova_tbl = day_table(:, [{'Age'}, measureVars]);
        anovaResults    = runMixedANOVA(anova_tbl, measureVars,'Platform');
        postHocResults  = runTukeyPostHocMixed(anova_tbl, measureVars,'Platform');
        anovaResults.Day  = repmat(dVal, height(anovaResults), 1);
        postHocResults.Day = repmat(dVal, height(postHocResults), 1);
        anovaResults.Properties.RowNames   = {};
        postHocResults.Properties.RowNames = {};
        all_anova = [all_anova; anovaResults];
        all_tukey = [all_tukey; postHocResults];
    end

    writetable(all_anova, fullfile(processed_dir, sprintf('3%s_PlatformGroups_ANOVA.csv', stratTitle)));
    writetable(all_tukey, fullfile(processed_dir, sprintf('3%s_PlatformGroups_Tukey.csv', stratTitle)));

    % Summary stats per Day × Platform × Level
    nDays = numel(uniqueDays);
    % meanGrp{li} = nDays × nPlat matrix for level li
    meanGrp  = cell(1, nLvls_S7);
    semGrp   = cell(1, nLvls_S7);
    dataGrp  = cell(nLvls_S7, nDays, nPlat);  % raw per-rat values

    for li = 1:nLvls_S7
        meanGrp{li} = nan(nDays, nPlat);
        semGrp{li}  = nan(nDays, nPlat);
    end

    for d = 1:nDays
        dVal = uniqueDays(d);
        day_tbl = master_table(master_table.Day == dVal, :);
        for p = 1:nPlat
            colName = measureVars{p};
            for li = 1:nLvls_S7
                lvl = levels_S7(li);
                vals = day_tbl{strcmp(day_tbl.Age, lvl), colName};
                meanGrp{li}(d,p) = mean(vals,'omitnan');
                semGrp{li}(d,p)  = std(vals,'omitnan')/sqrt(sum(~isnan(vals)));
                dataGrp{li,d,p}  = vals;
            end
        end
    end

    % Build grouped bar matrix: columns ordered as [Lvl1P1, Lvl2P1, ..., Lvl1P2, Lvl2P2, ...]
    % i.e. for each platform, one bar per level
    nCols = nLvls_S7 * nPlat;
    meanMatrix = nan(nDays, nCols);
    semMatrix  = nan(nDays, nCols);
    for p = 1:nPlat
        for li = 1:nLvls_S7
            col = (p-1)*nLvls_S7 + li;
            meanMatrix(:, col) = meanGrp{li}(:, p);
            semMatrix(:, col)  = semGrp{li}(:, p);
        end
    end

    fH = figure('Name', stratTitle, 'Color', [1 1 1]);
    sgtitle(stratTitle, 'FontSize', 14, 'FontWeight', 'bold');
    hold on;
    hBar = bar(meanMatrix, 'grouped');
    set(gca, 'XTick', 1:nDays, 'XTickLabel', arrayfun(@num2str, uniqueDays, 'UniformOutput', false));
    xlabel('Day'); ylabel('Mean Strategy Usage');

    for p = 1:nPlat
        for li = 1:nLvls_S7
            col = (p-1)*nLvls_S7 + li;
            shadeIdx = min(p, size(platformShades{li},1));
            hBar(col).FaceColor = platformShades{li}(shadeIdx,:);
            hBar(col).FaceAlpha = 0.75;
        end
    end

    drawnow;
    for c = 1:nCols
        xVals = hBar(c).XEndPoints;
        yVals = meanMatrix(:, c);
        eVals = semMatrix(:, c);
        errorbar(xVals, yVals, eVals, 'k', 'LineStyle', 'none', 'LineWidth', 1.2);
    end

    jitterAmount = 0.05;
    for d = 1:nDays
        for p = 1:nPlat
            for li = 1:nLvls_S7
                col = (p-1)*nLvls_S7 + li;
                pts = dataGrp{li,d,p};
                if isempty(pts), continue; end
                xCenter = hBar(col).XEndPoints(d);
                xJitter = xCenter + (rand(size(pts))-0.5)*jitterAmount;
                scatter(xJitter, pts, 12, 'MarkerFaceColor', hBar(col).FaceColor, ...
                    'MarkerEdgeColor','none', 'MarkerFaceAlpha',1);
            end
        end
    end

    % Significance markers from post hoc
    sigPairs_all = {}; sigPvals_all = []; sigColors_all = {};
    for iRow = 1:height(all_tukey)
        day_val = all_tukey.Day(iRow);
        dayIdx  = find(uniqueDays == day_val);
        if isempty(dayIdx), continue; end
        p1   = all_tukey.Platform_1(iRow);
        p2   = all_tukey.Platform_2(iRow);
        pVal = all_tukey.pValue(iRow);
        if pVal >= 0.05, continue; end
        compAge = all_tukey.Age{iRow};

        % Find matching level
        lvlIdx_sig = find(strcmp(levels_S7, compAge), 1);
        if ~isempty(lvlIdx_sig)
            col1 = (p1-1)*nLvls_S7 + lvlIdx_sig;
            col2 = (p2-1)*nLvls_S7 + lvlIdx_sig;
            x1 = hBar(col1).XEndPoints(dayIdx);
            x2 = hBar(col2).XEndPoints(dayIdx);
            sigColor = sigColors_lvl{lvlIdx_sig};
        else
            % 'both' or unrecognised — span all levels for both platforms
            x_candidates = [];
            for li = 1:nLvls_S7
                x_candidates(end+1) = hBar((p1-1)*nLvls_S7+li).XEndPoints(dayIdx);
                x_candidates(end+1) = hBar((p2-1)*nLvls_S7+li).XEndPoints(dayIdx);
            end
            x1 = min(x_candidates); x2 = max(x_candidates);
            sigColor = clrBothSig;
        end
        sigPairs_all{end+1} = [x1, x2];
        sigPvals_all(end+1) = pVal;
        sigColors_all{end+1} = sigColor;
    end
    if ~isempty(sigPairs_all)
        hS = sigstar(sigPairs_all, sigPvals_all);
        for k = 1:size(hS,1)
            set(hS(k), 'Color', sigColors_all{k});
        end
    end

    hold off;
    pubify_figure_axis_robust(14,14);
    saveas(fH, fullfile(fig_dir,'Platform', sprintf('3Platform_%s.png', stratTitle)));
    close(fH);
end

% ---- Per Sex×Age×APP combo figures for Section 7 ---- %
% One figure per strategy per combo: mean strategy use per platform per day.
uniqueDays_c7 = unique(data1.Day);
uniquePlatforms_c7 = unique(data1.pd);
uniquePlatforms_c7 = uniquePlatforms_c7(~isnan(uniquePlatforms_c7));
nPlat_c7 = numel(uniquePlatforms_c7);
nDays_c7 = numel(uniqueDays_c7);

for s = 1:numel(strategyNames)
    stratName_c7  = strategyNames{s};
    stratTitle_c7 = strategy_titles{s};

    for ci = 1:nCombos
        idx     = comboList(ci).mask;
        nRats_c = comboList(ci).n;
        rats_c  = unique(data1.x_TargetID(idx));

        % Mean per rat per day per platform
        meanByDay = nan(nDays_c7, nPlat_c7);
        semByDay  = nan(nDays_c7, nPlat_c7);

        for d = 1:nDays_c7
            dVal  = uniqueDays_c7(d);
            dMask = idx & (data1.Day == dVal);
            for p = 1:nPlat_c7
                pVal  = uniquePlatforms_c7(p);
                pMask = dMask & (data1.pd == pVal);
                % Per-rat means for this day×platform
                ratVals = nan(numel(rats_c),1);
                for r = 1:numel(rats_c)
                    rIdx = strcmp(data1.x_TargetID, rats_c{r}) & pMask;
                    if any(rIdx)
                        ratVals(r) = mean(data1.(stratName_c7)(rIdx));
                    end
                end
                meanByDay(d,p) = mean(ratVals,'omitnan');
                if nRats_c > 1
                    semByDay(d,p) = std(ratVals,'omitnan') ./ sqrt(sum(~isnan(ratVals)));
                else
                    semByDay(d,p) = 0;
                end
            end
        end

        fH_c7 = figure; hold on;
        hBar_c7 = bar(meanByDay, 'grouped');
        % Shade each platform bar using darker shades of the combo colour
        blendL = linspace(0.6, 0, nPlat_c7)';
        shades_c7 = blendL .* 1 + (1-blendL) .* comboList(ci).color;
        for p = 1:nPlat_c7
            hBar_c7(p).FaceColor = shades_c7(p,:);
            hBar_c7(p).FaceAlpha = 0.8;
        end
        drawnow;
        for p = 1:nPlat_c7
            errorbar(hBar_c7(p).XEndPoints, meanByDay(:,p), semByDay(:,p), ...
                'k', 'LineStyle','none', 'LineWidth', 1.2);
        end

        set(gca, 'XTick', 1:nDays_c7, ...
            'XTickLabel', arrayfun(@num2str, uniqueDays_c7, 'UniformOutput', false));
        title(sprintf('%s  —  %s  (n=%d)', stratTitle_c7, comboList(ci).label, nRats_c));
        xlabel('Day'); ylabel('Mean Strategy Usage');
        pubify_figure_axis_robust(14,14);
        hold off;

        saveas(fH_c7, fullfile(fig_dir, 'Platform', ...
            sprintf('3Platform_%s_%s.png', stratTitle_c7, comboList(ci).label)));
        close(fH_c7);
    end
end

%% 8) Group Strategy x Platform Group
% -----------------------------------------------------------------------
% NEW: Set groupBy to choose colouring/grouping factor for this section.
% -----------------------------------------------------------------------
groupBy_S8 = 'AgeGroup';   % <-- CHANGE THIS PER ANALYSIS

levels_S8 = unique(data1.(groupBy_S8));
nLvls_S8  = numel(levels_S8);

platformShades8 = cell(1, nLvls_S8);
for li = 1:nLvls_S8
    if li <= size(baseCols,1), baseC = baseCols(li,:);
    else, c=hsv(nLvls_S8); baseC=c(li,:); end
    blendLvls = linspace(0.6,0,4)';
    platformShades8{li} = blendLvls .* 1 + (1-blendLvls) .* baseC;
end
sigColors_lvl8 = cell(1,nLvls_S8);
for li=1:nLvls_S8
    if li<=size(baseCols,1), sigColors_lvl8{li}=baseCols(li,:);
    else, c=hsv(nLvls_S8); sigColors_lvl8{li}=c(li,:); end
end
clrBothSig8 = [0.3 0.3 0.3];

strategyGroups = {
    {'thigmotaxis', 'circling', 'randomPath'}, 'Non-Goal Oriented';
    {'scanning', 'chaining'}, 'Procedural';
    {'directedSearch', 'correctedPath', 'directPath','perseverance'}, 'Allocentric'
    };
groupNames8 = strategyGroups(:, 2);
nGroups8 = numel(groupNames8);

groupStrat8 = cell(1, nGroups8);
uniqueDays=unique(data1.Day);
uniqueRats=unique(data1.x_TargetID);
uniquePlatforms=unique(data1.pd);
uniquePlatforms = uniquePlatforms(~isnan(uniquePlatforms));

for g = 1:nGroups8
    currentStrategies = strategyGroups{g,1};
    groupProb = data1.(currentStrategies{1});
    for ii = 2:numel(currentStrategies)
        groupProb = groupProb + data1.(currentStrategies{ii});
    end
    groupStrat8{g} = groupProb;
    groupTitle = groupNames8{g};

    nPlat = numel(uniquePlatforms);
    master_data = {};

    for d = 1:numel(uniqueDays)
        dVal = uniqueDays(d);
        dayData = data1(data1.Day == dVal, :);
        groupProbDay = groupStrat8{g}(data1.Day == dVal);

        for r = 1:numel(uniqueRats)
            ratID = uniqueRats{r};
            isCurrentRat = strcmp(dayData.x_TargetID, ratID);
            if ~any(isCurrentRat), continue; end

            ratLvl = dayData.(groupBy_S8){find(isCurrentRat,1)};
            if iscell(ratLvl), ratLvl = ratLvl{1}; end

            meanUse = nan(1, nPlat);
            for p = 1:nPlat
                pVal = uniquePlatforms(p);
                idx = isCurrentRat & (dayData.pd == pVal);
                if any(idx)
                    meanUse(p) = mean(groupProbDay(idx));
                end
            end

            row = {ratID, ratLvl, dVal};
            for p = 1:nPlat, row = [row, {meanUse(p)}]; end
            master_data(end+1,:) = row;
        end
    end

    measureVars = arrayfun(@(x) sprintf('P%d', x), uniquePlatforms, 'UniformOutput', false);
    measureVars = measureVars(:)';
    varNames = [{'RatID','Age','Day'}, measureVars];
    master_table = cell2table(master_data, 'VariableNames', varNames);

    all_anova = table(); all_tukey = table();
    for d = 1:numel(uniqueDays)
        dVal = uniqueDays(d);
        day_table = master_table(master_table.Day == dVal, :);
        anova_tbl = day_table(:, [{'Age'}, measureVars]);
        anovaResults   = runMixedANOVA(anova_tbl, measureVars, 'Platform');
        postHocResults = runTukeyPostHocMixed(anova_tbl, measureVars, 'Platform');
        anovaResults.Day  = repmat(dVal, height(anovaResults), 1);
        disp(anovaResults)
        postHocResults.Day = repmat(dVal, height(postHocResults), 1);
        anovaResults.Properties.RowNames   = {};
        postHocResults.Properties.RowNames = {};
        all_anova = [all_anova; anovaResults];
        all_tukey = [all_tukey; postHocResults];
    end

    writetable(all_anova, fullfile(processed_dir, sprintf('3%s_PlatformGroups_ANOVA.csv', groupTitle)));
    writetable(all_tukey, fullfile(processed_dir, sprintf('3%s_PlatformGroups_Tukey.csv', groupTitle)));

    % Summary stats
    nDays8 = numel(uniqueDays);
    meanGrp  = cell(1, nLvls_S8);
    semGrp   = cell(1, nLvls_S8);
    dataGrp  = cell(nLvls_S8, nDays8, nPlat);

    for li = 1:nLvls_S8
        meanGrp{li} = nan(nDays8, nPlat);
        semGrp{li}  = nan(nDays8, nPlat);
    end

    for d = 1:nDays8
        dVal = uniqueDays(d);
        day_tbl = master_table(master_table.Day == dVal, :);
        for p = 1:nPlat
            colName = measureVars{p};
            for li = 1:nLvls_S8
                lvl  = levels_S8(li);
                vals = day_tbl{strcmp(day_tbl.Age, lvl), colName};
                meanGrp{li}(d,p) = mean(vals,'omitnan');
                semGrp{li}(d,p)  = std(vals,'omitnan')/sqrt(sum(~isnan(vals)));
                dataGrp{li,d,p}  = vals;
            end
        end
    end

    nCols = nLvls_S8 * nPlat;
    meanMatrix = nan(nDays8, nCols);
    semMatrix  = nan(nDays8, nCols);
    for p = 1:nPlat
        for li = 1:nLvls_S8
            col = (p-1)*nLvls_S8 + li;
            meanMatrix(:, col) = meanGrp{li}(:, p);
            semMatrix(:, col)  = semGrp{li}(:, p);
        end
    end

    fH = figure('Name', groupTitle, 'Color', [1 1 1]);
    sgtitle(groupTitle, 'FontSize', 14, 'FontWeight', 'bold');
    ylim([0 1.5])
    hold on;
    hBar = bar(meanMatrix, 'grouped');
    set(gca, 'XTick', 1:nDays8, 'XTickLabel', arrayfun(@num2str, uniqueDays, 'UniformOutput', false));
    xlabel('Day'); ylabel('Mean Aggregated Strategy Usage');

    for p = 1:nPlat
        for li = 1:nLvls_S8
            col = (p-1)*nLvls_S8 + li;
            shadeIdx = min(p, size(platformShades8{li},1));
            hBar(col).FaceColor = platformShades8{li}(shadeIdx,:);
            hBar(col).FaceAlpha = 0.75;
        end
    end

    drawnow;
    for c = 1:nCols
        xVals = hBar(c).XEndPoints;
        yVals = meanMatrix(:, c);
        eVals = semMatrix(:, c);
        errorbar(xVals, yVals, eVals, 'k', 'LineStyle', 'none', 'LineWidth', 1.2);
    end

    jitterAmount = 0.05;
    for d = 1:nDays8
        for p = 1:nPlat
            for li = 1:nLvls_S8
                col = (p-1)*nLvls_S8 + li;
                pts = dataGrp{li,d,p};
                if isempty(pts), continue; end
                xCenter = hBar(col).XEndPoints(d);
                xJitter = xCenter + (rand(size(pts))-0.5)*jitterAmount;
                scatter(xJitter, pts, 12, 'MarkerFaceColor', hBar(col).FaceColor, ...
                    'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 1);
            end
        end
    end

    sigPairs_all = {}; sigPvals_all = []; sigColors_all = {};
    for iRow = 1:height(all_tukey)
        day_val = all_tukey.Day(iRow);
        dayIdx  = find(uniqueDays == day_val);
        if isempty(dayIdx), continue; end
        p1   = all_tukey.Platform_1(iRow);
        p2   = all_tukey.Platform_2(iRow);
        pVal = all_tukey.pValue(iRow);
        if pVal >= 0.05, continue; end
        compAge = all_tukey.Age{iRow};

        lvlIdx_sig = find(strcmp(levels_S8, compAge), 1);
        if ~isempty(lvlIdx_sig)
            col1 = (p1-1)*nLvls_S8 + lvlIdx_sig;
            col2 = (p2-1)*nLvls_S8 + lvlIdx_sig;
            x1 = hBar(col1).XEndPoints(dayIdx);
            x2 = hBar(col2).XEndPoints(dayIdx);
            sigColor = sigColors_lvl8{lvlIdx_sig};
        else
            x_candidates = [];
            for li = 1:nLvls_S8
                x_candidates(end+1) = hBar((p1-1)*nLvls_S8+li).XEndPoints(dayIdx);
                x_candidates(end+1) = hBar((p2-1)*nLvls_S8+li).XEndPoints(dayIdx);
            end
            x1 = min(x_candidates); x2 = max(x_candidates);
            sigColor = clrBothSig8;
        end
        sigPairs_all{end+1} = [x1, x2];
        sigPvals_all(end+1) = pVal;
        sigColors_all{end+1} = sigColor;
    end
    if ~isempty(sigPairs_all)
        hS = sigstar(sigPairs_all, sigPvals_all);
        for k = 1:size(hS,1)
            set(hS(k), 'Color', sigColors_all{k});
        end
    end

    hold off;
    pubify_figure_axis_robust(14,14);
    saveas(fH, fullfile(fig_dir, 'Platform',sprintf('3Platform_%s.png', groupTitle)));
    close(fH);
end

% ---- Per Sex×Age×APP combo figures for Section 8 ---- %
strategyGroups_c8 = {
    {'thigmotaxis', 'circling', 'randomPath'}, 'Non-Goal Oriented';
    {'scanning', 'chaining'}, 'Procedural';
    {'directedSearch', 'correctedPath', 'directPath','perseverance'}, 'Allocentric'
    };
uniqueDays_c8     = unique(data1.Day);
uniquePlatforms_c8 = unique(data1.pd);
uniquePlatforms_c8 = uniquePlatforms_c8(~isnan(uniquePlatforms_c8));
nPlat_c8 = numel(uniquePlatforms_c8);
nDays_c8 = numel(uniqueDays_c8);

for g = 1:size(strategyGroups_c8,1)
    currentStrats_c8 = strategyGroups_c8{g,1};
    groupProb_c8 = data1.(currentStrats_c8{1});
    for ii = 2:numel(currentStrats_c8)
        groupProb_c8 = groupProb_c8 + data1.(currentStrats_c8{ii});
    end
    gName_c8 = strategyGroups_c8{g,2};

    for ci = 1:nCombos
        idx     = comboList(ci).mask;
        nRats_c = comboList(ci).n;
        rats_c  = unique(data1.x_TargetID(idx));

        meanByDay8 = nan(nDays_c8, nPlat_c8);
        semByDay8  = nan(nDays_c8, nPlat_c8);

        for d = 1:nDays_c8
            dVal  = uniqueDays_c8(d);
            dMask = idx & (data1.Day == dVal);
            for p = 1:nPlat_c8
                pVal  = uniquePlatforms_c8(p);
                pMask = dMask & (data1.pd == pVal);
                ratVals = nan(numel(rats_c),1);
                for r = 1:numel(rats_c)
                    rIdx = strcmp(data1.x_TargetID, rats_c{r}) & pMask;
                    if any(rIdx)
                        ratVals(r) = mean(groupProb_c8(rIdx));
                    end
                end
                meanByDay8(d,p) = mean(ratVals,'omitnan');
                if nRats_c > 1
                    semByDay8(d,p) = std(ratVals,'omitnan') ./ sqrt(sum(~isnan(ratVals)));
                else
                    semByDay8(d,p) = 0;
                end
            end
        end

        fH_c8 = figure; hold on;
        hBar_c8 = bar(meanByDay8, 'grouped');
        blendL8 = linspace(0.6, 0, nPlat_c8)';
        shades_c8 = blendL8 .* 1 + (1-blendL8) .* comboList(ci).color;
        for p = 1:nPlat_c8
            hBar_c8(p).FaceColor = shades_c8(p,:);
            hBar_c8(p).FaceAlpha = 0.8;
        end
        drawnow;
        for p = 1:nPlat_c8
            errorbar(hBar_c8(p).XEndPoints, meanByDay8(:,p), semByDay8(:,p), ...
                'k', 'LineStyle','none', 'LineWidth', 1.2);
        end

        set(gca, 'XTick', 1:nDays_c8, ...
            'XTickLabel', arrayfun(@num2str, uniqueDays_c8, 'UniformOutput', false));
        title(sprintf('%s  —  %s  (n=%d)', gName_c8, comboList(ci).label, nRats_c));
        xlabel('Day'); ylabel('Mean Aggregated Strategy Usage');
        ylim([0 1.5]);
        pubify_figure_axis_robust(14,14);
        hold off;

        saveas(fH_c8, fullfile(fig_dir, 'Platform', ...
            sprintf('3Platform_%s_%s.png', gName_c8, comboList(ci).label)));
        close(fH_c8);
    end
end

%% 9) Strategy Group Diff - per chosen grouping factor separately
% -----------------------------------------------------------------------
% NEW: Set groupBy to choose which factor defines the subplot panels.
%   Each unique level of groupBy_S9 gets its own figure.
% -----------------------------------------------------------------------
groupBy_S9 = 'AgeGroup';   % <-- CHANGE THIS PER ANALYSIS

levels_S9 = unique(data1.(groupBy_S9));
nLvls_S9  = numel(levels_S9);

strategyGroups = {
    {'thigmotaxis','circling','randomPath'}, 'NonGoal';
    {'scanning','chaining'},                 'Procedural';
    {'directedSearch','correctedPath','directPath','perseverance'}, 'Allocentric'
    };
groupLabels = strategyGroups(:,2);
groupNames9  = {'Non‑Goal Oriented','Procedural','Allocentric'};
nGroups9     = numel(groupLabels);
uniqueDays   = unique(data1.Day);

stratColors9 = [
    0.3961, 0.2627, 0.1294;
    1.0000, 0.7020, 0.4000;
    0,      0.5020, 0
    ];

rats  = unique(data1.x_TargetID);
nRats = numel(rats);
nDays = numel(uniqueDays);

ratMeans = nan(nRats, nGroups9 * nDays);
LvlCell  = cell(nRats,1);   % stores groupBy_S9 label per rat

for i = 1:nRats
    ridx        = strcmp(data1.x_TargetID, rats{i});
    LvlCell{i}  = data1.(groupBy_S9){find(ridx,1)};
    col = 1;
    for g = 1:nGroups9
        gp = zeros(sum(ridx),1);
        for s = strategyGroups{g,1}
            gp = gp + data1.(s{1})(ridx);
        end
        for d = 1:nDays
            ratMeans(i, col) = mean(gp(data1.Day(ridx)==uniqueDays(d)));
            col = col + 1;
        end
    end
end

varNames9 = {'RatID','Age'};
for g = 1:nGroups9
    for d = 1:nDays
        varNames9{end+1} = sprintf('%s_Day%d', groupLabels{g}, uniqueDays(d));
    end
end
ratTable9 = table(rats, LvlCell, ...
    ratMeans(:,1),  ratMeans(:,2),  ratMeans(:,3),  ratMeans(:,4), ...
    ratMeans(:,5),  ratMeans(:,6),  ratMeans(:,7),  ratMeans(:,8), ...
    ratMeans(:,9),  ratMeans(:,10), ratMeans(:,11), ratMeans(:,12), ...
    'VariableNames', varNames9);

DaysVec  = repmat(uniqueDays, nGroups9, 1);
StrVec   = repelem(groupLabels(:), nDays, 1);
withinDesign9 = table(DaysVec, categorical(StrVec), ...
    'VariableNames', {'Day','StrategyGroup'});
measureVars9 = ratTable9.Properties.VariableNames(3:end);

% Loop over each level of the chosen grouping factor (was 'young'/'old')
for li = 1:nLvls_S9
    ag = char(levels_S9(li));
    idxAge = strcmp(ratTable9.Age, ag);
    subTbl = ratTable9(idxAge, :);

    formula = sprintf('%s-%s ~ 1', measureVars9{1}, measureVars9{end});
    rm = fitrm(subTbl, formula, 'WithinDesign', withinDesign9);

    anovaResults9 = ranova(rm, 'WithinModel', 'Day*StrategyGroup');
    writetable(anovaResults9, fullfile(processed_dir, sprintf('ANOVA_%s.csv', ag)));

    postHocStrat = multcompare(rm, 'StrategyGroup', 'By', 'Day', 'ComparisonType', 'tukey-kramer');
    writetable(postHocStrat, fullfile(processed_dir, sprintf('PostHoc_Strat_%s.csv', ag)));
    postHocDay   = multcompare(rm, 'Day', 'By', 'StrategyGroup', 'ComparisonType', 'tukey-kramer');
    writetable(postHocDay,   fullfile(processed_dir, sprintf('PostHoc_Day_%s.csv', ag)));

    meanData = zeros(nDays, nGroups9);
    semData  = zeros(nDays, nGroups9);
    dataMeansPerGroup = cell(nGroups9, nDays);
    for g = 1:nGroups9
        for d = 1:nDays
            colName = sprintf('%s_Day%d', groupLabels{g}, uniqueDays(d));
            vals = subTbl.(colName);
            meanData(d,g) = mean(vals);
            semData(d,g)  = std(vals)/sqrt(numel(vals));
            dataMeansPerGroup{g,d} = vals;
        end
    end

    f = figure('Position',[95,100,1100,630]);
    hold on;
    hBar = bar(uniqueDays, meanData, 'grouped', 'BarWidth', 0.8);
    for g = 1:nGroups9
        hBar(g).FaceColor = stratColors9(g,:);
        hBar(g).FaceAlpha = 0.35;
    end
    drawnow;
    barCenters = nan(nDays, nGroups9);
    for g = 1:nGroups9, barCenters(:,g) = hBar(g).XEndPoints'; end

    for g = 1:nGroups9
        errorbar(barCenters(:,g), meanData(:,g), semData(:,g), 'k', ...
            'LineStyle','none','LineWidth',1.5);
    end

    jit = 0.08;
    for g = 1:nGroups9
        for d = 1:nDays
            pts = dataMeansPerGroup{g,d};
            xj  = barCenters(d,g) + (rand(size(pts))-0.5)*jit;
            scatter(xj, pts, 20, 'MarkerFaceColor', stratColors9(g,:), ...
                'MarkerEdgeColor','none','MarkerFaceAlpha',0.8);
        end
    end

    % Between-strategy sigstar
    sigPairs = {}; sigP = []; sigCols = {};
    cs = multcompare(rm, 'StrategyGroup', 'By', 'Day', 'ComparisonType', 'tukey-kramer');
    for rr = 1:height(cs)
        pval = cs.pValue(rr);
        if pval < 0.05
            dVal = cs.Day(rr);
            dIdx = find(uniqueDays == dVal, 1);
            lvl1 = char(cs.StrategyGroup_1(rr));
            lvl2 = char(cs.StrategyGroup_2(rr));
            g1 = find(strcmp(groupLabels, lvl1), 1);
            g2 = find(strcmp(groupLabels, lvl2), 1);
            if ~isempty(dIdx) && ~isempty(g1) && ~isempty(g2) && (g1 < g2)
                sigPairs{end+1} = [barCenters(dIdx, g1), barCenters(dIdx, g2)];
                sigP(end+1)     = pval;
                sigCols{end+1}  = [0.5,0.5,0.5];
            end
        end
    end
    if ~isempty(sigPairs)
        hSig = sigstar(sigPairs, sigP);
        for k = 1:size(hSig,1)
            set(hSig(k,1), 'Color', sigCols{k});
            set(hSig(k,2), 'Color', sigCols{k});
        end
    end

    % Within-strategy day sigstar (unchanged from original)
    barCenters = nan(nDays, nGroups9);
    for g = 1:nGroups9, barCenters(:,g) = hBar(g).XEndPoints'; end
    sigPairs = {}; sigP = []; sigCols = {};
    darkGray = [0.2 0.2 0.2];
    cd9 = multcompare(rm, 'Day','By','StrategyGroup','ComparisonType','tukey-kramer');
    for d = 1:nDays-1
        sel = cd9.Day_1==uniqueDays(d) & cd9.Day_2==uniqueDays(d+1);
        sub = cd9(sel,:);
        pv = nan(1,nGroups9);
        for rr = 1:height(sub)
            gidx = find(strcmp(groupLabels, char(sub.StrategyGroup(rr))),1);
            pv(gidx) = sub.pValue(rr);
        end
        sigIdx = find(pv<0.05);
        switch numel(sigIdx)
            case 0
            case 3
                sigPairs{end+1} = [mean(barCenters(d,:)), mean(barCenters(d+1,:))];
                sigP(end+1) = min(pv(sigIdx)); sigCols{end+1} = [0 0 0];
            case 2
                sigPairs{end+1} = [mean(barCenters(d,sigIdx)), mean(barCenters(d+1,sigIdx))];
                sigP(end+1) = min(pv(sigIdx)); sigCols{end+1} = darkGray;
            case 1
                g = sigIdx;
                sigPairs{end+1} = [barCenters(d,g), barCenters(d+1,g)];
                sigP(end+1) = pv(g); sigCols{end+1} = stratColors9(g,:);
        end
    end
    valid = cellfun(@(c)isnumeric(c)&&numel(c)==2, sigPairs);
    sigPairs=sigPairs(valid); sigP=sigP(valid); sigCols=sigCols(valid);
    if ~isempty(sigPairs)
        hSig = sigstar(sigPairs, sigP);
        for k = 1:size(hSig,1)
            set(hSig(k,1), 'Color', sigCols{k});
            set(hSig(k,2), 'Color', sigCols{k});
        end
    end

    title(sprintf('%s', upper(ag)));
    xlabel('Day'); ylabel('Probability for Strategy Group');
    xticks(uniqueDays);
    legend(groupNames9, 'Location','northeastoutside');
    ylim([0 1.2]);
    pubify_figure_axis_robust(16,16);
    hold off;

    saveas(f, fullfile(fig_dir, sprintf('DayStrategyUse_%s.png', ag)));
    close(f);
end

% ---- Per Sex×Age×APP combo figures for Section 9 ---- %
% One figure per combo: all three strategy groups as grouped bars across days.
% No stats — descriptive only, since each combo may have very few rats.
strategyGroups_c9 = {
    {'thigmotaxis','circling','randomPath'}, 'NonGoal',    'Non-Goal Oriented';
    {'scanning','chaining'},                 'Procedural',  'Procedural';
    {'directedSearch','correctedPath','directPath','perseverance'}, 'Allocentric', 'Allocentric'
    };
stratColors_c9 = [
    0.3961, 0.2627, 0.1294;
    1.0000, 0.7020, 0.4000;
    0,      0.5020, 0
    ];
nGroups_c9   = size(strategyGroups_c9, 1);
uniqueDays_c9 = unique(data1.Day);
nDays_c9      = numel(uniqueDays_c9);

for ci = 1:nCombos
    idx     = comboList(ci).mask;
    nRats_c = comboList(ci).n;
    rats_c  = unique(data1.x_TargetID(idx));

    meanData_c9 = nan(nDays_c9, nGroups_c9);
    semData_c9  = nan(nDays_c9, nGroups_c9);
    ratMeans_c9 = cell(nGroups_c9, nDays_c9);  % per-rat means for scatter

    for g = 1:nGroups_c9
        currentStrats_c9 = strategyGroups_c9{g,1};
        groupProb_c9 = data1.(currentStrats_c9{1});
        for ii = 2:numel(currentStrats_c9)
            groupProb_c9 = groupProb_c9 + data1.(currentStrats_c9{ii});
        end

        for d = 1:nDays_c9
            dVal  = uniqueDays_c9(d);
            ratV  = nan(numel(rats_c),1);
            for r = 1:numel(rats_c)
                rIdx = strcmp(data1.x_TargetID, rats_c{r}) & idx & (data1.Day == dVal);
                if any(rIdx)
                    ratV(r) = mean(groupProb_c9(rIdx));
                end
            end
            ratMeans_c9{g,d}   = ratV(~isnan(ratV));
            meanData_c9(d,g)   = mean(ratV,'omitnan');
            if nRats_c > 1
                semData_c9(d,g) = std(ratV,'omitnan') ./ sqrt(sum(~isnan(ratV)));
            else
                semData_c9(d,g) = 0;
            end
        end
    end

    f_c9 = figure('Position',[95,100,1100,630]); hold on;
    hBar_c9 = bar(uniqueDays_c9, meanData_c9, 'grouped', 'BarWidth', 0.8);
    for g = 1:nGroups_c9
        hBar_c9(g).FaceColor = stratColors_c9(g,:);
        hBar_c9(g).FaceAlpha = 0.35;
    end
    drawnow;
    barCenters_c9 = nan(nDays_c9, nGroups_c9);
    for g = 1:nGroups_c9
        barCenters_c9(:,g) = hBar_c9(g).XEndPoints';
    end
    % Error bars
    for g = 1:nGroups_c9
        errorbar(barCenters_c9(:,g), meanData_c9(:,g), semData_c9(:,g), ...
            'k', 'LineStyle','none', 'LineWidth', 1.5);
    end
    % Jittered per-rat scatter
    jit_c9 = 0.08;
    for g = 1:nGroups_c9
        for d = 1:nDays_c9
            pts = ratMeans_c9{g,d};
            if isempty(pts), continue; end
            xj = barCenters_c9(d,g) + (rand(size(pts))-0.5)*jit_c9;
            scatter(xj, pts, 20, 'MarkerFaceColor', stratColors_c9(g,:), ...
                'MarkerEdgeColor','none', 'MarkerFaceAlpha', 0.8);
        end
    end

    title(sprintf('%s  (n=%d)', comboList(ci).label, nRats_c));
    xlabel('Day'); ylabel('Probability for Strategy Group');
    xticks(uniqueDays_c9);
    legend({strategyGroups_c9{:,3}}, 'Location','northeastoutside');
    ylim([0 1.2]);
    pubify_figure_axis_robust(16,16);
    hold off;

    saveas(f_c9, fullfile(fig_dir, ...
        sprintf('DayStrategyUse_%s.png', comboList(ci).label)));
    close(f_c9);
end

%% 10) Individual strategy differences per day per grouping level
% -----------------------------------------------------------------------
% NEW: Set groupBy for this section.
% -----------------------------------------------------------------------
groupBy_S10 = 'AgeGroup';   % <-- CHANGE THIS PER ANALYSIS

levels_S10 = unique(data1.(groupBy_S10));
nLvls_S10  = numel(levels_S10);

strategyList = { ...
    'thigmotaxis', 'circling',      'randomPath', ...
    'scanning',    'chaining', ...
    'directedSearch','correctedPath','directPath' };

stratColors10 = [
    0.850, 0.325, 0.098;
    0.094, 0.509, 0.800;
    0.467, 0.675, 0.188;
    0.800, 0.475, 0.655;
    0.894, 0.766, 0.039;
    0.436, 0.600, 0.800;
    0.905, 0.161, 0.541;
    0.400, 0.400, 0.400
    ];

nStrats10   = numel(strategyList);
uniqueDays  = unique(data1.Day);
nDays       = numel(uniqueDays);
rats        = unique(data1.x_TargetID);
nRats       = numel(rats);

ratMeans10 = nan(nRats, nStrats10 * nDays);
LvlCell10  = cell(nRats,1);

for i = 1:nRats
    ridx          = strcmp(data1.x_TargetID, rats{i});
    LvlCell10{i}  = data1.(groupBy_S10){find(ridx,1)};
    col = 1;
    for sIdx = 1:nStrats10
        stratVec = data1.(strategyList{sIdx})(ridx);
        for d = 1:nDays
            dayMask = (data1.Day(ridx)==uniqueDays(d));
            ratMeans10(i,col) = mean(stratVec(dayMask));
            col = col + 1;
        end
    end
end

varNames10 = {'RatID','Age'};
for sIdx = 1:nStrats10
    for d = 1:nDays
        varNames10{end+1} = sprintf('%s_Day%d', strategyList{sIdx}, uniqueDays(d));
    end
end

ratTable10 = cell2table([rats, LvlCell10, num2cell(ratMeans10)], ...
    'VariableNames', varNames10);

DaysVec10  = repmat(uniqueDays, nStrats10, 1);
StrVec10   = repelem(strategyList(:), nDays, 1);
withinDesign10 = table(categorical(DaysVec10), categorical(StrVec10), ...
    'VariableNames', {'Day','Strategy'});
measureVars10 = ratTable10.Properties.VariableNames(3:end);

for li = 1:nLvls_S10
    ag = char(levels_S10(li));
    idxAge = strcmp(ratTable10.Age, ag);
    subTbl = ratTable10(idxAge, :);

    formula10 = sprintf('%s-%s ~ 1', measureVars10{1}, measureVars10{end});
    rm = fitrm(subTbl, formula10, 'WithinDesign', withinDesign10, ...
        'WithinModel', 'Day*Strategy');

    anovaResults10 = ranova(rm, 'WithinModel', 'Day*Strategy');
    writetable(anovaResults10, fullfile(processed_dir, sprintf('All_ANOVA_%s.csv', ag)));
    disp(anovaResults10)

    postHocStrat10 = multcompare(rm, 'Strategy', 'By', 'Day', 'ComparisonType', 'tukey-kramer');
    writetable(postHocStrat10, fullfile(processed_dir, sprintf('ALL_PostHoc_Strat_%s.csv', ag)));
    postHocDay10 = multcompare(rm, 'Day', 'By', 'Strategy', 'ComparisonType', 'tukey-kramer');
    writetable(postHocDay10, fullfile(processed_dir, sprintf('ALL_PostHoc_Day_%s.csv', ag)));

    meanData = zeros(nDays, nStrats10);
    semData  = zeros(nDays, nStrats10);
    dataMeansPerGroup = cell(nStrats10, nDays);
    for g = 1:nStrats10
        for d = 1:nDays
            colName = sprintf('%s_Day%d', strategyList{g}, uniqueDays(d));
            vals = subTbl.(colName);
            meanData(d,g) = mean(vals);
            semData(d,g)  = std(vals)/sqrt(numel(vals));
            dataMeansPerGroup{g,d} = vals;
        end
    end

    f = figure('Position',[95,100,1100,630]);
    hold on;
    hBar = bar(uniqueDays, meanData, 'grouped', 'BarWidth', 0.8);
    for g = 1:nStrats10
        hBar(g).FaceColor = stratColors10(g,:);
        hBar(g).FaceAlpha = 0.35;
    end
    drawnow;
    barCenters10 = nan(nDays, nStrats10);
    for g = 1:nStrats10, barCenters10(:,g) = hBar(g).XEndPoints'; end

    jit = 0.03;
    for g = 1:nStrats10
        for d = 1:nDays
            pts = dataMeansPerGroup{g,d};
            xj  = barCenters10(d,g) + (rand(size(pts))-0.5)*jit;
            scatter(xj, pts, 20, 'MarkerFaceColor', stratColors10(g,:), ...
                'MarkerEdgeColor','none','MarkerFaceAlpha',0.8);
        end
    end
    for g = 1:nStrats10
        errorbar(barCenters10(:,g), meanData(:,g), semData(:,g), 'k', ...
            'LineStyle','none','LineWidth',2);
    end

    title(sprintf('%s', upper(ag)));
    xlabel('Day'); ylabel('Probability for Strategy Group');
    xticks(uniqueDays);
    ylim([0 0.8]);
    legend(strategyList, 'Location','northeast');
    pubify_figure_axis_robust(16,16);
    hold off;

    exportgraphics(f, fullfile(fig_dir, sprintf('Day_ALLStrategies_%s.png', ag)), ...
        'Resolution', 450);
    close(f);
end

% ---- Per Sex×Age×APP combo figures for Section 10 ---- %
% One figure per combo: all individual strategies as grouped bars across days.
strategyList_c10 = { ...
    'thigmotaxis', 'circling',       'randomPath', ...
    'scanning',    'chaining', ...
    'directedSearch','correctedPath','directPath' };
stratColors_c10 = [
    0.850, 0.325, 0.098;
    0.094, 0.509, 0.800;
    0.467, 0.675, 0.188;
    0.800, 0.475, 0.655;
    0.894, 0.766, 0.039;
    0.436, 0.600, 0.800;
    0.905, 0.161, 0.541;
    0.400, 0.400, 0.400
    ];
nStrats_c10  = numel(strategyList_c10);
uniqueDays_c10 = unique(data1.Day);
nDays_c10      = numel(uniqueDays_c10);

for ci = 1:nCombos
    idx     = comboList(ci).mask;
    nRats_c = comboList(ci).n;
    rats_c  = unique(data1.x_TargetID(idx));

    meanData_c10 = nan(nDays_c10, nStrats_c10);
    semData_c10  = nan(nDays_c10, nStrats_c10);
    ratMeans_c10 = cell(nStrats_c10, nDays_c10);

    for s = 1:nStrats_c10
        stratVec_c10 = data1.(strategyList_c10{s});
        for d = 1:nDays_c10
            dVal = uniqueDays_c10(d);
            ratV = nan(numel(rats_c),1);
            for r = 1:numel(rats_c)
                rIdx = strcmp(data1.x_TargetID, rats_c{r}) & idx & (data1.Day == dVal);
                if any(rIdx)
                    ratV(r) = mean(stratVec_c10(rIdx));
                end
            end
            ratMeans_c10{s,d}   = ratV(~isnan(ratV));
            meanData_c10(d,s)   = mean(ratV,'omitnan');
            if nRats_c > 1
                semData_c10(d,s) = std(ratV,'omitnan') ./ sqrt(sum(~isnan(ratV)));
            else
                semData_c10(d,s) = 0;
            end
        end
    end

    f_c10 = figure('Position',[95,100,1100,630]); hold on;
    hBar_c10 = bar(uniqueDays_c10, meanData_c10, 'grouped', 'BarWidth', 0.8);
    for s = 1:nStrats_c10
        hBar_c10(s).FaceColor = stratColors_c10(s,:);
        hBar_c10(s).FaceAlpha = 0.35;
    end
    drawnow;
    barCenters_c10 = nan(nDays_c10, nStrats_c10);
    for s = 1:nStrats_c10
        barCenters_c10(:,s) = hBar_c10(s).XEndPoints';
    end
    % Jittered per-rat scatter
    jit_c10 = 0.03;
    for s = 1:nStrats_c10
        for d = 1:nDays_c10
            pts = ratMeans_c10{s,d};
            if isempty(pts), continue; end
            xj = barCenters_c10(d,s) + (rand(size(pts))-0.5)*jit_c10;
            scatter(xj, pts, 20, 'MarkerFaceColor', stratColors_c10(s,:), ...
                'MarkerEdgeColor','none', 'MarkerFaceAlpha', 0.8);
        end
    end
    % Error bars
    for s = 1:nStrats_c10
        errorbar(barCenters_c10(:,s), meanData_c10(:,s), semData_c10(:,s), ...
            'k', 'LineStyle','none', 'LineWidth', 2);
    end

    title(sprintf('%s  (n=%d)', comboList(ci).label, nRats_c));
    xlabel('Day'); ylabel('Probability of Strategy Use');
    xticks(uniqueDays_c10);
    ylim([0 0.8]);
    legend(strategyList_c10, 'Location','northeast');
    pubify_figure_axis_robust(16,16);
    hold off;

    exportgraphics(f_c10, fullfile(fig_dir, ...
        sprintf('Day_ALLStrategies_%s.png', comboList(ci).label)), ...
        'Resolution', 450);
    close(f_c10);
end

%% 11) Within strategy group ANOVA and comparisons - STATS
% -----------------------------------------------------------------------
% No groupBy needed here — this section loops within each strategy subgroup
% and already uses Age as the between-subjects factor via fitrm.
% If you want to substitute a different between-subjects factor, change
% the formula string in the fitrm call below (replace 'Age' with e.g. 'Sex').
% -----------------------------------------------------------------------
betweenFactor_S11 = 'Age';   % <-- CHANGE THIS PER ANALYSIS
%   'Age'      → uses the Age column of the rat table (= groupBy used in S3)
%   'Sex'      → add a Sex column to the rat table and use that
%   'APP'      → similarly

strategyGroups11 = { ...
    {'thigmotaxis','circling','randomPath'},'NonGoal'; ...
    {'scanning','chaining'},  'Procedural'; ...
    {'directedSearch','correctedPath','directPath'},  'Allocentric'};

groupPrettyNames11 = {'Non-Goal Oriented','Procedural','Allocentric'};

stratColors11 = [
    0.850 0.325 0.098;
    0.094 0.509 0.800;
    0.467 0.675 0.188;
    0.800 0.475 0.655;
    0.894 0.766 0.039;
    0.436 0.600 0.800;
    0.905 0.161 0.541;
    0.400 0.400 0.400];

outDir = fullfile(processed_dir,'Strategy_Day_Comparisons');
if ~exist(outDir,'dir'); mkdir(outDir); end

uniqueDays11 = unique(data1.Day); nDays11 = numel(uniqueDays11);
rats11       = unique(data1.x_TargetID);
nRats11      = numel(rats11);

%--- Descriptive Mean ± SD ---%
descTab = table();
for g = 1:size(strategyGroups11,1)
    stratList11 = strategyGroups11{g,1};
    subgroupTag = strategyGroups11{g,2};
    for dIdx = 1:nDays11
        dVal = uniqueDays11(dIdx);
        for ageTag = levels_S9'   % use levels from groupBy used in S9, or hardcode
            ageTag = char(ageTag);
            dayAgeMask = (data1.Day == dVal) & strcmp(data1.(groupBy_S9), ageTag);
            for sIdx = 1:numel(stratList11)
                stratName11 = stratList11{sIdx};
                vals = data1.(stratName11)(dayAgeMask);
                m  = mean(vals,'omitnan');
                sd = std(vals,'omitnan');
                descTab = [descTab;
                    table(categorical({subgroupTag}), categorical({ageTag}), dVal, ...
                    categorical({stratName11}), m, sd, ...
                    'VariableNames',{'Group','Age','Day','Strategy','Mean','SD'})];
            end
        end
    end
end
writetable(descTab, fullfile(outDir,'MeanSD_StrategyUse.csv'));

%--- Loop through each subgroup for RM-ANOVAs ---%
for g = 1:size(strategyGroups11,1)
    stratList11 = strategyGroups11{g,1};
    subgroupTag = strategyGroups11{g,2};
    prettyName  = groupPrettyNames11{g};
    nStrats11   = numel(stratList11);

    ratMeans11 = nan(nRats11, nStrats11*nDays11);
    AgeCell11  = cell(nRats11,1);

    for r = 1:nRats11
        rMask        = strcmp(data1.x_TargetID, rats11{r});
        AgeCell11{r} = data1.(groupBy_S9){find(rMask,1)};
        col = 1;
        for sIdx = 1:nStrats11
            vec    = data1.(stratList11{sIdx})(rMask);
            dayVec = data1.Day(rMask);
            for dIdx = 1:nDays11
                ratMeans11(r,col) = mean(vec(dayVec==uniqueDays11(dIdx)));
                col = col + 1;
            end
        end
    end

    varNames11 = {'RatID','Age'};
    for sIdx = 1:nStrats11
        for dIdx = 1:nDays11
            varNames11{end+1} = sprintf('%s_Day%d', stratList11{sIdx}, uniqueDays11(dIdx));
        end
    end
    T11 = cell2table([rats11, AgeCell11, num2cell(ratMeans11)], 'VariableNames',varNames11);

    DaysVec11 = repmat(uniqueDays11, nStrats11,1);
    StrVec11  = repelem(stratList11(:), nDays11,1);
    WD11      = table(categorical(DaysVec11), categorical(StrVec11), ...
        'VariableNames',{'Day','Strategy'});
    measVars11 = T11.Properties.VariableNames(3:end);

    rm11 = fitrm(T11, sprintf('%s-%s ~ Age', measVars11{1}, measVars11{end}), ...
        'WithinDesign', WD11, 'WithinModel', 'Day*Strategy');

    a11 = ranova(rm11,'WithinModel','Day*Strategy');
    writetable(a11, fullfile(outDir, sprintf('RM_ANOVA_%s.csv', subgroupTag)), ...
        'WriteRowNames', false);

    phS11 = multcompare(rm11,'Strategy','By','Day','ComparisonType','tukey-kramer');
    phS11.sigFlag = phS11.pValue < 0.05;
    writetable(phS11, fullfile(outDir,sprintf('PostHoc_Strategy_%s.csv',subgroupTag)));

    phD11 = multcompare(rm11,'Day','By','Strategy','ComparisonType','tukey-kramer');
    phD11.sigFlag = phD11.pValue < 0.05;
    writetable(phD11, fullfile(outDir,sprintf('PostHoc_Day_%s.csv',subgroupTag)));

    % Within-level (e.g. within-age) ANOVAs
    for li = 1:nLvls_S9
        ageTag11 = char(levels_S9(li));
        idx11    = strcmp(T11.Age, ageTag11);
        subT11   = T11(idx11,:);
        rmW11    = fitrm(subT11, sprintf('%s-%s~1',measVars11{1},measVars11{end}), ...
            'WithinDesign',WD11,'WithinModel','Day*Strategy');
        aW11     = ranova(rmW11,'WithinModel','Day*Strategy');
        writetable(aW11, fullfile(outDir, ...
            sprintf('WithinAge_ANOVA_%s_%s.csv',subgroupTag,upper(ageTag11))), ...
            'WriteRowNames', true);
    end

    % Verify line+errorbar plots per strategy group
    markers = {'o','s','^','d','v','>','<','p','h'};
    figure('Name',subgroupTag,'Color','w'); hold on;

    for li = 1:nLvls_S9
        ageTag11 = char(levels_S9(li));
        colC = clrMap_S3{min(li,numel(clrMap_S3))};
        ageShift = (li - (nLvls_S9+1)/2)*0.08;

        for sIdx = 1:nStrats11
            stratName11 = stratList11{sIdx};
            marker = markers{mod(sIdx-1,numel(markers))+1};
            meanVals = nan(size(uniqueDays11));
            semVals  = nan(size(uniqueDays11));
            for d = 1:nDays11
                mask = (data1.Day==uniqueDays11(d)) & strcmp(data1.(groupBy_S9), ageTag11);
                v    = data1.(stratName11)(mask);
                meanVals(d) = mean(v,'omitnan');
                semVals(d)  = std(v,'omitnan')/sqrt(sum(mask));
            end
            stratShift = (sIdx - (nStrats11+1)/2)*0.15;
            xPos = uniqueDays11 + ageShift + stratShift;
            errorbar(xPos, meanVals, semVals, '-o', ...
                'Color', colC, 'Marker', marker, ...
                'MarkerFaceColor', colC, 'MarkerSize', 6, 'LineWidth', 1.5, ...
                'DisplayName', sprintf('%s (%s)', stratName11, ageTag11));
        end
    end

    xlabel('Day'); ylabel('Mean Strategy Probability');
    title(sprintf('%s: Strategy × Day by Group', subgroupTag), 'Interpreter','none');
    xticks(uniqueDays11);
    legend('Location','bestoutside');
    hold off;
end

% ---- Per Sex×Age×APP combo figures for Section 11 ---- %
% One figure per strategy subgroup per combo: individual strategies as
% error-bar line plots across days. Descriptive only — no stats per combo.
strategyGroups_c11 = { ...
    {'thigmotaxis','circling','randomPath'}, 'NonGoal',   'Non-Goal Oriented'; ...
    {'scanning','chaining'},                 'Procedural', 'Procedural'; ...
    {'directedSearch','correctedPath','directPath'}, 'Allocentric', 'Allocentric'};

% Strategy colours consistent with Section 11
stratColors_c11 = [
    0.850 0.325 0.098;   % thigmotaxis / non-goal 1
    0.094 0.509 0.800;   % circling    / non-goal 2
    0.467 0.675 0.188;   % randomPath / non-goal 3
    0.800 0.475 0.655;   % scanning    / procedural 1
    0.894 0.766 0.039;   % chaining    / procedural 2
    0.436 0.600 0.800;   % directed    / allocentric 1
    0.905 0.161 0.541;   % corrected   / allocentric 2
    0.400 0.400 0.400];  % directPath / allocentric 3

markers_c11 = {'o','s','^','d','v','>','<','p'};

uniqueDays_c11 = unique(data1.Day);
nDays_c11      = numel(uniqueDays_c11);

% Accumulate colour index across subgroups so each strategy gets a unique colour
globalStratIdx = 0;

for g = 1:size(strategyGroups_c11,1)
    stratList_c11   = strategyGroups_c11{g,1};
    subgroupTag_c11 = strategyGroups_c11{g,2};
    prettyName_c11  = strategyGroups_c11{g,3};
    nStrats_c11     = numel(stratList_c11);

    % Colour indices for this subgroup's strategies
    stratColIdx = globalStratIdx + (1:nStrats_c11);
    globalStratIdx = globalStratIdx + nStrats_c11;

    for ci = 1:nCombos
        idx     = comboList(ci).mask;
        nRats_c = comboList(ci).n;
        rats_c  = unique(data1.x_TargetID(idx));

        f_c11 = figure('Name', sprintf('%s_%s', subgroupTag_c11, comboList(ci).label), ...
            'Color', 'w'); hold on;

        for sIdx = 1:nStrats_c11
            stratName_c11 = stratList_c11{sIdx};
            colIdx = min(stratColIdx(sIdx), size(stratColors_c11,1));
            colC   = stratColors_c11(colIdx,:);
            marker = markers_c11{mod(sIdx-1, numel(markers_c11))+1};

            meanVals_c11 = nan(nDays_c11,1);
            semVals_c11  = nan(nDays_c11,1);

            for d = 1:nDays_c11
                dVal = uniqueDays_c11(d);
                ratV = nan(numel(rats_c),1);
                for r = 1:numel(rats_c)
                    rIdx = strcmp(data1.x_TargetID, rats_c{r}) & idx & (data1.Day == dVal);
                    if any(rIdx)
                        ratV(r) = mean(data1.(stratName_c11)(rIdx));
                    end
                end
                meanVals_c11(d) = mean(ratV,'omitnan');
                if nRats_c > 1
                    semVals_c11(d) = std(ratV,'omitnan') ./ sqrt(sum(~isnan(ratV)));
                else
                    semVals_c11(d) = 0;
                end
            end

            % Slight horizontal jitter per strategy so lines don't overlap
            stratShift_c11 = (sIdx - (nStrats_c11+1)/2) * 0.12;
            xPos_c11 = uniqueDays_c11 + stratShift_c11;

            errorbar(xPos_c11, meanVals_c11, semVals_c11, '-o', ...
                'Color', colC, 'Marker', marker, ...
                'MarkerFaceColor', colC, 'MarkerEdgeColor', 'none', ...
                'MarkerSize', 6, 'LineWidth', 1.5, ...
                'DisplayName', stratName_c11);
        end

        title(sprintf('%s  —  %s  (n=%d)', prettyName_c11, comboList(ci).label, nRats_c), ...
            'Interpreter', 'none');
        xlabel('Day'); ylabel('Mean Strategy Probability');
        xticks(uniqueDays_c11);
        legend('Location','bestoutside');
        pubify_figure_axis_robust(14,14);
        hold off;

        saveas(f_c11, fullfile(fig_dir, ...
            sprintf('StrategyDay_%s_%s.png', subgroupTag_c11, comboList(ci).label)));
        close(f_c11);
    end
end