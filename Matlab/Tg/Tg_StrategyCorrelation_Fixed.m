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

% Create output subfolders if they don't exist yet
subfolders = {'CIPL_Strategy', 'StrategyUse', 'Entropy', 'Platform'};
for sf = subfolders
    p = fullfile(fig_dir, sf{1});
    if ~exist(p, 'dir'), mkdir(p); end
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

% SexGeno: e.g. "M-WT", "F-APP+" — slash stripped so it is safe for
% filenames and can be used directly as a groupBy option throughout the script.
data1.SexGeno = data1.Sex + "-" + strrep(data1.APP, '/', '');

% Build a combined Group label: AgeGroup-Sex-APP  e.g. "Young-M-WT"
% (This is used in sections that want all-combinations colour mapping)
data1.Group = data1.AgeGroup + "-" + data1.Sex + "-" + data1.APP;

% -----------------------------------------------------------------------
% NEW: Derive the full list of groups actually present in the data
%      and assign a colour to each one.
% -----------------------------------------------------------------------
grpList = unique(data1.Group);          % e.g. ["Mid-F-APP/+","Mid-F-WT", ...]
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

%% 1 ) CIPL Across Days - SexGeno Groups
% CIPL across days by SexGeno group

uniqueDays_cipl = unique(data1.x_Day);
sgGroups_cipl   = unique(data1.SexGeno);
nSG_cipl        = numel(sgGroups_cipl);

if nSG_cipl <= size(baseCols,1)
    sgCols_cipl = num2cell(baseCols(1:nSG_cipl,:), 2);
else
    c = hsv(nSG_cipl); sgCols_cipl = num2cell(c, 2);
end

% Build per-rat daily mean CIPL from data2 (which holds Platform_CIPL)
% We match via Animal ID and Day
uniqueRats_cipl = unique(data1.x_TargetID);
nRats_cipl      = numel(uniqueRats_cipl);
ratDayCIPL      = nan(nRats_cipl, numel(uniqueDays_cipl));
ratSexGeno      = strings(nRats_cipl, 1);

for r = 1:nRats_cipl
    ratID = uniqueRats_cipl{r};
    rIdx  = strcmp(data1.x_TargetID, ratID);
    ratSexGeno(r) = data1.SexGeno(find(rIdx,1));
    animalNum = str2double(ratID);

    for d = 1:numel(uniqueDays_cipl)
        dVal  = uniqueDays_cipl(d);
        % Get all trials for this rat on this day from data1
        trialIdx = rIdx & (data1.x_Day == dVal);
        trials   = data1.x_Trial(trialIdx);
        if isempty(trials), continue; end
        % Match to platformScores (already aligned to filtered data1)
        ratDayCIPL(r,d) = mean(platformScores(trialIdx), 'omitnan');
    end
end

f_cipl = figure; hold on;
for sgi = 1:nSG_cipl
    sg     = sgGroups_cipl(sgi);
    mask   = ratSexGeno == sg;
    vals   = ratDayCIPL(mask, :);
    nR     = sum(mask);
    meanV  = mean(vals, 1, 'omitnan');
    if nR > 1
        semV = std(vals, 0, 1, 'omitnan') ./ sqrt(nR);
    else
        semV = zeros(size(meanV));
    end
    errorbar(uniqueDays_cipl, meanV, semV, '-o', ...
        'Color', sgCols_cipl{sgi}, 'LineWidth', 3.5, ...
        'MarkerFaceColor', sgCols_cipl{sgi}, 'MarkerEdgeColor', 'none', ...
        'MarkerSize', 9, ...
        'DisplayName', char(sg));
end

title('Mean CIPL Score Across Days', ...
    'FontSize', 24, 'FontWeight', 'bold');
xlabel('Day', 'FontSize', 18, 'FontWeight', 'bold');
xlim([.9, 4.5]);
ylabel('CIPL Score (m·s)', 'FontSize', 18, 'FontWeight', 'bold');
xticks(uniqueDays_cipl);
legend('Location', 'northeast', 'FontSize', 12);
set(gca, 'FontSize', 18, 'FontWeight', 'bold', ...
    'LineWidth', 3, 'Box', 'off', 'TickDir', 'out');
hold off;

saveas(f_cipl, fullfile(fig_dir, 'CIPL_Strategy', 'CIPL_AcrossDays_SexGeno.png'));
close;



%% 2) CIPL vs Strategy Probabilities
% -----------------------------------------------------------------------
%   Set groupBy to choose which factor to colour by in this section.
%   Options: 'AgeGroup'  → Young / Mid / Old
%            'Sex'       → M / F  (or whatever values are in data1.Sex)
%            'APP'       → WT / APP/+
%            'Group'     → all combinations (AgeGroup_Sex_APP)
%            'SexGeno'   → all combinations (SexAPP)
%            'Age'       → individual numeric age values (one colour each)
% -----------------------------------------------------------------------
groupBy_S1 = 'SexGeno';   % <-- CHANGE THIS PER ANALYSIS

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

%% 3) CIPL vs Group Probabilities

groupBy_S2 = 'SexGeno';   % <-- CHANGE THIS PER ANALYSIS

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

    title(sprintf('%s', groupNames{s}));
    xlabel('CIPL Score (m.s)');
    ylabel('Probability');
    legend('show', 'Location', 'northeastoutside');
    xlim([0 60]);
    ylim([0 1]);
    pubify_figure_axis_robust(21,21)
    hold off;

    saveas(f, fullfile(fig_dir, 'CIPL_Strategy', sprintf('CIPL_%s',groupNames{s})),'png');
end


%% 4) CIPL vs Sex and Genotype
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
        pubify_figure_axis_robust(14,14)
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
        pubify_figure_axis_robust(14,14)
        hold off;

        saveas(f, fullfile(fig_dir, 'CIPL_Strategy', ...
            sprintf('CIPL_%s_%s.png', groupNames{s}, sgFileLabel)));
        close;
    end
end

% ------- Sex and Genotype - one plot ------

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
    pubify_figure_axis_robust(14,14)
    legend(uniqueSG, 'Location', 'best');
    hold off;

    saveas(f, fullfile(fig_dir, 'CIPL_Strategy', ...
        sprintf('CIPL_%s_SexGenotypeCombined.png', groupNames{s})));
    close;
end

%% 5) Strategy-Day Plots (individual and group)

groupBy_S3 = 'SexGeno';

apaTbl = table( ...
    strings(0,1), strings(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    'VariableNames', {'Strategy','Effect','df','SS','MS','F','p','eta2'});

levels_S3 = unique(data1.(groupBy_S3));
nLvls_S3  = numel(levels_S3);
if nLvls_S3 <= size(baseCols,1)
    clrMap_S3 = num2cell(baseCols(1:nLvls_S3,:), 2);
else
    c = hsv(nLvls_S3); clrMap_S3 = num2cell(c,2);
end

for s = 1:nStrategies

    stratProb  = data1.(strategyNames{s});
    uniqueDays = unique(data1.x_Day);
    uniqueRats = unique(data1.x_TargetID);

    % --- Build mean_strat without plotting individual rats --- %
    mean_strat  = {};
    validRatIdx = 0;

    for r = 1:numel(uniqueRats)
        ratID        = uniqueRats{r};
        isCurrentRat = strcmp(data1.x_TargetID, ratID);
        ratLvl       = data1.(groupBy_S3){find(isCurrentRat,1)};

        meanUse = nan(1, numel(uniqueDays));
        for d = 1:numel(uniqueDays)
            idx = isCurrentRat & (data1.x_Day == uniqueDays(d));
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
        mean_strat{validRatIdx, 2} = ratLvl;
        mean_strat(validRatIdx, 3:6) = num2cell(meanUse);
    end

    mean_strat_table = cell2table(mean_strat, 'VariableNames', ...
        {'RatID', 'Age', 'Day1', 'Day2', 'Day3', 'Day4'});

    anovaResults   = runMixedANOVA(mean_strat_table, {'Day1','Day2','Day3','Day4'});
    postHocResults = runTukeyPostHocMixed(mean_strat_table, {'Day1','Day2','Day3','Day4'});

    writetable(anovaResults,   fullfile(processed_dir, sprintf('Anova_%s.csv',         strategyNames{s})));
    writetable(postHocResults, fullfile(processed_dir, sprintf('PostHoc_Tukey_%s.csv', strategyNames{s})));

    % ---- Within-group day comparisons ---- %
    fprintf('\n========== Within-Group Day Comparisons: %s ==========\n', strategyNames{s});

    sgGroups_wg   = unique(data1.SexGeno);
    uniqueDays_wg = unique(data1.x_Day);
    nDays_wg      = numel(uniqueDays_wg);
    dayColNames_wg = arrayfun(@(d) sprintf('Day%d',d), uniqueDays_wg, 'UniformOutput', false);

    for sgi = 1:numel(sgGroups_wg)
        sg_wg   = sgGroups_wg(sgi);
        sgLabel = strrep(char(sg_wg), '_', '-');

        sgMask_wg = strcmp(data1.SexGeno, sg_wg);
        rats_wg   = unique(data1.x_TargetID(sgMask_wg));
        nRats_wg  = numel(rats_wg);
        if nRats_wg < 3, continue; end

        ratDay_wg = nan(nRats_wg, nDays_wg);
        for r = 1:nRats_wg
            rIdx = strcmp(data1.x_TargetID, rats_wg{r});
            for d = 1:nDays_wg
                dIdx = rIdx & (data1.x_Day == uniqueDays_wg(d));
                if any(dIdx)
                    ratDay_wg(r,d) = mean(stratProb(dIdx));
                end
            end
        end

        nanRows_wg = any(isnan(ratDay_wg), 2);
        ratDay_wg  = ratDay_wg(~nanRows_wg, :);
        if size(ratDay_wg,1) < 3, continue; end

        wTbl_wg  = array2table(ratDay_wg, 'VariableNames', dayColNames_wg);
        WD_wg    = table(uniqueDays_wg, 'VariableNames', {'Day'});
        frm_wg   = sprintf('%s-%s ~ 1', dayColNames_wg{1}, dayColNames_wg{end});
        rm_wg    = fitrm(wTbl_wg, frm_wg, 'WithinDesign', WD_wg);
        anova_wg = ranova(rm_wg);
        ph_day_wg = multcompare(rm_wg, 'Day', 'ComparisonType', 'tukey-kramer');

        fprintf('\n--- %s | %s (n=%d) ---\n', strategyNames{s}, sgLabel, size(ratDay_wg,1));
        disp(anova_wg)
        sig_wg = ph_day_wg(ph_day_wg.pValue < 0.05, :);
        if isempty(sig_wg)
            fprintf('  No significant day differences\n');
        else
            disp(sig_wg)
        end
    end

    % ---- APA table rows ---- %
    rowNames = string(anovaResults.Properties.RowNames);

    if ~strcmpi(strategyNames{s},'perseverance')
        wanted  = ["Age","(Intercept):Day","Age:Day"];
        tidyLab = ["Group","Day","Group×Day"];

        keepIdx = find(ismember(rowNames, wanted));
        if isempty(keepIdx), warning('%s has no mixed-ANOVA rows', strategyNames{s}); end
        L = numel(keepIdx);

        SS   = anovaResults.SumSq(keepIdx);
        df   = anovaResults.DF(keepIdx);
        MS   = anovaResults.MeanSq(keepIdx);
        Fval = anovaResults.F(keepIdx);
        pVal = anovaResults.pValue(keepIdx);

        eta2 = nan(L,1);
        for j = 1:L
            k    = keepIdx(j);
            eRow = find(startsWith(rowNames(k+1:end),"Error"),1,'first') + k;
            if ~isempty(eRow)
                eta2(j) = SS(j) / (SS(j) + anovaResults.SumSq(eRow));
            end
        end

        tmp = table( ...
            repmat(string(strategy_titles{s}), L, 1), ...
            tidyLab(ismember(wanted, rowNames(keepIdx))).', ...
            df, SS, MS, Fval, pVal, eta2, ...
            'VariableNames', {'Strategy','Effect','df','SS','MS','F','p','eta2'});

        apaTbl = [apaTbl; tmp];
    end

    % ---- Group mean per day plot with pairwise significance markers ---- %
    dataMeanGrps = cell(nLvls_S3, numel(uniqueDays));
    for li = 1:nLvls_S3
        lvl = levels_S3(li);
        for d = 1:numel(uniqueDays)
            colName = sprintf('Day%d', d);
            dataMeanGrps{li,d} = mean_strat_table{strcmp(mean_strat_table.Age, lvl), colName};
        end
    end

    % --- Compute pairwise t-tests between groups for each day --- %
    % Rows = days, stored as [groupA, groupB, p-value]
    sigResults = [];   % will hold [dayVal, liA, liB, pVal]
    for d = 1:numel(uniqueDays)
        for liA = 1:nLvls_S3
            for liB = (liA+1):nLvls_S3
                vA = dataMeanGrps{liA, d};
                vB = dataMeanGrps{liB, d};
                if numel(vA) > 1 && numel(vB) > 1
                    [~, pv] = ttest2(vA, vB, 'Vartype', 'unequal');
                    if pv < 0.05
                        sigResults(end+1, :) = [uniqueDays(d), liA, liB, pv];
                    end
                end
            end
        end
    end

    f3 = figure('Position', [100 100 650 520]);
    hold on;

    % Plot error bar lines for each group
    hLines = gobjects(nLvls_S3, 1);
    for li = 1:nLvls_S3
        meanVals = cellfun(@mean, dataMeanGrps(li,:));
        semVals  = cellfun(@(x) std(x)/sqrt(numel(x)), dataMeanGrps(li,:));
        hLines(li) = errorbar(uniqueDays, meanVals, semVals, '-o', ...
            'Color',           clrMap_S3{li}, ...
            'LineWidth',       3, ...
            'MarkerFaceColor', clrMap_S3{li}, ...
            'MarkerEdgeColor', 'none', ...
            'MarkerSize',      8, ...
            'DisplayName',     char(levels_S3(li)));
    end

    % Add significance markers
    % Stack multiple comparisons on the same day at increasing heights
    if ~isempty(sigResults)
        yTop    = 0.55;          % top of your ylim
        yStep   = 0.04;          % vertical spacing between stacked brackets
        barHt   = 0.015;         % height of the bracket tick marks
        xJitter = 0.04;          % slight x offset so brackets don't overlap

        % Count how many comparisons exist per day to set starting height
        for d = 1:numel(uniqueDays)
            dayRows = sigResults(sigResults(:,1) == uniqueDays(d), :);
            if isempty(dayRows), continue; end

            nComp = size(dayRows, 1);
            for ci = 1:nComp
                liA  = dayRows(ci, 2);
                liB  = dayRows(ci, 3);
                pv   = dayRows(ci, 4);
                dVal = uniqueDays(d);

                % y position for this bracket — stack upward
                yBracket = yTop - (nComp - ci) * yStep;

                % x positions: slightly offset per comparison pair
                xA = dVal - xJitter * (nComp - ci);
                xB = dVal + xJitter * (nComp - ci);

                % Draw bracket
                plot([xA xB], [yBracket yBracket], '-', ...
                    'Color', [0.3 0.3 0.3], 'LineWidth', 1.5, ...
                    'HandleVisibility', 'off');
                plot([xA xA], [yBracket - barHt, yBracket], '-', ...
                    'Color', [0.3 0.3 0.3], 'LineWidth', 1.5, ...
                    'HandleVisibility', 'off');
                plot([xB xB], [yBracket - barHt, yBracket], '-', ...
                    'Color', [0.3 0.3 0.3], 'LineWidth', 1.5, ...
                    'HandleVisibility', 'off');

                % Significance symbol
                if pv < 0.001
                    sigSym = '***';
                elseif pv < 0.01
                    sigSym = '**';
                else
                    sigSym = '*';
                end
                text(mean([xA xB]), yBracket + 0.005, sigSym, ...
                    'HorizontalAlignment', 'center', ...
                    'FontSize', 14, 'FontWeight', 'bold', ...
                    'Color', [0.3 0.3 0.3], ...
                    'HandleVisibility', 'off');
            end
        end
    end

    xlabel('Day');
    ylabel('Probability of Strategy Use');
    title(strategy_titles{s});
    legend('Location', 'northeastoutside');
    ylim([0 0.62]);    % slightly above 0.55 to give room for brackets
    xlim([0.75 4.25]);
    xticks(uniqueDays);
    pubify_figure_axis_robust(16, 18);
    hold off;

    exportgraphics(f3, fullfile(fig_dir, 'StrategyUse', ...
        sprintf('1_RatMeanPerDay_%s.png', strategyNames{s})), 'Resolution', 450);
    close

end

writetable(apaTbl, fullfile(processed_dir, 'ANOVA_APA_AllStrategies.csv'));
disp(apaTbl)
saveTablePNG_APA(apaTbl, ...
    fullfile(fig_dir, 'ANOVA_APA_AllStrategies.png'), ...
    'Title', 'Repeated Measures ANOVA for Strategies');

%% 6) Strategy by Group - Non-spatial/Procedural/Allocentric

groupBy_S4 = 'SexGeno';   % <-- CHANGE THIS PER ANALYSIS

strategyGroups = {
    {'thigmotaxis', 'circling', 'randomPath'}, 'Platform-Independent';
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
uniqueDays = unique(data1.x_Day);

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
            dayIdx = isCurrentRat & (data1.x_Day == uniqueDays(d));
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
    pubify_figure_axis_robust(14,14)
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
        set(gca, 'FontSize', 14, 'FontWeight', 'bold', ...
            'LineWidth', 1.5, 'Box', 'off', 'TickDir', 'out');        hold off;
        saveas(f3, fullfile(fig_dir, 'StrategyUse',sprintf('RatMeanPerDay_%s', groupNames{g})), 'png');
        close;
    else
        f3 = figure; hold on;
        for li = 1:nLvls_S4
            meanVals = cellfun(@mean, dataMeanGrps(li,:));
            semVals  = cellfun(@(x) std(x)/sqrt(numel(x)), dataMeanGrps(li,:));
            errorbar(uniqueDays, meanVals, semVals, '-o', ...
                'Color', clrMap_S4{li}, 'LineWidth', 5, ...
                'DisplayName', char(levels_S4(li)));
        end
        xlabel('Day'); ylabel('Probability of Strategy Use');
        title(sprintf('%s: Rat Mean Per Day', groupNames{g}), 'FontSize', 16);
        legend('Location','northeastoutside');
        ylim([0 1.25]);
        xticks(uniqueDays);
        set(gca, 'FontSize', 14, 'FontWeight', 'bold', ...
            'LineWidth', 1.5, 'Box', 'off', 'TickDir', 'out');
        hold off;
        saveas(f3, fullfile(fig_dir, 'StrategyUse',sprintf('RatMeanPerDay_%s', groupNames{g})), 'png');
        close;
    end
end

%% 7) Entropy Calculation

groupBy_S5 = 'SexGeno';   % <-- CHANGE THIS PER ANALYSIS

% Derive levels and colours for chosen factor
levels_S5 = unique(data1.(groupBy_S5));
nLvls_S5  = numel(levels_S5);
if nLvls_S5 <= size(baseCols,1)
    clrMap_S5 = num2cell(baseCols(1:nLvls_S5,:), 2);
else
    c = hsv(nLvls_S5); clrMap_S5 = num2cell(c,2);
end

%------------------ Compute Entropy for Each Trial ------------------%
% Index strategy columns by name rather than position (position shifts
% when new columns like Sex/APP/AgeGroup are added to data1)
stratCols = {'thigmotaxis','circling','randomPath','scanning',...
    'chaining','directedSearch','correctedPath','directPath','perseverance'};
pMat = data1{:, stratCols};
entropy_vals = -sum(pMat .* log2(pMat + eps), 2);
data1.entropy = entropy_vals;

uniqueRats = unique(data1.x_TargetID);
uniqueDays = unique(data1.x_Day);
nRats = numel(uniqueRats);

meanDayEntropy = nan(nRats, numel(uniqueDays));
AgeCell = cell(nRats,1);   % stores the groupBy_S5 label per rat

for r = 1:nRats
    ratID = uniqueRats{r};
    idxRat = strcmp(data1.x_TargetID, ratID);
    AgeCell{r} = data1.(groupBy_S5){find(idxRat,1)};
    for d = 1:numel(uniqueDays)
        idxDay = idxRat & (data1.x_Day == uniqueDays(d));
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
        'MarkerEdgeColor', 'none', 'HandleVisibility', 'off');
end

for li = 1:nLvls_S5
    ageMask = strcmp(AgeCell, levels_S5(li));
    grpMean = mean(meanDayEntropy(ageMask,:), 1, 'omitnan');
    plot(uniqueDays, grpMean, '-o', ...
        'Color', clrMap_S5{li}, 'LineWidth', lwMean, ...
        'MarkerSize', 6, 'MarkerFaceColor', clrMap_S5{li}, ...
        'DisplayName', strrep(char(levels_S5(li)), '_', '-'));
end

xlabel('Day'); ylabel('Entropy');
title('Entropy Across All 8 Strategy Types');
xticks(uniqueDays);
legend('Location','Northeast');
pubify_figure_axis_robust(14,14);
hold off;
saveas(fe3, fullfile(fig_dir, 'Entropy', 'Entropy_All8_MeanRat'), 'png');
close(fe3);


% ---- Sex×Genotype entropy figure for Section 5 ---- %
% One figure with one line per Sex×Genotype combo (M-WT, M-APP/+, F-WT, F-APP/+).
uniqueDays_sg5 = unique(data1.x_Day);

% Build Sex×Genotype label per rat
if ~ismember('SexGenotype', data1.Properties.VariableNames)
    data1.SexGenotype = data1.Sex + "-" + data1.APP;
end
[uniqueRats_sg5, ia_sg5] = unique(data1.x_TargetID);
sgLabel5 = data1.SexGenotype(ia_sg5);
[uniqueSG5, ~, sgIdx5] = unique(sgLabel5);
nSG5 = numel(uniqueSG5);

if nSG5 <= size(baseCols,1)
    sgCols5 = num2cell(baseCols(1:nSG5,:), 2);
else
    c = hsv(nSG5); sgCols5 = num2cell(c, 2);
end

% Compute per-rat mean entropy per day for each Sex×Genotype group
sgRatMeans5 = cell(1, nSG5);
for sgi = 1:nSG5
    rats_sg  = uniqueRats_sg5(sgIdx5 == sgi);
    nRats_sg = numel(rats_sg);
    ratE     = nan(nRats_sg, numel(uniqueDays_sg5));
    for r = 1:nRats_sg
        rIdx = strcmp(data1.x_TargetID, rats_sg{r});
        for d = 1:numel(uniqueDays_sg5)
            dIdx = rIdx & (data1.x_Day == uniqueDays_sg5(d));
            if any(dIdx)
                ratE(r,d) = mean(data1.entropy(dIdx));
            end
        end
    end
    sgRatMeans5{sgi} = ratE;
end

% Single figure: all Sex×Genotype groups as separate lines
fe_sg5 = figure; hold on;
for sgi = 1:nSG5
    ratE     = sgRatMeans5{sgi};
    nRats_sg = size(ratE, 1);
    meanE_sg = mean(ratE, 1, 'omitnan');
    if nRats_sg > 1
        semE_sg = std(ratE, 0, 1, 'omitnan') ./ sqrt(nRats_sg);
    else
        semE_sg = zeros(size(meanE_sg));
    end

    errorbar(uniqueDays_sg5, meanE_sg, semE_sg, '-o', ...
        'Color', sgCols5{sgi}, 'LineWidth', 2, ...
        'MarkerFaceColor', sgCols5{sgi}, 'MarkerEdgeColor', 'none', ...
        'MarkerSize', 7, ...
        'DisplayName', sprintf('%s (n=%d)', char(uniqueSG5(sgi)), nRats_sg));
end

title('Mean Entropy by Sex × Genotype', 'FontSize', 16, 'FontWeight', 'bold');
xlabel('Day', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Entropy', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'northeast');
xticks(uniqueDays_sg5);
set(gca, 'FontSize', 14, 'FontWeight', 'bold', ...
    'LineWidth', 1.5, 'Box', 'off', 'TickDir', 'out');
hold off;

saveas(fe_sg5, fullfile(fig_dir, 'Entropy', 'Entropy_MeanRat_SexGenotype.png'));
close;

%% 8) Strategy Group Probability Bar Plots by Sex×Genotype across Days

strategyGroups_sg = {
    {'thigmotaxis','circling','randomPath'},                          'Platform-Independent';
    {'scanning','chaining'},                                          'Procedural';
    {'directedSearch','correctedPath','directPath','perseverance'},   'Allocentric'
    };
uniqueDays_sg = unique(data1.x_Day);
nDays_sg      = numel(uniqueDays_sg);

% Get SexGeno groups
sgGroups_bar  = unique(data1.SexGeno);
nSG_bar       = numel(sgGroups_bar);
if nSG_bar <= size(baseCols,1)
    sgCols_bar = baseCols(1:nSG_bar,:);
else
    sgCols_bar = hsv(nSG_bar);
end

for s = 1:size(strategyGroups_sg,1)
    % Compute group probability for this strategy group
    currentStrats = strategyGroups_sg{s,1};
    groupProb_sg  = data1.(currentStrats{1});
    for ii = 2:numel(currentStrats)
        groupProb_sg = groupProb_sg + data1.(currentStrats{ii});
    end
    sName = strategyGroups_sg{s,2};

    % Build per-rat daily means
    uniqueRats_sg = unique(data1.x_TargetID);
    nRats_sg      = numel(uniqueRats_sg);
    ratDayProb    = nan(nRats_sg, nDays_sg);
    ratSG         = strings(nRats_sg,1);
    for r = 1:nRats_sg
        rIdx = strcmp(data1.x_TargetID, uniqueRats_sg{r});
        ratSG(r) = data1.SexGeno(find(rIdx,1));
        for d = 1:nDays_sg
            dIdx = rIdx & (data1.x_Day == uniqueDays_sg(d));
            if any(dIdx)
                ratDayProb(r,d) = mean(groupProb_sg(dIdx));
            end
        end
    end

    % Compute mean and SEM per SexGeno per day
    meanData_sg = nan(nDays_sg, nSG_bar);
    semData_sg  = nan(nDays_sg, nSG_bar);
    ratPts_sg   = cell(nDays_sg, nSG_bar);
    for sgi = 1:nSG_bar
        mask = ratSG == sgGroups_bar(sgi);
        nR   = sum(mask);
        for d = 1:nDays_sg
            vals = ratDayProb(mask, d);
            vals = vals(~isnan(vals));
            meanData_sg(d,sgi) = mean(vals);
            semData_sg(d,sgi)  = std(vals) / sqrt(numel(vals));
            ratPts_sg{d,sgi}   = vals;
        end
    end

    % Plot
    f_sg = figure('Position', [95, 100, 1100, 600]); hold on;
    hBar_sg = bar(uniqueDays_sg, meanData_sg, 'grouped', 'BarWidth', 1);
    for sgi = 1:nSG_bar
        hBar_sg(sgi).FaceColor = sgCols_bar(sgi,:);
        hBar_sg(sgi).FaceAlpha = 0.5;
        hBar_sg(sgi).DisplayName = strrep(char(sgGroups_bar(sgi)), '_', '-');
    end
    drawnow;

    % Bar centers for error bars and scatter
    barCtrs_sg = nan(nDays_sg, nSG_bar);
    for sgi = 1:nSG_bar
        barCtrs_sg(:,sgi) = hBar_sg(sgi).XEndPoints';
    end

    % Error bars
    for sgi = 1:nSG_bar
        errorbar(barCtrs_sg(:,sgi), meanData_sg(:,sgi), semData_sg(:,sgi), ...
            'k', 'LineStyle', 'none', 'LineWidth', 3, 'HandleVisibility', 'off');
    end

    % Jittered individual rat points
    jit_sg = 0.06;
    for sgi = 1:nSG_bar
        for d = 1:nDays_sg
            pts = ratPts_sg{d,sgi};
            if isempty(pts), continue; end
            xj = barCtrs_sg(d,sgi) + (rand(size(pts))-0.5)*jit_sg;
            scatter(xj, pts, 20, 'MarkerFaceColor', sgCols_bar(sgi,:), ...
                'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.7, ...
                'HandleVisibility', 'off');
        end
    end

    title(sprintf('%s', sName), ...
        'FontSize', 96, 'FontWeight', 'bold');
    xlabel('Day');
    ylabel('Mean Probability');
    xticks(uniqueDays_sg);
    ylim([0 1.1]);
    legend(hBar_sg, 'Location', 'northeastoutside', 'Interpreter', 'none', 'FontSize', 18);    set(gca, 'FontSize', 18, 'FontWeight', 'bold', ...
        'LineWidth', 3, 'Box', 'off', 'TickDir', 'out');
    pubify_figure_axis_robust(21,21)
    hold off;

    % ---- Significance lines: WT vs APP+ within each sex, per day ---- %
    idx_MWT  = find(sgGroups_bar == "M-WT",  1);
    idx_MAPP = find(sgGroups_bar == "M-APP+", 1);
    idx_FWT  = find(sgGroups_bar == "F-WT",  1);
    idx_FAPP = find(sgGroups_bar == "F-APP+", 1);

    sigPairs_bar = {};
    sigP_bar     = [];
    comparisons_bar = {idx_MWT, idx_MAPP; idx_FWT, idx_FAPP};
    for ci = 1:size(comparisons_bar,1)
        ia = comparisons_bar{ci,1};
        ib = comparisons_bar{ci,2};
        if isempty(ia) || isempty(ib), continue; end
        for d = 1:nDays_sg
            g1 = ratPts_sg{d, ia};
            g2 = ratPts_sg{d, ib};
            if numel(g1) > 1 && numel(g2) > 1
                [~, pv] = ttest2(g1, g2, 'Vartype', 'unequal');
                if pv < 0.05
                    sigPairs_bar{end+1} = [barCtrs_sg(d,ia), barCtrs_sg(d,ib)];
                    sigP_bar(end+1)     = pv;
                end
            end
        end
    end

    if ~isempty(sigPairs_bar)
        sigstar(sigPairs_bar, sigP_bar);
    end

    hold off;

    saveas(f_sg, fullfile(fig_dir, 'StrategyUse', ...
        sprintf('StrategyProb_%s_bySexGeno.png', strrep(sName,' ','-'))));
    close(f_sg);

    % ---- Stats: 3-way mixed ANOVA (Sex × Genotype × Day) ---- %
    % Build wide table: one row per rat, Day1..Day4 as within factor,
    % Sex and Genotype as between factors.
    wideCell = cell(nRats_sg, 2 + nDays_sg);
    for r = 1:nRats_sg
        rIdx_w = strcmp(data1.x_TargetID, uniqueRats_sg{r});
        wideCell{r,1} = char(data1.Sex(find(rIdx_w,1)));
        wideCell{r,2} = strrep(char(data1.APP(find(rIdx_w,1))), '/', '');
        for d = 1:nDays_sg
            wideCell{r,2+d} = ratDayProb(r,d);
        end
    end
    dayColNames = arrayfun(@(d) sprintf('Day%d',d), uniqueDays_sg, ...
        'UniformOutput', false);
    allVarNames = [{'Sex','Genotype'}, dayColNames(:)'];
    wTbl = cell2table(wideCell, 'VariableNames', allVarNames);
    wTbl.Sex      = string(wTbl.Sex);
    wTbl.Genotype = string(wTbl.Genotype);

    % Remove rows with any NaN day values
    nanRows = any(isnan(wTbl{:, dayColNames}), 2);
    wTbl = wTbl(~nanRows, :);

    WD_sg   = table(uniqueDays_sg, 'VariableNames', {'Day'});
    frm_sg  = sprintf('%s-%s ~ Sex * Genotype', dayColNames{1}, dayColNames{end});
    rm_sg   = fitrm(wTbl, frm_sg, 'WithinDesign', WD_sg);
    anova_sg = ranova(rm_sg, 'WithinModel', 'Day');

    % Post-hoc comparisons
    ph_geno_sg = multcompare(rm_sg, 'Genotype', 'By', 'Day', ...
        'ComparisonType', 'tukey-kramer');
    ph_sex_sg  = multcompare(rm_sg, 'Sex', 'By', 'Day', ...
        'ComparisonType', 'tukey-kramer');
    ph_day_sg  = multcompare(rm_sg, 'Day', 'ComparisonType', 'tukey-kramer');

    % Print to command window
    sNameClean = strrep(sName, ' ', '-');
    fprintf('========== %s: Mixed ANOVA (Sex × Genotype × Day) ==========', sName);
    disp(anova_sg)
    fprintf('--- Significant Genotype × Day post-hoc ---');
    disp(ph_geno_sg(ph_geno_sg.pValue < 0.05, :))
    fprintf('--- Significant Sex × Day post-hoc ---');
    disp(ph_sex_sg(ph_sex_sg.pValue < 0.05, :))

    % Save CSVs
    writetable(anova_sg, fullfile(processed_dir, ...
        sprintf('StrategyProb_%s_ANOVA.csv', sNameClean)), 'WriteRowNames', true);
    writetable(ph_geno_sg, fullfile(processed_dir, ...
        sprintf('StrategyProb_%s_PostHoc_Genotype.csv', sNameClean)));
    writetable(ph_sex_sg, fullfile(processed_dir, ...
        sprintf('StrategyProb_%s_PostHoc_Sex.csv', sNameClean)));
    writetable(ph_day_sg, fullfile(processed_dir, ...
        sprintf('StrategyProb_%s_PostHoc_Day.csv', sNameClean)));

    fprintf('Stats saved for %s\n', sName);
end
%% 9) Stats for biggest strat differences figs

% ---- Stats: Individual Strategy Use by Sex × Genotype × Day ---- %
% 3-way mixed ANOVA for correctedPath, directedSearch, thigmotaxis,
% chaining
individualStrats = {'correctedPath', 'directedSearch', 'thigmotaxis', 'chaining'};
stratLabels_stat = {'Corrected Path', 'Directed Search', 'Thigmotaxis', 'Chaining'};

uniqueDays_is  = unique(data1.x_Day);
nDays_is       = numel(uniqueDays_is);
uniqueRats_is  = unique(data1.x_TargetID);
nRats_is       = numel(uniqueRats_is);

for si = 1:numel(individualStrats)
    stratName  = individualStrats{si};
    stratLabel = stratLabels_stat{si};

    % Per-rat daily mean probability for this strategy
    ratDayProb_is = nan(nRats_is, nDays_is);
    for r = 1:nRats_is
        rIdx = strcmp(data1.x_TargetID, uniqueRats_is{r});
        for d = 1:nDays_is
            dIdx = rIdx & (data1.x_Day == uniqueDays_is(d));
            if any(dIdx)
                ratDayProb_is(r,d) = mean(data1.(stratName)(dIdx));
            end
        end
    end

    % Build wide table with Sex and Genotype as between factors
    dayColNames_is = arrayfun(@(d) sprintf('Day%d',d), uniqueDays_is, 'UniformOutput', false);
    wideCell_is = cell(nRats_is, 2 + nDays_is);
    for r = 1:nRats_is
        rIdx_w = strcmp(data1.x_TargetID, uniqueRats_is{r});
        wideCell_is{r,1} = char(data1.Sex(find(rIdx_w,1)));
        wideCell_is{r,2} = strrep(char(data1.APP(find(rIdx_w,1))), '/', '');
        for d = 1:nDays_is
            wideCell_is{r,2+d} = ratDayProb_is(r,d);
        end
    end
    allVarNames_is = [{'Sex','Genotype'}, dayColNames_is(:)'];
    wTbl_is = cell2table(wideCell_is, 'VariableNames', allVarNames_is);
    wTbl_is.Sex      = string(wTbl_is.Sex);
    wTbl_is.Genotype = string(wTbl_is.Genotype);

    % Remove rows with any NaN
    nanRows_is = any(isnan(wTbl_is{:, dayColNames_is}), 2);
    wTbl_is = wTbl_is(~nanRows_is, :);

    % Fit repeated-measures model
    WD_is  = table(uniqueDays_is, 'VariableNames', {'Day'});
    frm_is = sprintf('%s-%s ~ Sex * Genotype', dayColNames_is{1}, dayColNames_is{end});
    rm_is  = fitrm(wTbl_is, frm_is, 'WithinDesign', WD_is);
    anova_is = ranova(rm_is, 'WithinModel', 'Day');

    % Post-hoc
    ph_geno_is = multcompare(rm_is, 'Genotype', 'By', 'Day', 'ComparisonType', 'tukey-kramer');
    ph_sex_is  = multcompare(rm_is, 'Sex',      'By', 'Day', 'ComparisonType', 'tukey-kramer');
    ph_day_is  = multcompare(rm_is, 'Day',             'ComparisonType', 'tukey-kramer');

    % Print to command window
    fprintf('\n========== %s: Mixed ANOVA (Sex x Genotype x Day) ==========\n', stratLabel);
    disp(anova_is)
    fprintf('--- Significant Genotype x Day post-hoc ---\n');
    sig_geno = ph_geno_is(ph_geno_is.pValue < 0.05, :);
    if isempty(sig_geno), fprintf('  none\n'); else, disp(sig_geno); end
    fprintf('--- Significant Sex x Day post-hoc ---\n');
    sig_sex = ph_sex_is(ph_sex_is.pValue < 0.05, :);
    if isempty(sig_sex), fprintf('  none\n'); else, disp(sig_sex); end
    fprintf('--- Significant Day post-hoc ---\n');
    sig_day = ph_day_is(ph_day_is.pValue < 0.05, :);
    if isempty(sig_day), fprintf('  none\n'); else, disp(sig_day); end

    % Save CSVs
    stratFile = strrep(stratLabel, ' ', '-');
    writetable(anova_is,   fullfile(processed_dir, sprintf('IndivStrat_%s_ANOVA.csv',            stratFile)), 'WriteRowNames', true);
    writetable(ph_geno_is, fullfile(processed_dir, sprintf('IndivStrat_%s_PostHoc_Genotype.csv', stratFile)));
    writetable(ph_sex_is,  fullfile(processed_dir, sprintf('IndivStrat_%s_PostHoc_Sex.csv',      stratFile)));
    writetable(ph_day_is,  fullfile(processed_dir, sprintf('IndivStrat_%s_PostHoc_Day.csv',      stratFile)));

    fprintf('Stats saved for %s\n', stratLabel);
end

%% 10) Group Strategy Based Entropy

groupBy_S6 = 'SexGeno';   % <-- CHANGE THIS PER ANALYSIS

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
uniqueDays = unique(data1.x_Day);
nRats = numel(uniqueRats);

meanDayEntropy = nan(nRats, numel(uniqueDays));
AgeCell = cell(nRats,1);

for r = 1:nRats
    ratID = uniqueRats{r};
    idxRat = strcmp(data1.x_TargetID, ratID);
    AgeCell{r} = data1.(groupBy_S6){find(idxRat,1)};
    for d = 1:numel(uniqueDays)
        idxDay = idxRat & (data1.x_Day == uniqueDays(d));
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
set(gca, 'FontSize', 14, 'FontWeight', 'bold', ...
    'LineWidth', 1.5, 'Box', 'off', 'TickDir', 'out');
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
    set(gca, 'FontSize', 14, 'FontWeight', 'bold', ...
        'LineWidth', 1.5, 'Box', 'off', 'TickDir', 'out');
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
    set(gca, 'FontSize', 14, 'FontWeight', 'bold', ...
        'LineWidth', 1.5, 'Box', 'off', 'TickDir', 'out');
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
set(gca, 'FontSize', 14, 'FontWeight', 'bold', ...
    'LineWidth', 1.5, 'Box', 'off', 'TickDir', 'out');
xticks(uniqueDays);
hold off;
saveas(fe3, fullfile(fig_dir, 'Entropy', 'GroupStrat_Entropy_RatChanges'), 'png');
close;


% -----------------------------------------------------------------------
% NEW: Set groupBy to choose which factor defines the subplot panels.
%   Each unique level of groupBy_S9 gets its own figure.
% -----------------------------------------------------------------------
groupBy_S9 = 'SexGeno';   % <-- CHANGE THIS PER ANALYSIS

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
uniqueDays   = unique(data1.x_Day);

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
            ratMeans(i, col) = mean(gp(data1.x_Day(ridx)==uniqueDays(d)));
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
            dVal = double(cs.x_Day(rr));
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
        sel = double(cd9.x_Day_1)==uniqueDays(d) & double(cd9.x_Day_2)==uniqueDays(d+1);
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

    title(sprintf('%s', upper(ag)), 'FontSize', 16, 'FontWeight', 'bold');
    xlabel('Day', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Probability for Strategy Group', 'FontSize', 14, 'FontWeight', 'bold');
    xticks(uniqueDays);
    ylim([0 1.2]);
    set(gca, 'FontSize', 14, 'FontWeight', 'bold', ...
        'LineWidth', 1.5, 'Box', 'off', 'TickDir', 'out');
    legend(groupNames9, 'Location', 'northeastoutside', 'FontSize', 12);
    hold off;

    saveas(f, fullfile(fig_dir, sprintf('DayStrategyUse_%s.png', ag)));
    close(f);
end

% ---- Per Sex×Age×APP combo figures for Section 9 ---- %
% One figure per combo: all three strategy groups as grouped bars across days.
% No stats — descriptive only, since each combo may have very few rats.
strategyGroups_c9 = {
    {'thigmotaxis','circling','random_path'}, 'NonGoal',    'Non-Goal Oriented';
    {'scanning','chaining'},                 'Procedural',  'Procedural';
    {'directed_search','corrected_search','direct_path','perseverance'}, 'Allocentric', 'Allocentric'
    };
stratColors_c9 = [
    0.3961, 0.2627, 0.1294;
    1.0000, 0.7020, 0.4000;
    0,      0.5020, 0
    ];
nGroups_c9   = size(strategyGroups_c9, 1);
uniqueDays_c9 = unique(data1.x_Day);
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
                rIdx = strcmp(data1.x_TargetID, rats_c{r}) & idx & (data1.x_Day == dVal);
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

%% 11) Individual strategy differences per day per grouping level

groupBy_S10 = 'SexGeno';   % <-- CHANGE THIS PER ANALYSIS

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
uniqueDays  = unique(data1.x_Day);
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
            dayMask = (data1.x_Day(ridx)==uniqueDays(d));
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


% -----------------------------------------------------------------------
% No groupBy needed here — this section loops within each strategy subgroup
% and already uses Age as the between-subjects factor via fitrm.
% If you want to substitute a different between-subjects factor, change
% the formula string in the fitrm call below (replace 'Age' with e.g. 'Sex').
% -----------------------------------------------------------------------
betweenFactor_S11 = 'SexGeno';   % <-- CHANGE THIS PER ANALYSIS
%   'Age'      → uses the Age column of the rat table (= groupBy used in S3)
%   'Sex'      → add a Sex column to the rat table and use that
%   'APP'      → similarly

strategyGroups11 = { ...
    {'thigmotaxis','circling','randomPath'},'NonGoal'; ...
    {'scanning','chaining'},  'Procedural'; ...
    {'directedSearch','correctedPath','directPath'},  'Allocentric'};

groupPrettyNames11 = {'Platform-Independent','Procedural','Allocentric'};

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

uniqueDays11 = unique(data1.x_Day); nDays11 = numel(uniqueDays11);
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
            dayAgeMask = (data1.x_Day == dVal) & strcmp(data1.(groupBy_S9), ageTag);
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
            dayVec = data1.x_Day(rMask);
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
                mask = (data1.x_Day==uniqueDays11(d)) & strcmp(data1.(groupBy_S9), ageTag11);
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
    {'thigmotaxis','circling','randomPath'}, 'NonGoal',   'Platform-Independent'; ...
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

uniqueDays_c11 = unique(data1.x_Day);
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
                    rIdx = strcmp(data1.x_TargetID, rats_c{r}) & idx & (data1.x_Day == dVal);
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
        set(gca, 'FontSize', 14, 'FontWeight', 'bold', ...
            'LineWidth', 1.5, 'Box', 'off', 'TickDir', 'out');
        hold off;

        saveas(f_c11, fullfile(fig_dir, ...
            sprintf('StrategyDay_%s_%s.png', subgroupTag_c11, comboList(ci).label)));
        close(f_c11);
    end
end