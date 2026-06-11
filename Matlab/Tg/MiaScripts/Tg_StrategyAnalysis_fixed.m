% Code for Morris Water Maze - Tg version
% Loads data from cohort xlsx files and writes combined CSVs.
% Robustly handles numeric/string/cell Animal IDs and fills
% Age, Sex, and APP from the Key sheet.

%% GET DATA
base_dir = '/Users/miasponseller/Desktop/Lab/Rtrack/Tg/MatlabFiles';
proc_dir = fullfile(base_dir,'Processed');
if ~exist(proc_dir,"dir"), mkdir(proc_dir); end

files = dir(fullfile(base_dir,'*.xlsx'));

AllSpatial = table();
AllProbe   = table();

for i = 1:numel(files)
    fname = fullfile(base_dir, files(i).name);

    %--- Parse cohort label AND sex from filename: e.g. "Coh10M_Spatial.xlsx"
    %    Cohort is stored as the full string "10M" / "10F" (not just the number).
    %    Pattern: Coh<number><M or F>
    matchCoh = regexp(files(i).name, 'Coh(\d+[MmFf])', 'tokens', 'once');
    if ~isempty(matchCoh)
        cohortStr   = upper(matchCoh{1});           % e.g. '10M', '1F'
        sexFromFile = upper(cohortStr(end));         % last char is M or F
    else
        % Fallback: just grab the number with no sex letter
        matchNum = regexp(files(i).name, 'Coh(\d+)', 'tokens', 'once');
        if ~isempty(matchNum)
            cohortStr = matchNum{1};
        else
            cohortStr = 'Unknown';
        end
        sexFromFile = '';
        fprintf('\n  WARNING: could not parse sex letter from "%s".', files(i).name);
    end

    %=== Read "Key" sheet
    keyT = readtable(fname,'Sheet','Key','VariableNamesRange',1);

    % Barnes/Cowen ID column
    bnFields = intersect({'BarnesID','Barnes ID','CowenID','Cowen ID'}, ...
                         keyT.Properties.VariableNames);
    if isempty(bnFields)
        warning('No Barnes/Cowen ID column found in %s\nAvailable: %s', ...
            fname, strjoin(keyT.Properties.VariableNames,', '));
        continue   % skip this file — can't match animals without an ID
    end

    % Sex and APP columns (optional — silently absent if not in Key)
    sexField = intersect({'Sex','sex','Gender','gender'}, ...
                         keyT.Properties.VariableNames, 'stable');
    appField = intersect({'APP','App','app','Genotype','genotype'}, ...
                         keyT.Properties.VariableNames, 'stable');

    % Helper: convert any BN value to a plain double
    toNum = @(x) str2double(string(x));

    %=== Read "Spatial" sheet
    %--- Read Spatial sheet then rename mangled column headers before selecting.
    %    Excel headers like "Duration (s)" become "Duration_s_" after MATLAB
    %    sanitises them. We rename them back to plain names right away.
    spT_raw = readtable(fname, 'Sheet', 'Spatial', 'VariableNamesRange', 1, ...
        'VariableNamingRule', 'preserve');
    % Build a safe copy of variable names for programmatic use
    spT = spT_raw;
    spT.Properties.VariableNames = matlab.lang.makeValidName(spT.Properties.VariableNames);

    % Rename mangled names → canonical names used throughout the script
    renamePairs = {
        'Duration_s_',   'Duration';
        'Distance_m_',   'Distance';
        'Duration_s_1',  'Duration';   % in case of duplicate suffix
        'Distance_m_1',  'Distance';
    };
    for rr = 1:size(renamePairs,1)
        if ismember(renamePairs{rr,1}, spT.Properties.VariableNames) && ...
           ~ismember(renamePairs{rr,2}, spT.Properties.VariableNames)
            spT = renamevars(spT, renamePairs{rr,1}, renamePairs{rr,2});
        end
    end

    spCols = {'Test','Animal','Trial','Duration','Distance','PathEfficiency','Platform_CIPL'};
    spKeep = intersect(spCols, spT.Properties.VariableNames);
    spT    = spT(:, spKeep);

    % -----------------------------------------------------------------
    % Standardise Animal to numeric BEFORE the age-filling loop so that
    % BN == Animal comparisons are always type-safe.
    % -----------------------------------------------------------------
    if iscell(spT.Animal)
        spT.Animal = cellfun(@(x) str2double(string(x)), spT.Animal);
    elseif isstring(spT.Animal) || ischar(spT.Animal)
        spT.Animal = str2double(string(spT.Animal));
    end
    % spT.Animal is now a numeric double vector

    % Initialise metadata columns
    spT.Age    = nan(height(spT), 1);
    spT.Sex    = repmat("", height(spT), 1);
    spT.APP    = repmat("", height(spT), 1);
    spT.Cohort = repmat(string(cohortStr), height(spT), 1);

    % -----------------------------------------------------------------
    % Fill Age / Sex / APP from Key sheet
    % -----------------------------------------------------------------
    for k = 1:height(keyT)
        % Get BN as a plain number regardless of Key sheet storage format
        rawBN = keyT.(bnFields{1})(k);
        if iscell(rawBN), rawBN = rawBN{1}; end
        BN = toNum(rawBN);

        matchIdx = (spT.Animal == BN);
        if ~any(matchIdx), continue; end

        % Age — strip units ("mo", "months", etc.) and keep the number
        ageStr = string(keyT.Age(k));
        ageVal = str2double(regexp(ageStr, '\d+\.?\d*', 'match', 'once'));
        spT.Age(matchIdx) = ageVal;

        % Sex — prefer Key sheet value; fall back to sex parsed from filename
        if ~isempty(sexField)
            rawSex = keyT.(sexField{1})(k);
            if iscell(rawSex), rawSex = rawSex{1}; end
            spT.Sex(matchIdx) = string(rawSex);
        elseif ~isempty(sexFromFile)
            spT.Sex(matchIdx) = sexFromFile;
        end

        % APP / genotype
        if ~isempty(appField)
            rawApp = keyT.(appField{1})(k);
            if iscell(rawApp), rawApp = rawApp{1}; end
            spT.APP(matchIdx) = string(rawApp);
        end
    end

    % -----------------------------------------------------------------
    % Stack into AllSpatial — align columns across files
    % -----------------------------------------------------------------
    if isempty(AllSpatial)
        AllSpatial = spT;
    else
        % Add columns present in AllSpatial but missing in spT
        for c = setdiff(AllSpatial.Properties.VariableNames, spT.Properties.VariableNames)
            if isnumeric(AllSpatial.(c{1}))
                spT.(c{1}) = nan(height(spT),1);
            else
                spT.(c{1}) = repmat("", height(spT),1);
            end
        end
        % Add columns present in spT but missing in AllSpatial
        for c = setdiff(spT.Properties.VariableNames, AllSpatial.Properties.VariableNames)
            if isnumeric(spT.(c{1}))
                AllSpatial.(c{1}) = nan(height(AllSpatial),1);
            else
                AllSpatial.(c{1}) = repmat("", height(AllSpatial),1);
            end
        end
        % Reorder spT to match AllSpatial column order then concatenate
        spT        = spT(:, AllSpatial.Properties.VariableNames);
        AllSpatial = [AllSpatial; spT];
    end

    fprintf('\nCohort %s added: %d Spatial rows | Age NaN: %d | Sex filled: %d | APP filled: %d', ...
        cohortStr, height(spT), ...
        sum(isnan(spT.Age)), ...
        sum(spT.Sex ~= ""), ...
        sum(spT.APP ~= ""));
end

fprintf('\n\nDone. Total Spatial rows: %d\n', height(AllSpatial));

% Quick NaN report before writing
fprintf('\n--- NaN / empty summary for AllSpatial ---\n');
for col = {'Age','Cohort','Distance','Duration','Platform_CIPL','Sex','APP'}
    c = col{1};
    if ismember(c, AllSpatial.Properties.VariableNames)
        if isnumeric(AllSpatial.(c))
            n = sum(isnan(AllSpatial.(c)));
            fprintf('  %s: %d NaN  (%.1f%%)\n', c, n, 100*n/height(AllSpatial));
        else
            n = sum(AllSpatial.(c) == "");
            fprintf('  %s: %d empty  (%.1f%%)\n', c, n, 100*n/height(AllSpatial));
        end
    end
end

writetable(AllSpatial, fullfile(proc_dir,'AllMorrisWaterMazeData_Spatial.csv'));
fprintf('\nCSV written to: %s\n', proc_dir);

%% Spatial Trials
fname = fullfile(proc_dir,'AllMorrisWaterMazeData_Spatial.csv');
T = readtable(fname);
T = sortrows(T, ["Animal","Trial"]);
T.Group = repmat("Old", height(T), 1);
T.Group(T.Age < 15) = "Young";
T.Day = floor((T.Trial-1)/6) + 1;

%% First Trial Comparison — Sex x Genotype Groups
first_cipl = T.Platform_CIPL(T.Trial == 1);
first_sex  = T.Sex(T.Trial == 1);
first_app  = T.APP(T.Trial == 1);

% Define groups
groups = {"M","WT"; "M","APP/+"; "F","WT"; "F","APP/+"};
labels = {'Male-WT','Male-APP/+','Female-WT','Female-APP/+'};
colors = [0.1216, 0.4706, 0.7059;   % blue   - Male-WT
          0.8392, 0.1529, 0.1569;   % red    - Male-APP
          0.4157, 0.1059, 0.6039;  % purple   - Female-WT
          0.2196, 0.5569, 0.2353]; % green - Female-APP

groupData = cell(4,1);
for g = 1:4
    mask = strcmpi(first_sex, groups{g,1}) & strcmpi(first_app, groups{g,2});
    groupData{g} = first_cipl(mask);
    fprintf('Group %s: n=%d\n', labels{g}, sum(mask));
end

figure;
hold on
for g = 1:4
    boxchart(g * ones(size(groupData{g})), groupData{g}, ...
        'BoxFaceColor', colors(g,:), 'MarkerColor', colors(g,:), 'LineWidth', 3);
end

% Pairwise t-tests (WT vs APP within each sex)
[~, p_male]   = ttest2(groupData{1}, groupData{2}, 'Vartype','unequal');
[~, p_female] = ttest2(groupData{3}, groupData{4}, 'Vartype','unequal');
fprintf('\nMale   WT vs APP: p = %f\n', p_male);
fprintf('Female WT vs APP: p = %f\n', p_female);

title('CIPL Score on First Trial')
xticks(1:4)
xticklabels(labels)
ylabel('Platform CIPL')
set(gca,'FontSize',18,'FontWeight','bold','LineWidth',3,'Box','off','TickDir','out');

%% Repeated Measures ANOVA — Sex x Genotype x Day (mean CIPL per animal per day)

% --- Step 1: Compute per-animal daily mean CIPL ---
T.SexGeno = strcat(T.Sex, "-", T.APP);   % e.g. "M-WT", "F-APP/+"

animals   = unique(T.Animal);
days      = unique(T.Day);
nAnimals  = numel(animals);
nDays     = numel(days);

% Pre-allocate summary table
animalMeans = table();
animalMeans.Animal  = animals;
animalMeans.SexGeno = repmat("", nAnimals, 1);

dayMatrix = nan(nAnimals, nDays);   % rows = animals, cols = days

for a = 1:nAnimals
    idx = T.Animal == animals(a);
    animalMeans.SexGeno(a) = unique(T.SexGeno(idx));   % should be one value per animal

    for d = 1:nDays
        trialIdx = idx & T.Day == days(d);
        if any(trialIdx)
            dayMatrix(a, d) = mean(T.Platform_CIPL(trialIdx), 'omitnan');
        end
    end
end

% Attach day columns to the summary table
for d = 1:nDays
    animalMeans.(sprintf('Day%d', days(d))) = dayMatrix(:, d);
end

fprintf('\n=== Animals per group ===\n');
grps = unique(animalMeans.SexGeno);
for g = 1:numel(grps)
    fprintf('  %s: n=%d\n', grps(g), sum(animalMeans.SexGeno == grps(g)));
end

% --- Step 2: Remove animals with any missing day (rm-ANOVA requires complete cases) ---
completeCases = all(~isnan(dayMatrix), 2);
fprintf('\n%d animals removed due to missing days\n', sum(~completeCases));
animalMeans  = animalMeans(completeCases, :);
dayMatrix    = dayMatrix(completeCases, :);

% --- Step 3: Run repeated measures ANOVA ---
dayVarNames = arrayfun(@(d) sprintf('Day%d', d), days, 'UniformOutput', false);
withinTable = table(days, 'VariableNames', {'Day'});

rm = fitrm(animalMeans, ...
    sprintf('%s-%s ~ SexGeno', dayVarNames{1}, dayVarNames{end}), ...
    'WithinDesign', withinTable);

fprintf('\n=== Repeated Measures ANOVA ===\n');
ranovatbl = ranova(rm, 'WithinModel', 'Day');
disp(ranovatbl)

fprintf('\n=== Between-subjects effect (SexGeno) ===\n');
disp(rm.Coefficients)

% --- Step 4: Mauchly's test for sphericity ---
fprintf('\n=== Mauchly''s Test for Sphericity ===\n');
disp(mauchly(rm))

% --- Step 5: Plot group means +/- SEM across days ---
groupList  = {"M-WT","M-APP/+","F-WT","F-APP/+"};
plotLabels = {'Male-WT','Male-APP/+','Female-WT','Female-APP/+'};
colors     = [0.2 0.4 0.8; 0.1 0.6 0.9; 0.85 0.3 0.5; 0.6 0.0 0.4];

figure; hold on;
for g = 1:4
    idx  = animalMeans.SexGeno == groupList{g};
    vals = dayMatrix(idx, :);
    m    = mean(vals, 1, 'omitnan');
    sem  = std(vals, 0, 1, 'omitnan') ./ sqrt(sum(idx));

    errorbar(days, m, sem, '-o', ...
        'Color', colors(g,:), 'MarkerFaceColor', colors(g,:), ...
        'LineWidth', 1.8, 'MarkerSize', 6, 'DisplayName', plotLabels{g});
end

xlabel('Day'); ylabel('Mean Platform CIPL');
title('Spatial Learning — Sex x Genotype x Day');
legend('Location','best');
xticks(days);
set(gca,'FontSize',14,'FontWeight','bold','LineWidth',1.5,'Box','off','TickDir','out');

%% Pairwise Welch's t-tests — all group combinations
groupPairs = nchoosek(1:4, 2);  % all 6 pairwise combinations

fprintf('\n=== Pairwise Welch''s t-tests (First Trial CIPL) ===\n');
for i = 1:size(groupPairs, 1)
    g1 = groupPairs(i,1);
    g2 = groupPairs(i,2);
    [~, p, ~, stats] = ttest2(groupData{g1}, groupData{g2}, 'Vartype','unequal');
    fprintf('%s vs %s:  t(%0.1f) = %.3f,  p = %.4f\n', ...
        labels{g1}, labels{g2}, stats.df, stats.tstat, p);
end

%% Rest of code
varList = {'Platform_CIPL','Duration','Distance','PathEfficiency'};

for v = 1:numel(varList)
    varName = varList{v};
    disp(['Processing variable: ', varName]);

    G         = groupsummary(T, ["Group","Day"], {"mean","std","nnz"}, varName);
    RatDayAvg = groupsummary(T, ["Animal","Group","Day"], "mean", varName);

    % Group-level average plot
    fig1 = figure('Name',['Average ',varName],'Position',[95,100,800,630]);
    hold on
    grpList = ["Young","Old"];
    clrMap  = {[0 0.6 0],[0.5 0 0.5]};
    for gi = 1:numel(grpList)
        grp    = grpList(gi);
        tmp    = G(G.Group==grp,:);
        semVal = tmp.(['std_',varName]) ./ sqrt(tmp.(['nnz_',varName]));
        errorbar(tmp.Day, tmp.(['mean_',varName]), semVal, '-o', ...
            'Color',clrMap{gi},'LineWidth',3,'DisplayName',grp);
    end
    xlabel('Day','FontSize',12,'FontWeight','bold');
    ylabel(['Average ',varName],'FontSize',12,'FontWeight','bold','Interpreter','none');
    xticks([1 2 3 4]); xlim([0.8 4.2]);
    ymax = max(G.(['mean_',varName]) + 5);
    if ~isnan(ymax), ylim([0 ymax]); end
    legend('FontSize',12);
    title(['Average ',varName,' Scores across Days (Spatial)'],'FontSize',14,'Interpreter','none');
    set(gca,'FontSize',16,'FontWeight','bold','LineWidth',1.5,'Box','off','TickDir','out');
    hold off
    saveas(fig1, fullfile(proc_dir,['AvgMorrisWaterMaze_',varName]),'png');

    % All-animals plot
    fig2 = figure('Name',['All ',varName],'Position',[95,100,800,630]);
    hold on
    ratsYoung = unique(RatDayAvg.Animal(RatDayAvg.Group=="Young"));
    ratsOld   = unique(RatDayAvg.Animal(RatDayAvg.Group=="Old"));
    for r = 1:numel(ratsYoung)
        d = RatDayAvg(RatDayAvg.Animal==ratsYoung(r) & RatDayAvg.Group=="Young",:);
        d = sortrows(d,"Day");
        plot(d.Day, d.(['mean_',varName]),'Color',[clrMap{1},0.175],'LineWidth',2);
    end
    for r = 1:numel(ratsOld)
        d = RatDayAvg(RatDayAvg.Animal==ratsOld(r) & RatDayAvg.Group=="Old",:);
        d = sortrows(d,"Day");
        plot(d.Day, d.(['mean_',varName]),'Color',[clrMap{2},0.175],'LineWidth',2);
    end
    y = plot(G.Day(G.Group=="Young"), G.(['mean_',varName])(G.Group=="Young"), ...
        '-o','Color',clrMap{1},'LineWidth',8,'DisplayName','Young');
    o = plot(G.Day(G.Group=="Old"),   G.(['mean_',varName])(G.Group=="Old"), ...
        '-o','Color',clrMap{2},'LineWidth',8,'DisplayName','Old');
    xlabel('Day','FontSize',14,'FontWeight','bold');
    ylabel('CIPL Score (m\cdots)','FontSize',14,'FontWeight','bold','Interpreter','tex');
    xticks([1 2 3 4]); xlim([0.95 4.05]);
    ymax2 = max(RatDayAvg.(['mean_',varName]));
    if ~isnan(ymax2), ylim([0 ymax2]); end
    legend([y;o],{'Young','Old'},'FontSize',14,'FontWeight','bold');
    set(gca,'FontSize',16,'FontWeight','bold','LineWidth',1.5,'Box','off','TickDir','out');
    title('Performance on Morris watermaze','FontSize',20);
    hold off
    saveas(fig2, fullfile(proc_dir,['MorrisWaterMaze_',varName,'_All']),'png');
end

%% Mixed-design ANOVA + post-hoc for each variable
for v = 1:numel(varList)
    varName = varList{v};
    disp(['Mixed ANOVA for variable: ', varName]);

    RatDayAvg = groupsummary(T, ["Animal","Group","Day"], "mean", varName);
    RatDayAvg.Properties.VariableNames{end} = varName;
    wideTbl = unstack(RatDayAvg, varName, "Day");
    wideTbl.Properties.VariableNames{'Animal'} = 'Subject';
    wideTbl.Properties.VariableNames{'Group'}   = 'Age';

    allCols     = wideTbl.Properties.VariableNames;
    measureVars = setdiff(allCols, {'Subject','Age'}, 'stable');

    anovaResults = runMixedANOVA(wideTbl, measureVars, 'Day');
    writetable(anovaResults, fullfile(base_dir,[varName,'_MixedANOVA_Results.csv']));

    tukeyTbl = runTukeyPostHocMixed(wideTbl, measureVars, 'Day');
    writetable(tukeyTbl, fullfile(base_dir,[varName,'_PostHoc_Tukey.csv']));

    nLvls     = numel(measureVars);
    WithinDes = table((1:nLvls)', 'VariableNames', {'Day'});
    frm       = sprintf('%s-%s ~ Age', measureVars{1}, measureVars{end});
    rmModel   = fitrm(wideTbl, frm, 'WithinDesign', WithinDes);
    bfDay     = multcompare(rmModel, 'Day', 'By', 'Age', 'ComparisonType','bonferroni');

    postHocDay = array2table(bfDay, ...
        'VariableNames',{'Group','Level1','Level2','LowerCI','MeanDiff','UpperCI','pValue'});
    writetable(postHocDay, fullfile(base_dir,[varName,'_PostHoc_Bonferroni_Day.csv']));
end