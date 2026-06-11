% Code for Morris Water Maze
% This script runs all analysis for the Morris WaterMaze Code and Plotting 

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

% Find only unique barnesIds
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
allAccum = {table(), table(), table()};
missingRats = [];

% Loop through and fill tables (doing all three at once)

for i = 1:numel(files)
    fname = fullfile(base_dir, files(i).name);

    matchCoh = regexp(files(i).name,'Coh(\d+)_','tokens','once');
    if ~isempty(matchCoh)
        cohortNum = str2double(matchCoh{1});
    else
        cohortNum = NaN;
    end

    spT = readtable(fname,'Sheet','Spatial','VariableNamesRange',1,'VariableNamingRule','preserve');
    spCols = {'Test','Animal','Trial','Duration','Distance','PathEfficiency','Platform_CIPL'};
    for c = 1:numel(spCols)
        if ~ismember(spCols{c}, spT.Properties.VariableNames)
            spT.(spCols{c}) = nan(height(spT), 1);
        end
    end
    spT = spT(:, spCols);

    pbT = readtable(fname,'Sheet','Probe','VariableNamesRange',1,'VariableNamingRule','preserve');
    pbCols = {'Animal','Trial','Duration','Distance','Q1_Time','Platform_Entries', ...
        'Platform_PathEfficiencyToEntry','Platform_CIPL'};
    for c = 1:numel(pbCols)
        if ~ismember(pbCols{c}, pbT.Properties.VariableNames)
            pbT.(pbCols{c}) = nan(height(pbT), 1);
        end
    end
    pbT = pbT(:, pbCols);

    vsT = readtable(fname,'Sheet','Cued','VariableNamesRange',1,'VariableNamingRule','preserve');
    vsCols = {'Animal','Duration','Distance','Platform_CIPL'};
    for c = 1:numel(vsCols)
        if ~ismember(vsCols{c}, vsT.Properties.VariableNames)
            vsT.(vsCols{c}) = nan(height(vsT), 1);
        end
    end
    vsT = vsT(:, vsCols);

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

    for tIdx = 1:3
        newT = trialTables{tIdx};
        accT = allAccum{tIdx};
        if isempty(accT)
            allAccum{tIdx} = newT;
        else
            missingInNew = setdiff(accT.Properties.VariableNames, newT.Properties.VariableNames);
            for c = 1:numel(missingInNew)
                col = missingInNew{c};
                if ismember(col, {'Sex','Genotype'})
                    newT.(col) = strings(height(newT), 1);
                else
                    newT.(col) = nan(height(newT), 1);
                end
            end
            missingInAcc = setdiff(newT.Properties.VariableNames, accT.Properties.VariableNames);
            for c = 1:numel(missingInAcc)
                col = missingInAcc{c};
                if ismember(col, {'Sex','Genotype'})
                    accT.(col) = strings(height(accT), 1);
                else
                    accT.(col) = nan(height(accT), 1);
                end
            end
            newT = newT(:, accT.Properties.VariableNames);
            allAccum{tIdx} = [accT; newT];
        end
    end

    fprintf('\nMorris WaterMaze data added for cohort: %d', cohortNum);
end

AllSpatial = allAccum{1};
AllProbe   = allAccum{2};
AllVisual  = allAccum{3};

missingRats = unique(missingRats, 'stable');
if ~isempty(missingRats)
    missingTbl = table(missingRats, 'VariableNames', {'BarnesID'});
    writetable(missingTbl, fullfile(proc_dir, 'Tg_MissingFromAgingSheet.csv'));
    fprintf('\nWarning: %d Barnes ID(s) in spatial files were not found in the aging sheet.\n', numel(missingRats));
end

writetable(AllSpatial, fullfile(proc_dir,'Tg_AllMorrisWaterMazeData_Spatial.csv'));
writetable(AllProbe,   fullfile(proc_dir,'Tg_AllMorrisWaterMazeData_Probe.csv'));
writetable(AllVisual,  fullfile(proc_dir,'Tg_AllMorrisWaterMazeData_Visual.csv'));

%% Spatial Trials 
% This runs the analysis for calculation of ANOVA and Plotting by Day for
% the two age-groups for MWM

% For Spatial data
fname=fullfile(proc_dir,'Tg_AllMorrisWaterMazeData_Spatial.csv');
T = readtable(fname); % adjust path as needed
T = sortrows(T, ["Animal","Trial"]); 
% Categorize each rat as young or old where <15 is young and and >15 is old
% (Works for the distribution of our rats)
T.Group = repmat("Old",height(T),1);
T.Group(T.Age < 15) = "Young";
T.Day = floor((T.Trial-1)/6) + 1; % Separates trials by which day (6 trials per day)
%% First Trial COmparison
first_cipl=T.Platform_CIPL(T.Trial==1);
first_age=T.Group(T.Trial==1);
youngData=first_cipl(first_age=='Young');
oldData=first_cipl(first_age=='Old');
figure;
boxchart(ones(size(youngData)), youngData, ...
     'BoxFaceColor',[0 0.6 0], 'LineWidth',1.5); % green
hold on
boxchart(2*ones(size(oldData)), oldData, ...
     'BoxFaceColor',[0.5 0 0.5], 'LineWidth',1.5); % purple
% T-test (Welch's)
[~, p] = ttest2(youngData, oldData, 'Vartype','unequal');
fprintf('\n The P value for the first trial is %f \n',p);
title('First Trial CIPL Performance')
xticks([1 2])
pubify_figure_axis_robust(14,14)
xticklabels(['Young' 'Old'])
%% Rest of code
%Variables to interate over
varList = {'Platform_CIPL','Duration','Distance','PathEfficiency'};

%Iterate through variables
for v = 1:numel(varList)
    varName = varList{v}; 
    disp(['Processing variable: ', varName]);

    % ========== Group-level summary (Avg + Error Bars) ==========
    G = groupsummary(T, ["Group","Day"], {"mean","std","nnz"}, varName);
    RatDayAvg = groupsummary(T, ["Animal","Group","Day"], "mean", varName);

    % Figure for group averages
    fig1 = figure('Name',['Average ', varName],'Position',[95,100,800,630]);
    hold on
    grpList = ["Young","Old"];
    clrMap  = {[0 0.6 0],[0.5 0 0.5]}; % green, purple
    for gi = 1:numel(grpList)
        grp = grpList(gi);
        tmp = G(G.Group==grp,:);
        semVal = tmp.(['std_',varName]) ./ sqrt(tmp.(['nnz_',varName]));
        errorbar(tmp.Day, tmp.(['mean_',varName]), semVal, '-o', ...
                 'Color', clrMap{gi}, 'LineWidth',3, 'DisplayName', grp);
    end
    xlabel('Day','FontSize',12,'FontWeight','bold'); 
    ylabel(['Average ', varName],'FontSize',12,'FontWeight','bold','Interpreter','none'); 
    xticks([1 2 3 4])
    xlim([0.8 4.2])
    ymax = max(G.(['mean_',varName]) + 5);
    if ~isnan(ymax), ylim([0 ymax]); end
    legend('FontSize',12); 
    title(['Average ', varName,' Scores across Days (Spatial)'],'FontSize',14,'Interpreter','none');
    pubify_figure_axis_robust(16,16)
    hold off
    saveas(fig1, fullfile(proc_dir,['AvgMorrisWaterMaze_',varName]),'png');

    % ========== Plot All Animals ==========
    fig2 = figure('Name',['All ', varName],'Position',[95,100,800,630]);
    hold on
    ratsYoung = unique(RatDayAvg.Animal(RatDayAvg.Group=="Young"));
    ratsOld   = unique(RatDayAvg.Animal(RatDayAvg.Group=="Old"));

    % Plot faint lines for Young
    for r = 1:numel(ratsYoung)
        d = RatDayAvg(RatDayAvg.Animal==ratsYoung(r) & RatDayAvg.Group=="Young",:);
        d = sortrows(d,"Day");
        plot(d.Day, d.(['mean_',varName]), 'Color',[clrMap{1}, 0.175],'LineWidth',2);
    end
    % Plot faint lines for Old
    for r = 1:numel(ratsOld)
        d = RatDayAvg(RatDayAvg.Animal==ratsOld(r) & RatDayAvg.Group=="Old",:);
        d = sortrows(d,"Day");
        plot(d.Day, d.(['mean_',varName]), 'Color',[clrMap{2}, 0.175],'LineWidth',2);
    end

    % Bold group means on top
    y = plot(G.Day(G.Group=="Young"), G.(['mean_',varName])(G.Group=="Young"), ...
             '-o','Color',clrMap{1},'LineWidth',8,'DisplayName','Young');
    o = plot(G.Day(G.Group=="Old"),   G.(['mean_',varName])(G.Group=="Old"), ...
             '-o','Color',clrMap{2},'LineWidth',8,'DisplayName','Old');

    xlabel('Day','FontSize',14,'FontWeight','bold'); 
    %ylabel(['Average ', varName],'FontSize',12,'FontWeight','bold','Interpreter','none'); 
    ylabel(['CIPL Score (m\cdots)'],'FontSize',14,'FontWeight','bold','Interpreter','tex'); 
    xticks([1 2 3 4])
    xlim([0.95 4.05])
    ymax2 = max(RatDayAvg.(['mean_',varName]));
    if ~isnan(ymax2), ylim([0 ymax2]); end
    legend([y;o],{'Young','Old'},'FontSize',14,'FontWeight','bold');
    pubify_figure_axis_robust(16,16)
    title(['Performance on Morris watermaze'],'FontSize',20);
    
    hold off
    saveas(fig2, fullfile(proc_dir,['MorrisWaterMaze_',varName,'_All']),'png');

    
end

%% ========== MIXED‑DESIGN ANOVA + POST‑HOC FOR EACH VARIABLE ==========
for v = 1:numel(varList)
    varName = varList{v}; 
    disp(['Mixed ANOVA for variable: ', varName]);

    % 1) Build subject×day means table
    RatDayAvg = groupsummary(T, ["Animal","Group","Day"], "mean", varName);
    % rename the summary column to match varName
    RatDayAvg.Properties.VariableNames{end} = varName;
    % pivot so each row=1 animal, columns=Day1,Day2,...
    wideTbl = unstack(RatDayAvg, varName, "Day");
    % rename grouping vars to match your helper‑functions
    wideTbl.Properties.VariableNames{'Animal'} = 'Subject';
    wideTbl.Properties.VariableNames{'Group'}   = 'Age';

    % identify the repeated‑measure columns
    allCols    = wideTbl.Properties.VariableNames;
    measureVars = setdiff(allCols, {'Subject','Age'}, 'stable');

    % 2) run the mixed ANOVA
    anovaResults = runMixedANOVA(wideTbl, measureVars, 'Day');
    writetable(anovaResults, fullfile(base_dir, ...
        [varName,'_MixedANOVA_Results.csv']));

    % 3) run Tukey‑Kramer post‑hoc on Day within each Age
    tukeyTbl = runTukeyPostHocMixed(wideTbl, measureVars, 'Day');
    writetable(tukeyTbl, fullfile(base_dir, ...
        [varName,'_PostHoc_Tukey.csv']));

    % 4) Bonferroni‑corrected pairwise Day comparisons within each Age
    nLvls       = numel(measureVars);
    WithinDes   = table((1:nLvls)', 'VariableNames', {'Day'});
    frm         = sprintf('%s-%s ~ Age', measureVars{1}, measureVars{end});
    rmModel     = fitrm(wideTbl, frm, 'WithinDesign', WithinDes);
    bfDay       = multcompare(rmModel, 'Day', 'By', 'Age', ...
                              'ComparisonType','bonferroni');

    postHocDay = array2table(bfDay, ...
    'VariableNames', {'Group','Level1','Level2','LowerCI','MeanDiff','UpperCI','pValue'});

      writetable(postHocDay, fullfile(base_dir, ...
        [varName,'_PostHoc_Bonferroni_Day.csv']));
end






%% Probe trials 
%Get Probe trial data
fname=fullfile(proc_dir,'Tg_AllMorrisWaterMazeData_Probe.csv');
TProbe = readtable(fname); % adjust path as needed
% Categorize each rat as young or old where <15 is young and and >15 is old
TProbe.Platform_PathEfficiencyToEntry(isnan(TProbe.Platform_PathEfficiencyToEntry)) = 0;
TProbe.Group = repmat("Old",height(TProbe),1);
TProbe.Group(TProbe.Age < 15) = "Young";

%List of variables to iterate over
varList = {'Distance','Q1_Time','Platform_Entries','Platform_PathEfficiencyToEntry','Platform_CIPL'};

% Preallocate table to store p-values
pTable = table('Size',[numel(varList), 2], ...
               'VariableTypes', ["string","double"], ...
               'VariableNames', ["Variable","pValue"]);

% Iterate over variables
for v = 1:numel(varList)
    varName = varList{v};
    % Create figure
    fig = figure('Name', varName, 'Units','normalized','Position',[0.2 0.2 0.5 0.5]);

    % Separate Young vs Old
    youngData = TProbe.(varName)(TProbe.Group=="Young");
    oldData   = TProbe.(varName)(TProbe.Group=="Old");

    % Box chart
    boxchart(ones(size(youngData)), youngData, ...
             'BoxFaceColor',[0 0.6 0], 'LineWidth',1.5); % green
    hold on
    boxchart(2*ones(size(oldData)), oldData, ...
             'BoxFaceColor',[0.5 0 0.5], 'LineWidth',1.5); % purple
    hold off

    xlabel('Group'); 
    ylabel(varName, 'Interpreter','none');
    title(['Probe Trials: ', varName],'Interpreter','none','FontSize',14,'FontWeight','bold');
    pubify_figure_axis_robust(14,14)
    set(gca,'XTick',[1 2],'XTickLabel',{'Young','Old'},'FontSize',14);
   
    saveas(fig, fullfile(proc_dir, ['Probe_', varName, '.png']));

    % T-test (Welch's)
    [~, p] = ttest2(youngData, oldData, 'Vartype','unequal');
    pTable.Variable(v) = varName;
    pTable.pValue(v)   = p;
end

% Save all p-values in talbe
writetable(pTable, fullfile(base_dir, 'Probe_TTest_pValues.csv'));
disp(pTable)
