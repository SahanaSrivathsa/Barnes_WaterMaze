% Code for Morris Water Maze
% This script runs all analysis for the Morris WaterMaze Code and Plotting 
% Assumes all files as in the format 'Coh21_01.01.23_Spatial_W-Maze.xlsx
% Within each file the first sheet should be labelled Key and the Age
% column should contain the age of each rat followed by mo. 
% It also needs the BarnesID for animal number

%% GET DATA
% Parameters
base_dir='C:\DATA\WaterMaze\Mia_Proj\RawFiles';
proc_dir=fullfile(base_dir,'Processed');
% Mk dir 
if ~exist(proc_dir, "dir")
    mkdir(proc_dir);
end

% Get all .xlsx files in base_dir
files = dir(fullfile(base_dir,'*.xlsx'));
    
 % Preallocate tables
AllSpatial = table();
AllProbe   = table();
AllVisual  = table();

for i =1:numel(files)
    fname = fullfile(base_dir, files(i).name);

    %--- Parse cohort from filename: e.g. "Coh21_..."
    matchCoh = regexp(files(i).name,'Coh(\d+)_','tokens','once');
    if ~isempty(matchCoh)
        cohortNum = str2double(matchCoh{1});
    else
        cohortNum = NaN;
    end

    %=== Read "Key" sheet (for RatName and Age)
    keyT = readtable(fname,'Sheet','Key','VariableNamesRange',1);
    % Define Barnes ID - two different header names so using this field
    bnFields = intersect({'BarnesID','Barnes ID'}, keyT.Properties.VariableNames);
    if isempty(bnFields)
        warning("Barnes ID column not found in %s", fname);
        bnFields = {''};
    end
    
    %=== Read "Spatial" sheet
    spT = readtable(fname,'Sheet','Spatial','VariableNamesRange',1);
    spCols = {'Test','Animal','Trial','Duration','Distance','PathEfficiency','Platform_CIPL'};
    spKeep = intersect(spCols, spT.Properties.VariableNames);
    spT = spT(:,spKeep);
    %--- Add Age column if missing
    if ~ismember('Age', spT.Properties.VariableNames)
        spT.Age = nan(height(spT), 1);
    end
    %--- For each row in the Key table, fill the matching Animal's Age
    for k = 1:height(keyT)
        BN = keyT.(bnFields{1})(k);  % Barnes ID
        matchIdx = (spT.Animal == BN);% Find matching Animal rows in Spatial
        ageStr = string(keyT.Age(k));
        ageVal = str2double(regexp(ageStr,'\d+','match','once'));
        if any(matchIdx)
            spT.Age(matchIdx) = repmat(ageVal, sum(matchIdx), 1);
        end
    end
    spT.Cohort  = repmat(cohortNum,height(spT),1);
    AllSpatial  = [AllSpatial; spT];
    
    %=== Read "Probe" sheet
    pbT = readtable(fname,'Sheet','Probe','VariableNamesRange',1);
    pbCols = {'Animal','Trial','Duration','Distance','Q1_Time',...
        'Platform_Entries'};%,'Platform_CIPL' 0 Currenlty all dont have CIPL
    pbKeep = intersect(pbCols, pbT.Properties.VariableNames);
    pbT = pbT(:,pbKeep);%--- For each row in the Key table, fill the matching Animal's Age
    for k = 1:height(keyT)
        BN = keyT.(bnFields{1})(k);  % Barnes ID
        matchIdx = (pbT.Animal == BN);% Find matching Animal rows in Spatial
        ageStr = string(keyT.Age(k));
        ageVal = str2double(regexp(ageStr,'\d+','match','once'));
        if any(matchIdx)
            pbT.Age(matchIdx) = repmat(ageVal, sum(matchIdx), 1);
        end
    end
    pbT.Cohort  = repmat(cohortNum,height(pbT),1);
    %If you want to have some of the params from one file change them above
    %and uncomment the following code
    % missingCols = setdiff(AllProbe.Properties.VariableNames, pbT.Properties.VariableNames);
    % %  Create each missing column in pbT with NaN
    % for c = 1:numel(missingCols)
    %     pbT.(missingCols{c}) = NaN(height(pbT),1);
    % end
    % % Reorder pbT columns to match AllProbe
    % pbT = pbT(:, AllProbe.Properties.VariableNames);
    AllProbe    = [AllProbe; pbT];
    
    %=== Read "Visual" sheet
    vsT = readtable(fname,'Sheet','Visual','VariableNamesRange',1);
    vsCols = {'Animal','Duration','Distance''Platform_CIPL'};%,'PathEfficiency',
    vsKeep = intersect(vsCols, vsT.Properties.VariableNames);
    vsT = vsT(:,vsKeep);
    for k = 1:height(keyT)
        BN = keyT.(bnFields{1})(k);  % Barnes ID
        matchIdx = (vsT.Animal == BN);% Find matching Animal rows in Spatial
        ageStr = string(keyT.Age(k));
        ageVal = str2double(regexp(ageStr,'\d+','match','once'));
        if any(matchIdx)
            vsT.Age(matchIdx) = repmat(ageVal, sum(matchIdx), 1);
        end
    end
    vsT.Cohort  = repmat(cohortNum,height(vsT),1);
    AllVisual   = [AllVisual; vsT];

    fprintf('\n Morris WaterMaze Data Added for Cohort: %d',cohortNum);
end

% Write the final tables to CSV in the same folder
writetable(AllSpatial, fullfile(proc_dir,'AllMorrisWaterMazeData_Spatial.csv'));
writetable(AllProbe,   fullfile(proc_dir,'AllMorrisWaterMazeData_Probe.csv'));
writetable(AllVisual,  fullfile(proc_dir,'AllMorrisWaterMazeData_Visual.csv'));

%% Spatial Trials 
% This runs the analysis for calculation of ANOVA and Plotting by Day for
% the two age-groups for MWM

% For Spatial data
fname=fullfile(proc_dir,'AllMorrisWaterMazeData_Spatial.csv');
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
fname=fullfile(proc_dir,'AllMorrisWaterMazeData_Probe.csv');
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

