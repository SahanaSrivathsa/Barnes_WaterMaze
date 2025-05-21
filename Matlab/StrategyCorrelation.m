%% ------------- MUST RUN SECTION FOR ALL FURTHER ANALYSIS ------------- %%
base_dir='C:\DATA\WaterMaze\Mia_Proj';
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
data1 = readtable(fullfile(base_dir,'MWM_results_04-04-2025.xlsx'));  % From RTrack, contains Track_ID, Strategy, Age
data2 = readtable(fullfile(base_dir,'AllMorrisWaterMazeData_Spatial.csv'));  % From Matlab Analysis contains Test_No, Cohort, Platform_CIPL
% Analysis Parameters that are constant - mostly colours for groups
grpList = ["Young","Old"];
clrMap  = {[0.2196,0.5569,0.2353],[0.4157,0.1059,0.6039]}; % green, purple
% oldColor=[0.4157,0.1059,0.6039];
% youngColor=[0.2196,0.5569,0.2353];
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

% Get Platform Scores for CIPL - make sure the trial numbers align
platformScores = zeros(height(data2),1);

% Loop over each row in Strategy Data to find corresponding CIPL from
% second data set
for i = 1:height(data1)
    %Get Track ID which has cohort number and trial no.
    tID = data1.Track_ID{i};
    % Use regexp to extract numeric parts: cohort and test number
    tokens = regexp(tID, 'Coh(\d+)_test(\d+)', 'tokens');
    if ~isempty(tokens)
        token = tokens{1};
        cohortNum = str2double(token{1});
        testNum = str2double(token{2});
        % Find the matching row in data2 using cohort and test numbers
        idx = find(data2.Cohort == cohortNum & data2.Test == testNum, 1);
        idx2=find(data2.Animal==str2double(data1.x_TargetID{i}) &...
            data2.Trial == data1.x_Trial(i));
        if(idx~=idx2)
            platformScores(i)=NaN;
        elseif ~isempty(idx)
            platformScores(i) = data2.Platform_CIPL(idx);
        else
            platformScores(i)=NaN;
        end
    end
end
% Remove all the values where platform scores are nan
% Identify valid rows (non-NaN platform scores)
validIdx = ~isnan(platformScores);

% Filter data1 accordingly
data1 = data1(validIdx, :);
platformScores = platformScores(validIdx);  % also update platformScores to match

% Get corresponding animal-trial combinations from filtered data1
validAnimals = str2double(data1.x_TargetID);
validTrials = data1.x_Trial;

% Use these to find matching rows in data2
keepIdx = false(height(data2), 1);
for i = 1:length(validAnimals)
    match = data2.Animal == validAnimals(i) & data2.Trial == validTrials(i);
    if any(match)
        keepIdx = keepIdx | match;
    end
end

% Filter data2 to keep only matching animal-trial rows
data2 = data2(keepIdx, :);
%platformScores(platformScores<=0)=0;
%% 1) CIPL-Strategy  Plots- All trials
polyOrder = 2;  % Change as needed

for s = 1:nStrategies
    % Extract the probability vector for the current strategy
    stratProb = data1.(strategyNames{s});

    % Separate indices for young and old groups
    isYoung = strcmp(data1.Age, 'young');
    isOld   = strcmp(data1.Age, 'old');

    % All data (for scatter)
    x_young = platformScores(isYoung);
    y_young = stratProb(isYoung);
    x_old   = platformScores(isOld);
    y_old   = stratProb(isOld);

    % Filter out non-positive x values for polynomial fitting only
    x_young_fit = x_young(x_young > 0);
    y_young_fit = y_young(x_young > 0);
    x_old_fit   = x_old(x_old > 0);
    y_old_fit   = y_old(x_old > 0);

    f = figure; hold on;

    % Scatter plot (including zero values)
    scatter(x_young_fit, y_young_fit, 20, clrMap{1}, 'filled',...
        'MarkerFaceAlpha', 0.5, 'DisplayName', 'YOUNG');
    scatter(x_old_fit,   y_old_fit,   20, clrMap{2}, 'filled',...
        'MarkerFaceAlpha', 0.5, 'DisplayName', 'OLD');

    % % --- Polynomial Fit: Young group ---
    % if ~isempty(x_young_fit)
    %     pY = polyfit(x_young_fit, y_young_fit, polyOrder);
    %     xFit_y = linspace(min(x_young_fit), max(x_young_fit), 100);
    %     yFit_y = polyval(pY, xFit_y);
    %     plot(xFit_y, yFit_y, 'Color', clrMap{1}, 'LineWidth', 2, 'DisplayName', 'YOUNG');
    % end
    %
    % % --- Polynomial Fit: Old group ---
    % if ~isempty(x_old_fit)
    %     pO = polyfit(x_old_fit, y_old_fit, polyOrder);
    %     xFit_o = linspace(min(x_old_fit), max(x_old_fit), 100);
    %     yFit_o = polyval(pO, xFit_o);
    %     plot(xFit_o, yFit_o, 'Color', clrMap{2}, 'LineWidth', 3, 'DisplayName', 'OLD');
    % end

    title(sprintf('%s: CIPL', strategy_titles{s}));

    xlabel('CIPL Score (m.s)');
    ylabel('Probability of Strategy Use');
    legend('show');
    xlim([0 60]);
    ylim([0 1]);
    pubify_figure_axis_robust(14, 14);
    hold off;

    saveas(f, fullfile(fig_dir, 'CIPL_Strategy', sprintf('CIPL_%s', strategyNames{s})), 'png');
end

%% 2) CIPL vs group Probabilities
polyOrder = 2;  % CHANGE if needed based on fit - 1,-linear, 2- quadratic and 3- cubic.
strategyGroups = {
    {'thigmotaxis', 'circling', 'random_path'}, 'Non-Goal Oriented';
    {'scanning', 'chaining'}, 'Procedural';
    {'directed_search', 'corrected_search', 'direct_path','perseverance'}, 'Allocentric'
    };
groupNames = strategyGroups(:, 2);  % Extract group labels
for s = 1:numel(groupNames)
    % Get the strategies in the current group
    currentStrategies = strategyGroups{s, 1};
    groupProb=data1.(currentStrategies{1});
    for ii=2:numel(currentStrategies)
        groupProb=groupProb+data1.(currentStrategies{ii});
    end

    % Extract the probability vector for the current strategy
    stratProb = groupProb;

    % Separate indices for young and old groups based on the Age column
    isYoung = strcmp(data1.Age, 'young');
    isOld = strcmp(data1.Age, 'old');

    % Get data for each
    x_young = platformScores(isYoung);
    y_young = stratProb(isYoung);
    n_young = numel(x_young);

    % Data for old
    x_old = platformScores(isOld);
    y_old = stratProb(isOld);
    n_old = numel(x_old);

    f=figure;
    hold on;

    % Scatter plot (only scatter handles used in the legend)
    scatter(x_young, y_young, 20, clrMap{1}, 'filled',...
        'MarkerFaceAlpha',0.5, 'DisplayName', 'YOUNG');
    scatter(x_old, y_old, 20, clrMap{2}, 'filled',...
        'MarkerFaceAlpha',0.5, 'DisplayName', 'OLD');

    %  % NEW FIT
    % % Fit a polynomial of specified order
    %  p_young = polyfit(x_young, y_young, polyOrder);
    %  % Generate fitted curve
    %  xFit_y = linspace(min(x_young), max(x_young), 100);
    %  yFit_y = polyval(p_young, xFit_y);
    %  h1 =plot(xFit_y, yFit_y, 'Color', clrMap{1}, 'LineWidth', 2,'DisplayName','YOUNG');
    %
    %  % Compute R^2 for the young fit
    %  yHat_young = polyval(p_young, x_young);
    %  SSE_y = sum((y_young - yHat_young).^2);
    %  SST_y = sum((y_young - mean(y_young)).^2);
    %  R2_young = 1 - SSE_y / SST_y;
    %  r_young = sqrt(max(R2_young, 0));  % effective correlation from r^2
    %
    %  %Old fit
    %  p_old = polyfit(x_old, y_old, polyOrder);
    %  % Generate fitted curve
    %  xFit_o = linspace(min(x_old), max(x_old), 100);
    %  yFit_o = polyval(p_old, xFit_o);
    %  h2 =plot(xFit_o, yFit_o, 'Color', clrMap{2}, 'LineWidth', 3,'DisplayName','OLD');
    %
    %  % Compute r^2 for the old fit
    %  yHat_old = polyval(p_old, x_old);
    %  SSE_o = sum((y_old - yHat_old).^2);
    %  SST_o = sum((y_old - mean(y_old)).^2);
    %  R2_old = 1 - SSE_o / SST_o;
    %  r_old = sqrt(max(R2_old, 0));
    %
    %   % Fishers z transformation to compare correlation coefficients
    %  z_young = atanh(r_young);
    %  z_old = atanh(r_old);
    %  % For large-sample approximation, use n - polyOrder - 1
    %  se_diff = sqrt(1/(n_young - 1) + 1/(n_old - 1));
    %  z_stat  = (z_young - z_old) / se_diff;
    %  p_value = 2 * (1 - normcdf(abs(z_stat)));

    % Title with correlation coefficients
    % title(sprintf('%s: Young r=%.2f, Old r=%.2f p=%.3f', groupNames{s},...
    %     r_young, r_old,p_value),'FontSize',18,'FontWeight','bold');
    title(sprintf('%s', groupNames{s}),'FontSize',18,'FontWeight','bold');
    xlabel('CIPL Score (m.s)','FontSize',14,'FontWeight','bold');
    ylabel('Probability','FontSize',14,'FontWeight','bold');
    legend('show');
    xlim([0 60]);
    ylim([0 1]);
    pubify_figure_axis_robust(14,14)
    hold off;
    % Save figure

    saveas(f, fullfile(fig_dir, 'CIPL_Strategy', sprintf('CIPL_%s',groupNames{s})),'png');
end

%% 3) Strategy-Day Plots (individual and group)
apaTbl = table( ...
    strings(0,1) , strings(0,1) , ...   % Strategy  , Effect
    zeros (0,1)  , ...                  % df
    zeros (0,1)  , zeros (0,1) , ...    % SS , MS
    zeros (0,1)  , zeros (0,1)  , ...   % F  , p
    zeros (0,1)  , ...                  % eta2
    'VariableNames', ...
    {'Strategy','Effect','df','SS','MS','F','p','eta2'});


for s = 1:nStrategies

    %-------------- Individual Rat Plot --------------%
    stratProb = data1.(strategyNames{s});
    % Get idx for young and old groups based on the age
    isYoung = strcmp(data1.Age, 'young');
    isOld   = strcmp(data1.Age, 'old');

    % Unique rats and days
    uniqueDays = unique(data1.Day);
    uniqueRats = unique(data1.x_TargetID);

    %Init fig
    f1 = figure;
    hold on;

    % Initialize cell array to save mean strategy
    mean_strat = {};
    validRatIdx = 0;
    firstYoung = true;
    firstOld   = true;

    % Loop over each rat
    for r = 1:numel(uniqueRats)
        ratID = uniqueRats{r};
        isCurrentRat = strcmp(data1.x_TargetID, ratID);

        % Determine if the rat is young or old and choose color
        if ismember(ratID, data1.x_TargetID(isYoung))
            ageGroup = 'young';
            color = clrMap{1};
        else
            ageGroup = 'old';
            color = clrMap{2};
        end

        % Preallocate mean use for each day
        meanUse = nan(1, numel(uniqueDays));
        % Loop over each day
        for d = 1:numel(uniqueDays)
            idx = isCurrentRat & (data1.Day == uniqueDays(d));
            if any(idx)
                meanUse(d) = mean(stratProb(idx));
            end
        end

        % Skip rats with missing data for any day -sanity check
        if any(isnan(meanUse))
            fprintf('Missing days for %s in strategy %s\n', string(ratID), strategyNames{s});
            continue;
        end

        % Save the datafor ANOVA
        validRatIdx = validRatIdx + 1;
        mean_strat{validRatIdx, 1} = ratID;
        mean_strat{validRatIdx, 2} = ageGroup;
        mean_strat(validRatIdx, 3:6) = num2cell(meanUse); % Day1-Day4

        % Plot individual points with jitter using swarmchart
        hSwarm = swarmchart(uniqueDays, meanUse, 20, 'filled', ...
            'MarkerFaceColor', color, 'MarkerFaceAlpha', 0.5);
        hSwarm.Annotation.LegendInformation.IconDisplayStyle = 'off';

        % Plot connecting lines -display name only for the first rat of
        % each group
        if strcmp(ageGroup, 'young')
            if firstYoung
                plot(uniqueDays, meanUse, '-', 'Color', color, 'LineWidth', 1, 'DisplayName', 'Young');
                firstYoung = false;
            else
                plot(uniqueDays, meanUse, '-', 'Color', color, 'LineWidth', 1, 'HandleVisibility','off');
            end
        else
            if firstOld
                plot(uniqueDays, meanUse, '-', 'Color', color, 'LineWidth', 1, 'DisplayName', 'Old');
                firstOld = false;
            else
                plot(uniqueDays, meanUse, '-', 'Color', color, 'LineWidth', 1, 'HandleVisibility','off');
            end
        end
    end

    % Figure details for individual rat plot
    title(strategy_titles{s});
    xlabel('Day');
    ylabel('Probability of Strategy Use');
    xticks(uniqueDays);
    legend('Location', 'Northeast');
    ylim([0 0.5]);
    xlim([0.75 4.25]);
    pubify_figure_axis_robust(14,14);
    hold off;
    % Save fig
    saveas(f1, fullfile(fig_dir,'StrategyUse', sprintf('IndividualRats_%s', strategyNames{s})), 'png');
    close;

    % Convert cell array to table for ANOVA and post hoc tests

    mean_strat_table = cell2table(mean_strat, 'VariableNames', ...
        {'RatID', 'Age', 'Day1', 'Day2', 'Day3', 'Day4'});
    % Run mixed-design ANOVA (within-subject-Day &between-subject-age)
    % fprintf('%s   N = %d (young=%d, old=%d)\n', ...
    %     strategyNames{s}, ...
    %     height(mean_strat_table), ...
    %     sum(strcmp(mean_strat_table.Age,'young')), ...
    %     sum(strcmp(mean_strat_table.Age,'old')));

    anovaResults = runMixedANOVA(mean_strat_table, {'Day1','Day2','Day3','Day4'});
    % Run Tukey post hoc tests for the Day factor (by Age)
    postHocResults = runTukeyPostHocMixed(mean_strat_table, {'Day1','Day2','Day3','Day4'});

    % Save ANOVA results
    writetable(anovaResults, fullfile(processed_dir, sprintf('Anova_%s.csv', strategyNames{s})));
    writetable(postHocResults, fullfile(processed_dir, sprintf('PostHoc_Tukey_%s.csv', strategyNames{s})));


    rowNames = string(anovaResults.Properties.RowNames);

    % ----------  APA-style rows for this strategy ----------
    if ~strcmpi(strategyNames{s},'perseverance')

        rowNames = string(anovaResults.Properties.RowNames);

        % Map true row names → tidy labels
        wanted   = ["Age","(Intercept):Day","Age:Day"];     % what ranova returns
        tidyLab  = ["Age","Day","Age×Day"];                 % how you want them shown

        keepIdx  = find(ismember(rowNames, wanted));
        if isempty(keepIdx),  warning('%s has no mixed-ANOVA rows',strategyNames{s});  return; end
        L        = numel(keepIdx);                          % 2 or 3 rows

        % Pull columns once
        SS   = anovaResults.SumSq(keepIdx);
        df   = anovaResults.DF   (keepIdx);
        MS   = anovaResults.MeanSq(keepIdx);
        Fval = anovaResults.F    (keepIdx);
        pVal = anovaResults.pValue(keepIdx);

        % η²:  effect SS / (effect SS + error SS)  --------------------------
        eta2 = nan(L,1);
        for j = 1:L
            k    = keepIdx(j);                              % row of effect
            eRow = find(startsWith(rowNames(k+1:end),"Error"),1,'first') + k;
            if ~isempty(eRow)
                eta2(j) = SS(j) / ( SS(j) + anovaResults.SumSq(eRow) );
            end
        end

        % Build little table (every column = L×1)  --------------------------
        tmp = table( ...
            repmat(string(strategy_titles{s}), L,1)  , ...   % Strategy
            tidyLab(ismember(wanted,rowNames(keepIdx))).' , ... % Effect
            df   , SS , MS , Fval , pVal , eta2 , ...
            'VariableNames',{'Strategy','Effect','df','SS','MS','F','p','eta2'});

        apaTbl = [apaTbl ; tmp];                             % append
    end

    % -------------------------------------------------------------------------

    %-------------- Overall Bar Plot with Mean, SEM, and Significance Markers --------------%
    % Compute means and SEM for young and old per day
    meanYoung = arrayfun(@(d) mean(stratProb(isYoung & (data1.Day == d))), uniqueDays);
    semYoung  = arrayfun(@(d) std(stratProb(isYoung & (data1.Day == d))) / sqrt(sum(isYoung & (data1.Day == d))), uniqueDays);
    meanOld   = arrayfun(@(d) mean(stratProb(isOld & (data1.Day == d))), uniqueDays);
    semOld    = arrayfun(@(d) std(stratProb(isOld & (data1.Day == d))) / sqrt(sum(isOld & (data1.Day == d))), uniqueDays);

    % Extract individual data points for jitter plotting (cell arrays, one cell per day)
    dataYoung = arrayfun(@(d) stratProb(isYoung & (data1.Day == d)), uniqueDays, 'UniformOutput', false);
    dataOld   = arrayfun(@(d) stratProb(isOld & (data1.Day == d)), uniqueDays, 'UniformOutput', false);

    % % Create bar plot using the separate function plot_bar_sem
    % f2 = figure;
    % plot_bar_sem_WaterMaze(uniqueDays, meanYoung, semYoung, meanOld, semOld,...
    %     dataYoung, dataOld, postHocResults);
    % %PlotParams
    % title(sprintf('%s: All Trials',strategy_titles{s}),'FontSize',16);
    %
    % xlabel('Day');
    % ylabel('Probability of Strategy Use');
    % legend({'Young', 'Old'}, 'Location', 'Northeast');
    % ylim([0 0.75]);
    % pubify_figure_axis_robust(14,14);
    % hold off;
    %
    %
    % % Save the overall bar plot
    % saveas(f2, fullfile(fig_dir, 'StrategyUse',sprintf('ALLTRIALS_%s', strategyNames{s})), 'png');
    % close;
    %Figure for mean of each trial per day
    dataMeanYoung = cell(numel(uniqueDays), 1);
    dataMeanOld   = cell(numel(uniqueDays), 1);
    for d = 1:numel(uniqueDays)
        colName = sprintf('Day%d', d);
        dataMeanYoung{d} = mean_strat_table{strcmp(mean_strat_table.Age, 'young'), colName};
        dataMeanOld{d}   = mean_strat_table{strcmp(mean_strat_table.Age, 'old'), colName};
    end
    meanYoung = cellfun(@mean, dataMeanYoung);
    semYoung  = cellfun(@(x) std(x)/sqrt(numel(x)), dataMeanYoung);
    meanOld   = cellfun(@mean, dataMeanOld);
    semOld    = cellfun(@(x) std(x)/sqrt(numel(x)), dataMeanOld);

    f3=figure;
    plot_bar_sem_WaterMaze(uniqueDays, meanYoung, semYoung, meanOld, semOld, dataMeanYoung,...
        dataMeanOld, postHocResults);
    %PlotParams
    % title(sprintf('%s: Rat Mean Per Day',strategy_titles{s}),'FontSize',16);
    %xlabel('Day');
    %ylabel('Probability of Strategy Use');
    %legend({'Young', 'Old'}, 'Location', 'Northeast');
    ylim([0 0.55]);
    pubify_figure_axis_robust(16,16);
    hold off;

    % Save the overall bar plot
    exportgraphics(f3, fullfile(fig_dir, 'StrategyUse',...
        sprintf('1_RatMeanPerDay_%s.png', strategyNames{s})),'Resolution', 450 );
    close
end

% -------- build one tidy APA table ---------
writetable(apaTbl, fullfile(processed_dir,'ANOVA_APA_AllStrategies.csv'));
disp(apaTbl)
saveTablePNG_APA(apaTbl, ...
    fullfile(fig_dir,'ANOVA_APA_AllStrategies.png'), ...
    'Title','Repeated Measures ANOVA for Strategies');
%% 4) Strategy by Group - Non-spatial/Procedural/Allocentric
% Define strategy groups and their names
strategyGroups = {
    {'thigmotaxis', 'circling', 'random_path'}, 'Non-Goal Oriented';
    {'scanning', 'chaining'}, 'Procedural';
    {'directed_search', 'corrected_search', 'direct_path','perseverance'}, 'Allocentric'
    };
groupNames = strategyGroups(:, 2);  % Extract group labels
nGroups = numel(groupNames);        % Number of groups

groupStrat=cell(1,3);
% Extract unique days
uniqueDays = unique(data1.Day);
%---------------------- GET MEAN and SEM per group ----------------------%
% Loop over each strategy group
for g = 1:nGroups
    % Initialize figure for individual rat plots
    f1 = figure;
    hold on;

    % Get the strategies in the current group
    currentStrategies = strategyGroups{g, 1};
    groupProb=data1.(currentStrategies{1});
    for ii=2:numel(currentStrategies)
        groupProb=groupProb+data1.(currentStrategies{ii});
    end
    groupStrat{g}=groupProb;

    % Get idx for young and old groups based on the age
    isYoung = strcmp(data1.Age, 'young');
    isOld   = strcmp(data1.Age, 'old');

    % Initialize cell array to save mean strategy use per rat
    mean_strat = {};
    validRatIdx = 0;
    firstYoung = true;
    firstOld = true;

    % Loop over each rat
    uniqueRats = unique(data1.x_TargetID);
    for r = 1:numel(uniqueRats)
        ratID = uniqueRats{r};
        isCurrentRat = strcmp(data1.x_TargetID, ratID);

        % Determine if the rat is young or old and choose color
        if ismember(ratID, data1.x_TargetID(strcmp(data1.Age, 'young')))
            ageGroup = 'young';
            color = clrMap{1};
        else
            ageGroup = 'old';
            color = clrMap{2};
        end

        % Preallocate mean use for each day
        meanUse = nan(1, numel(uniqueDays));

        % Loop over each day
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

        % Save the data for ANOVA
        validRatIdx = validRatIdx + 1;
        mean_strat{validRatIdx, 1} = ratID;
        mean_strat{validRatIdx, 2} = ageGroup;
        mean_strat(validRatIdx, 3:6) = num2cell(meanUse); % Day1-Day4

        % Plot individual points with jitter using swarmchart
        hSwarm = swarmchart(uniqueDays, meanUse, 20, 'filled', ...
            'MarkerFaceColor', color, 'MarkerFaceAlpha', 0.5);
        hSwarm.Annotation.LegendInformation.IconDisplayStyle = 'off';

        % Plot connecting lines - display name only for the first rat of each group
        if strcmp(ageGroup, 'young')
            if firstYoung
                plot(uniqueDays, meanUse, '-', 'Color', color, 'LineWidth', 1, 'DisplayName', 'Young');
                firstYoung = false;
            else
                plot(uniqueDays, meanUse, '-', 'Color', color, 'LineWidth', 1, 'HandleVisibility', 'off');
            end
        else
            if firstOld
                plot(uniqueDays, meanUse, '-', 'Color', color, 'LineWidth', 1, 'DisplayName', 'Old');
                firstOld = false;
            else
                plot(uniqueDays, meanUse, '-', 'Color', color, 'LineWidth', 1, 'HandleVisibility', 'off');
            end
        end
    end

    % Figure details for individual rat plot
    title(groupNames{g});
    xlabel('Day');
    ylabel('Probability of Strategy Use');
    xticks(uniqueDays);
    legend('Location', 'northeastoutside');
    ylim([0 1]);
    xlim([0.75 4.25]);
    pubify_figure_axis_robust(14, 14);
    hold off;

    % Save individual rat plot
    saveas(f1, fullfile(fig_dir, 'StrategyUse',sprintf('IndividualRats_%s', groupNames{g})), 'png');
    close;

    % Convert cell array to table for ANOVA and post hoc tests
    mean_strat_table = cell2table(mean_strat, 'VariableNames', ...
        {'RatID', 'Age', 'Day1', 'Day2', 'Day3', 'Day4'});

    % Run mixed-design ANOVA (within-subject-Day & between-subject-age)
    anovaResults = runMixedANOVA(mean_strat_table, {'Day1', 'Day2', 'Day3', 'Day4'});
    % Run Tukey post hoc tests for the Day factor (by Age)
    postHocResults = runTukeyPostHocMixed(mean_strat_table, {'Day1', 'Day2', 'Day3', 'Day4'});

    % Save ANOVA results
    writetable(anovaResults, fullfile(processed_dir, sprintf('Anova_%s.csv', groupNames{g})));
    writetable(postHocResults, fullfile(processed_dir, sprintf('PostHoc_Tukey_%s.csv', groupNames{g})));

    %-------------- Overall Bar Plot with Mean, SEM, and Significance Markers --------------%

    % Compute means and SEM for young and old per day using all trials
    meanYoung = arrayfun(@(d) mean(groupProb(isYoung & (data1.Day == d))), uniqueDays);
    semYoung  = arrayfun(@(d) std(groupProb(isYoung & (data1.Day == d))) / sqrt(sum(isYoung & (data1.Day == d))), uniqueDays);
    meanOld   = arrayfun(@(d) mean(groupProb(isOld & (data1.Day == d))), uniqueDays);
    semOld    = arrayfun(@(d) std(groupProb(isOld & (data1.Day == d))) / sqrt(sum(isOld & (data1.Day == d))), uniqueDays);

    % Extract individual data points for jitter plotting (cell arrays, one cell per day)
    dataYoung = arrayfun(@(d) groupProb(isYoung & (data1.Day == d)), uniqueDays, 'UniformOutput', false);
    dataOld   = arrayfun(@(d) groupProb(isOld & (data1.Day == d)), uniqueDays, 'UniformOutput', false);

    % % Create bar plot using the separate function plot_bar_sem
    % f2 = figure;
    % plot_bar_sem_WaterMaze(uniqueDays, meanYoung, semYoung, meanOld, semOld,...
    %     dataYoung, dataOld, postHocResults);
    % %PlotParams
    % title(sprintf('%s: All Trials',groupNames{g}),'FontSize',16);
    %
    % xlabel('Day');
    % ylabel('Probability of Strategy Use');
    % %legend({'Young', 'Old'}, 'Location', 'northeastoutside');
    % ylim([0 1.35]);
    % pubify_figure_axis_robust(14,14);
    % hold off;
    %
    %  % Save fig
    % saveas(f2, fullfile(fig_dir, sprintf('%s_AllTrials', groupNames{g})), 'png');
    % close;

    %-------------- Figure for Mean Per Rat Per Day --------------%

    dataMeanYoung = cell(numel(uniqueDays), 1);
    dataMeanOld   = cell(numel(uniqueDays), 1);

    for d = 1:numel(uniqueDays)
        colName = sprintf('Day%d', d);
        dataMeanYoung{d} = mean_strat_table{strcmp(mean_strat_table.Age, 'young'), colName};
        dataMeanOld{d}   = mean_strat_table{strcmp(mean_strat_table.Age, 'old'), colName};
    end

    meanYoung = cellfun(@mean, dataMeanYoung);
    semYoung  = cellfun(@(x) std(x) / sqrt(numel(x)), dataMeanYoung);
    meanOld   = cellfun(@mean, dataMeanOld);
    semOld    = cellfun(@(x) std(x) / sqrt(numel(x)), dataMeanOld);

    f3 = figure;
    plot_bar_sem_WaterMaze(uniqueDays, meanYoung, semYoung, meanOld, semOld, dataMeanYoung, dataMeanOld, postHocResults);

    % Plot parameters
    title(sprintf('%s: Rat Mean Per Day', groupNames{g}), 'FontSize', 16);
    xlabel('Day');
    ylabel('Probability of Strategy Use');
    %legend({'Young', 'Old'}, 'Location', 'northeastoutside');
    ylim([0 1.25]);
    pubify_figure_axis_robust(14, 14);
    hold off;

    % Save the overall bar plot for mean per rat per day
    saveas(f3, fullfile(fig_dir, 'StrategyUse',sprintf('RatMeanPerDay_%s', groupNames{g})), 'png');
    close;
end





%% 5) Entropy Calculation
 oldColor=[0.4157,0.1059,0.6039]; %green
 youngColor=[0.2196,0.5569,0.2353];% purple

%------------------ Compute Entropy for Each Trial ------------------%
% Entropy formula: H = -sum(p .* log2(p + eps))
pMat = data1{:, 11:19};
entropy_vals = -sum(pMat .* log2(pMat + eps), 2);
data1.entropy = entropy_vals;

% Make Entropy Table for stats
uniqueRats = unique(data1.x_TargetID);
uniqueDays = unique(data1.Day);
nRats = numel(uniqueRats);

% Preallocate cell arrays
meanDayEntropy = nan(nRats, numel(uniqueDays));
AgeCell = cell(nRats,1);

for r = 1:nRats
    ratID = uniqueRats{r};
    idxRat = strcmp(data1.x_TargetID, ratID);
    % Assume age is consistent for a rat.
    AgeCell{r} = data1.Age{find(idxRat,1)};
    for d = 1:numel(uniqueDays)
        idxDay = idxRat & (data1.Day == uniqueDays(d));
        meanDayEntropy(r,d) = mean(data1.entropy(idxDay));
    end
end

%Generate Table
entropyTable = table(uniqueRats, AgeCell, meanDayEntropy(:,1), meanDayEntropy(:,2), meanDayEntropy(:,3), meanDayEntropy(:,4), ...
    'VariableNames', {'RatID', 'Age', 'Day1', 'Day2', 'Day3', 'Day4'});


% Run mixed-design ANOVA (within-subject-Day &between-subject-age)
anovaResults = runMixedANOVA(entropyTable, {'Day1','Day2','Day3','Day4'});
postHocResults = runTukeyPostHocMixed(entropyTable, {'Day1','Day2','Day3','Day4'});
writetable(anovaResults, fullfile(processed_dir, 'Entropy_ANOVA.csv'));
writetable(postHocResults, fullfile(processed_dir, 'Entropy_PostHoc_Tukey.csv'));


%------------------Entropy for Each Trial for Each Rat ------------------%
fe1=figure;
hold on;
uniqueRats = unique(data1.x_TargetID);  % Unique rat IDs
jitterAmount = 0.1;  % Adjust the jitter range as needed

for r = 1:length(uniqueRats)
    ratID = uniqueRats{r};
    idx = strcmp(data1.x_TargetID, ratID);

    % Get Trials
    trialNumbers = data1.x_Trial(idx);

    % Add jitter to the trial numbers
    jitteredTrials_Y = trialNumbers - (rand(size(trialNumbers))-0.25)*jitterAmount -0.1;
    jitteredTrials_O = trialNumbers + (rand(size(trialNumbers))+0.25)*jitterAmount +0.1;

    % Extract entropy values for this rat
    entropy_r = data1.entropy(idx);

    % Determine age
    if strcmp(data1.Age{find(idx,1)}, 'young')
        scatter(jitteredTrials_Y, entropy_r, 36, 'MarkerFaceColor', youngColor, 'MarkerEdgeColor', 'none');
    else
        scatter(jitteredTrials_O, entropy_r, 36, 'MarkerFaceColor', oldColor, 'MarkerEdgeColor', 'none');
    end
end
% Figure Properties
xlabel('Trial');
ylabel('Entropy');
title('Entropy: All Trials and All Rats');
pubify_figure_axis_robust(14,14);
hold off;
saveas(fe1, fullfile(fig_dir,  'Entropy','EntropyAllTrials'), 'png');
close;



%------------------ Calculate Mean Entropy per Rat per Day ------------------%
uniqueAges = {'young','old'};


% For each day and each age group, compute the group mean and SEM from the per-rat means.
meanEntropy = zeros(numel(uniqueDays), numel(uniqueAges));
semEntropy = zeros(numel(uniqueDays), numel(uniqueAges));
% Also build cell arrays holding the per-rat means for each day (for scatter overlay)
dataEntropyYoung = cell(numel(uniqueDays),1);
dataEntropyOld   = cell(numel(uniqueDays),1);

for d = 1:numel(uniqueDays)
    for j = 1:numel(uniqueAges)
        idx = strcmp(entropyTable.Age, uniqueAges{j});
        % Extract the per-rat mean for day d
        vals = meanDayEntropy(idx, d);
        meanEntropy(d,j) = mean(vals,'omitnan');
        semEntropy(d,j)  = std(vals, 'omitnan') / sqrt(sum(~isnan(vals)));
        % Save the per-rat means into cell arrays
        if strcmp(uniqueAges{j}, 'young')
            dataEntropyYoung{d} = vals;
        elseif strcmp(uniqueAges{j}, 'old')
            dataEntropyOld{d} = vals;
        end
    end
end


%---------------- PLOT PER RAT PER DAY MEAN ENTROPY ---------------------------%
fe2=figure;
hold on;
plot_bar_sem_WaterMaze(uniqueDays, meanEntropy(:,1), semEntropy(:,1), ...
    meanEntropy(:,2), semEntropy(:,2), dataEntropyYoung, dataEntropyOld, postHocResults);
title('Mean Entropy ');
xlabel('Day');
ylabel('Entropy');
pubify_figure_axis_robust(14,14);
xticks(uniqueDays);
hold off;
saveas(fe2, fullfile(fig_dir, 'Entropy', 'Entropy_MeanRat'), 'png');
close;
% Plot all rats as line plot
% ---------------- PLOT PER-RAT TRAJECTORIES + GROUP MEAN ---------------- %
fe3 = figure;  clf
hold on;
 oldColor=[0.4157,0.1059,0.6039]; %green
 youngColor=[0.2196,0.5569,0.2353];% purple

%  Plot each rat with light transparency %
alphaIndiv = 0.2;          % transparency for individual lines
lwIndiv    = 1.5;           % thin
lwMean     = 6;           % thick mean line

for r = 1:nRats
    if strcmpi(AgeCell{r}, 'young')
        thisCol = [youngColor alphaIndiv];  
    else
        thisCol = [oldColor   alphaIndiv];
    end
    plot(uniqueDays, meanDayEntropy(r,:), '-o', ...
         'Color', thisCol,  ...
         'LineWidth', lwIndiv, ...
         'MarkerSize', 4, ...
         'MarkerFaceColor', thisCol(1:3), ...
         'MarkerEdgeColor', 'none');
end

% --- 2. Compute & draw group mean lines  %
for j = 1:numel(uniqueAges)
    ageTag = uniqueAges{j};
    ageMask = strcmpi(AgeCell, ageTag);
    grpMean = mean(meanDayEntropy(ageMask, :), 1, 'omitnan');
    if strcmpi(ageTag,'young')
        plot(uniqueDays, grpMean, '-o', ...
            'Color', youngColor, 'LineWidth', lwMean, ...
            'MarkerSize', 6, 'MarkerFaceColor', youngColor);
    else
        plot(uniqueDays, grpMean, '-o', ...
            'Color', oldColor,   'LineWidth', lwMean, ...
            'MarkerSize', 6, 'MarkerFaceColor', oldColor);
    end
end

xlabel('Day');
ylabel('Entropy');
title('Entropy Across All 8 Strategy Types');
xticks(uniqueDays);
pubify_figure_axis_robust(14,14);
hold off;

saveas(fe3, fullfile(fig_dir, 'Entropy', 'Entropy_All8_MeanRat'), 'png');
close(fe3);

% fe3 = figure;
% hold on;
% nRats = size(dayEntropy,1);  % Number of rats
% for r = 1:nRats
%     if strcmp(AgeCell{r}, 'young')
%         plot(uniqueDays, dayEntropy(r,:), '-o', 'Color', youngColor, 'LineWidth', 1.5);
%     else
%         plot(uniqueDays, dayEntropy(r,:), '-o', 'Color', oldColor, 'LineWidth', 1.5);
%     end
% end
% xlabel('Day');
% ylabel('Entropy');
% title('Entropy Per Rat');
% pubify_figure_axis_robust(14,14);
% xticks(uniqueDays);
% hold off;
% saveas(fe3, fullfile(fig_dir, 'Entropy_RatChanges'), 'png');
% close;

%--------------- PLOT PER RAT PER DAY MEAN ENTROPY ----------------------%
% fe4=figure;
% hold on;
% % Calculate dataEntropy for all trials
% dataEntropyYoung = cell(numel(uniqueDays),1);
% dataEntropyOld   = cell(numel(uniqueDays),1);
% for d = 1:numel(uniqueDays)
%     dataEntropyYoung{d}=entropy_vals(uniqueDays(d) & strcmp(data1.Age,'young'));
%     dataEntropyOld{d}=entropy_vals(data1.Day == uniqueDays(d) & strcmp(data1.Age,'old'));
% end
% plot_bar_sem_WaterMaze(uniqueDays, meanEntropy(:,1), semEntropy(:,1), ...
%     meanEntropy(:,2), semEntropy(:,2), dataEntropyYoung, dataEntropyOld, postHocResults);
% title('Mean Entropy ');
% xlabel('Day');
% ylabel('Entropy');
% pubify_figure_axis_robust(14,14);
% xticks(uniqueDays)
% hold off;
% saveas(fe4, fullfile(fig_dir, 'Entropy','Entropy_Mean_AllTrials'), 'png');
% close;
%% 6) Group Strategy Based Entropy
% NEED TO HAVE RUN GROUP STRAT BEFORE
%------------------ Compute Entropy for Each Trial ------------------%
% Entropy formula: H = -sum(p .* log2(p + eps))
pMat = [groupStrat{1}';groupStrat{2}';groupStrat{3}']';
entropy_vals = -sum(pMat .* log2(pMat + eps), 2);
data1.entropy = entropy_vals;

% Make Entropy Table for stats
uniqueRats = unique(data1.x_TargetID);
uniqueDays = unique(data1.Day);
nRats = numel(uniqueRats);

% Preallocate cell arrays
meanDayEntropy = nan(nRats, numel(uniqueDays));
AgeCell = cell(nRats,1);

for r = 1:nRats
    ratID = uniqueRats{r};
    idxRat = strcmp(data1.x_TargetID, ratID);
    % Assume age is consistent for a rat.
    AgeCell{r} = data1.Age{find(idxRat,1)};
    for d = 1:numel(uniqueDays)
        idxDay = idxRat & (data1.Day == uniqueDays(d));
        meanDayEntropy(r,d) = mean(data1.entropy(idxDay));
    end
end

%Generate Table
entropyTable = table(uniqueRats, AgeCell, meanDayEntropy(:,1), meanDayEntropy(:,2), meanDayEntropy(:,3), meanDayEntropy(:,4), ...
    'VariableNames', {'RatID', 'Age', 'Day1', 'Day2', 'Day3', 'Day4'});


% Run mixed-design ANOVA (within-subject-Day &between-subject-age)
anovaResults = runMixedANOVA(entropyTable, {'Day1','Day2','Day3','Day4'});
disp(anovaResults)
postHocResults = runTukeyPostHocMixed(entropyTable, {'Day1','Day2','Day3','Day4'});
writetable(anovaResults, fullfile(processed_dir, 'GroupStrat_Entropy_ANOVA.csv'));
writetable(postHocResults, fullfile(processed_dir, 'GroupStrat_Entropy_PostHoc_Tukey.csv'));


%------------------Entropy for Each Trial for Each Rat ------------------%
fe1=figure;
hold on;
uniqueRats = unique(data1.x_TargetID);  % Unique rat IDs
jitterAmount = 0.1;  % Adjust the jitter range as needed

for r = 1:length(uniqueRats)
    ratID = uniqueRats{r};
    idx = strcmp(data1.x_TargetID, ratID);

    % Get Trials
    trialNumbers = data1.x_Trial(idx);

    % Add jitter to the trial numbers
    jitteredTrials_Y = trialNumbers - (rand(size(trialNumbers))-0.25)*jitterAmount -0.1;
    jitteredTrials_O = trialNumbers + (rand(size(trialNumbers))+0.25)*jitterAmount +0.1;

    % Extract entropy values for this rat
    entropy_r = data1.entropy(idx);

    % Determine age
    if strcmp(data1.Age{find(idx,1)}, 'young')
        scatter(jitteredTrials_Y, entropy_r, 36, 'MarkerFaceColor', youngColor, 'MarkerEdgeColor', 'none');
    else
        scatter(jitteredTrials_O, entropy_r, 36, 'MarkerFaceColor', oldColor, 'MarkerEdgeColor', 'none');
    end
end
% Figure Properties
xlabel('Trial');
ylabel('Entropy');
title('Entropy: All Trials  - Grouped By Strategy Type');
pubify_figure_axis_robust(14,14);
hold off;
saveas(fe1, fullfile(fig_dir, 'Entropy', 'GroupStrat_EntropyAllTrials'), 'png');
close;



%------------------ Calculate Mean Entropy per Rat per Day ------------------%
uniqueAges = {'young','old'};


% For each day and each age group, compute the group mean and SEM from the per-rat means.
meanEntropy = zeros(numel(uniqueDays), numel(uniqueAges));
semEntropy = zeros(numel(uniqueDays),numel(uniqueAges));
% Also build cell arrays holding the per-rat means for each day (for scatter overlay)
dataEntropyYoung = cell(numel(uniqueDays),1);
dataEntropyOld   = cell(numel(uniqueDays),1);

for d = 1:numel(uniqueDays)
    for j = 1:numel(uniqueAges)
        idx = strcmp(entropyTable.Age, uniqueAges{j});
        % Extract the per-rat mean for day d
        vals = meanDayEntropy(idx, d);
        meanEntropy(d,j) = mean(vals,'omitnan');
        semEntropy(d,j)  = std(vals, 'omitnan') / sqrt(sum(~isnan(vals)));
        % Save the per-rat means into cell arrays
        if strcmp(uniqueAges{j}, 'young')
            dataEntropyYoung{d} = vals;
        elseif strcmp(uniqueAges{j}, 'old')
            dataEntropyOld{d} = vals;
        end
    end
end


%---------------- PLOT PER RAT PER DAY MEAN ENTROPY ---------------------------%
fe2=figure;
hold on;
plot_bar_sem_WaterMaze(uniqueDays, meanEntropy(:,1), semEntropy(:,1), ...
    meanEntropy(:,2), semEntropy(:,2), dataEntropyYoung, dataEntropyOld, postHocResults);
title('Entropy (Rats Day Means)  - Grouped By Strategy Type');
xlabel('Day');
ylabel('Entropy');
pubify_figure_axis_robust(14,14);
xticks(uniqueDays);
hold off;
saveas(fe2, fullfile(fig_dir,'Entropy', 'GroupStrat_Entropy_MeanRat'), 'png');
close;

% Plot all rats as line plot
fe3 = figure;
hold on;
nRats = size(meanDayEntropy,1);  % Number of rats
for r = 1:nRats
    if strcmp(AgeCell{r}, 'young')
        plot(uniqueDays, meanDayEntropy(r,:), '-o', 'Color', youngColor, 'LineWidth', 1.5);
    else
        plot(uniqueDays, meanDayEntropy(r,:), '-o', 'Color', oldColor, 'LineWidth', 1.5);
    end
end
xlabel('Day');
ylabel('Entropy');
title('Entropy Per Rat');
pubify_figure_axis_robust(14,14);
xticks(uniqueDays);
hold off;
saveas(fe3, fullfile(fig_dir, 'Entropy', 'GroupStrat_Entropy_RatChanges'), 'png');
close;
%--------------- PLOT PER RAT PER DAY MEAN ENTROPY ----------------------%
fe4=figure;
hold on;
% Calculate dataEntropy for all trials
dataEntropyYoung = cell(numel(uniqueDays),1);
dataEntropyOld   = cell(numel(uniqueDays),1);
for d = 1:numel(uniqueDays)
    dataEntropyYoung{d}=entropy_vals(uniqueDays(d) & strcmp(data1.Age,'young'));
    dataEntropyOld{d}=entropy_vals(data1.Day == uniqueDays(d) & strcmp(data1.Age,'old'));
end
plot_bar_sem_WaterMaze(uniqueDays, meanEntropy(:,1), semEntropy(:,1), ...
    meanEntropy(:,2), semEntropy(:,2), dataEntropyYoung, dataEntropyOld, postHocResults);
title(' Entropy (All Trials) - Grouped By Strategy Type');
xlabel('Day');
ylabel('Entropy');
pubify_figure_axis_robust(14,14);
xticks(uniqueDays)
hold off;
saveas(fe4, fullfile(fig_dir, 'Entropy', 'GroupStrat_Entropy_Mean_AllTrials'), 'png');
close;



%% Change platform groups -test
% 3 platform groups
data1.pd(data1.pd==4)=3; % Add platform 5 to the same group as 4&6
% 2 Platform groups
data1.pd(data1.pd==2)=1; % Group 2,8, 4,7 in one group

%% 7) Individual Strategy x Platform Group

%Unique variablers
uniqueDays=unique(data1.Day);
uniqueRats=unique(data1.x_TargetID);
uniquePlatforms=unique(data1.pd);
uniquePlatforms = uniquePlatforms(~isnan(uniquePlatforms)); % Remove NaNs to avoid duplicate names

% Colors for plotting
% Base colors
youngColor = [0.2196, 0.5569, 0.2353]; % green
oldColor   = [0.4157, 0.1059, 0.6039]; % purple

% Generate 4 shades from light to dark by blending with white and darkening
blendLevels = linspace(0.6, 0, 4)';  % from light to base color

% Blend with white
youngShades = blendLevels .* 1 + (1 - blendLevels) .* youngColor;
oldShades   = blendLevels .* 1 + (1 - blendLevels) .* oldColor;
% Base colors
clrYoungSig = [0.2196, 0.5569, 0.2353]; % green
clrOldSig   = [0.4157, 0.1059, 0.6039]; % purple
clrBothSig  = [0.3 0.3 0.3];% gray

%Loop over strategies
for s = 1:numel(strategyNames)
    stratName = strategyNames{s};
    stratTitle = strategy_titles{s};

    % Determine number of platforms dynamically
    nPlat = numel(uniquePlatforms);
    % Build measure variable names as cell array, e.g. {'P1','P2','P3','P4'}.
    measureVars = arrayfun(@(x) sprintf('P%d', x), uniquePlatforms, 'UniformOutput', false);

    % Initialize master cell array: each row will be {RatID, Age, Day, P1, P2, ... PnPlat}
    master_data = {};
    % Loop over each day
    for d = 1:numel(uniqueDays)
        dVal = uniqueDays(d);
        dayData = data1(data1.Day == dVal, :);

        % Loop over each rat (subject) present in dayData
        for r = 1:numel(uniqueRats)
            ratID = uniqueRats{r};
            isCurrentRat = strcmp(dayData.x_TargetID, ratID);
            if ~any(isCurrentRat)
                continue; % rat not present on this day
            end
            % Get the rat's Age (assumed constant within rat)
            ratAge = unique(dayData.Age(isCurrentRat));
            if iscell(ratAge)
                ratAge = ratAge{1};
            end

            % For each platform, compute the mean strategy usage for this rat on this day.
            meanUse = nan(1, nPlat);
            for p = 1:nPlat
                pVal = uniquePlatforms(p);
                idx = isCurrentRat & (dayData.pd == pVal);
                if any(idx)
                    meanUse(p) = mean(dayData.(stratName)(idx));
                end
            end

            % Create a row cell with 3 fixed columns then nPlat columns.
            row = {ratID, ratAge, dVal};
            for p = 1:nPlat
                row = [row, {meanUse(p)}];
            end
            master_data(end+1,:) = row;
        end
    end

    % Build variable names: fixed names and then the dynamic measureVars
    measureVars = arrayfun(@(x) sprintf('P%d', x), uniquePlatforms, 'UniformOutput', false);
    measureVars = measureVars(:)';  % force row vector
    varNames = [{'RatID','Age','Day'}, measureVars];
    master_table = cell2table(master_data, 'VariableNames', varNames);

    %-------- Run ANOVA & Tukey Post Hoc for each day and combine results
    all_anova = table();
    all_tukey = table();
    for d = 1:numel(uniqueDays)
        dVal = uniqueDays(d);
        day_table = master_table(master_table.Day == dVal, :);
        % Create a reduced table with only Age and the repeated measures variables.
        anova_tbl = day_table(:, [{'Age'}, measureVars]);

        % Run repeated measures ANOVA and post hoc tests on this reduced table.
        anovaResults = runMixedANOVA(anova_tbl, measureVars,'Platform');
        postHocResults = runTukeyPostHocMixed(anova_tbl, measureVars,'Platform');


        % Add a column to indicate Day
        anovaResults.Day = repmat(dVal, height(anovaResults), 1);
        postHocResults.Day = repmat(dVal, height(postHocResults), 1);
        % Remove row names for concatenation
        anovaResults.Properties.RowNames = {};
        postHocResults.Properties.RowNames = {};

        all_anova = [all_anova; anovaResults];
        all_tukey = [all_tukey; postHocResults];
    end

    % Save combined ANOVA & Tukey results for this strategy.
    anovaFile = fullfile(processed_dir, sprintf('3%s_PlatformGroups_ANOVA.csv', stratTitle));
    tukeyFile = fullfile(processed_dir, sprintf('3%s_PlatformGroups_Tukey.csv', stratTitle));
    writetable(all_anova, anovaFile);
    writetable(all_tukey, tukeyFile);

    %-----------------Compute summary statistics for plotting
    nDays = numel(uniqueDays);
    % Preallocate matrices for young and old (nDays x nPlat)
    meanYoung = nan(nDays, nPlat);
    semYoung  = nan(nDays, nPlat);
    meanOld   = nan(nDays, nPlat);
    semOld    = nan(nDays, nPlat);
    dataYoung = cell(nDays, nPlat);
    dataOld   = cell(nDays, nPlat);

    for d = 1:nDays
        dVal = uniqueDays(d);
        day_tbl = master_table(master_table.Day == dVal, :);
        for p = 1:nPlat
            colName = measureVars{p};  % e.g., 'P1', 'P2', etc.
            % For young rats:
            yVals = day_tbl{strcmp(day_tbl.Age, 'young'), colName};
            meanYoung(d, p) = mean(yVals, 'omitnan');
            semYoung(d, p) = std(yVals, 'omitnan') / sqrt(sum(~isnan(yVals)));
            dataYoung{d, p} = yVals;
            % For old rats:
            oVals = day_tbl{strcmp(day_tbl.Age, 'old'), colName};
            meanOld(d, p) = mean(oVals, 'omitnan');
            semOld(d, p) = std(oVals, 'omitnan') / sqrt(sum(~isnan(oVals)));
            dataOld{d, p} = oVals;
        end
    end

    % Build grouped bar matrices: for each day, create 2*nPlat bars (odd: young, even: old)
    meanMatrix = nan(nDays, 2*nPlat);
    semMatrix = nan(nDays, 2*nPlat);
    for p = 1:nPlat
        meanMatrix(:, 2*p-1) = meanYoung(:, p);
        meanMatrix(:, 2*p)   = meanOld(:, p);
        semMatrix(:, 2*p-1) = semYoung(:, p);
        semMatrix(:, 2*p)   = semOld(:, p);
    end

    %-------------Plotting the Grouped Bar Chart with Sigstar
    fH = figure('Name', stratTitle, 'Color', [1 1 1]);
    sgtitle(stratTitle, 'FontSize', 14, 'FontWeight', 'bold');
    hold on;
    hBar = bar(meanMatrix, 'grouped');
    set(gca, 'XTick', 1:nDays, 'XTickLabel', arrayfun(@num2str, uniqueDays, 'UniformOutput', false));
    xlabel('Day'); ylabel('Mean Strategy Usage');

    % Assign colors: for each platform, young bars get youngShades; old bars get oldShades.
    for p = 1:nPlat
        colYoung = 2*p-1;
        colOld = 2*p;
        idx = min(p, size(youngShades,1));  % if nPlat > number of defined colors, use the last available color
        hBar(colYoung).FaceColor = youngShades(idx, :);
        hBar(colYoung).FaceAlpha = 0.75;
        hBar(colOld).FaceColor = oldShades(idx, :);
        hBar(colOld).FaceAlpha = 0.75;
    end


    % Add error bars.
    drawnow; % ensure bar positions are updated
    for c = 1:size(meanMatrix,2)
        xVals = hBar(c).XEndPoints;
        yVals = meanMatrix(:, c);
        eVals = semMatrix(:, c);
        errorbar(xVals, yVals, eVals, 'k', 'LineStyle', 'none', 'LineWidth', 1.2);
    end

    % Add scatter points (jittered) for individual rat values.
    jitterAmount = 0.05;
    [nDays, nCols] = size(meanMatrix); % nCols = 2*nPlat
    for d = 1:nDays
        for c = 1:nCols
            if mod(c,2)==1
                % Odd col => Young
                pIdx = (c+1)/2; % e.g. c=1 => pIdx=1, c=3 => pIdx=2
                pts  = dataYoung{d, pIdx};
            else
                % Even col => Old
                pIdx = c/2;
                pts  = dataOld{d, pIdx};
            end
            if isempty(pts), continue; end

            xCenter = hBar(c).XEndPoints(d);
            xJitter = xCenter + (rand(size(pts))-0.5)*jitterAmount;

            scatter(xJitter, pts, 12, ...
                'MarkerFaceColor', hBar(c).FaceColor, ...
                'MarkerEdgeColor','none', 'MarkerFaceAlpha',1);
        end
    end

    % Add significance markers using sigstar.
    % Here we assume that the post hoc results (all_tukey) contain fields that indicate which platforms are compared.
    % For example, assume they have fields 'Platform1' and 'Platform2' (as indices 1:nPlat),
    % a field 'pValue', and a field 'Age' indicating 'young', 'old', or 'both'.
    sigPairs_all = {};
    sigPvals_all = [];
    sigColors_all = {};
    for iRow = 1:height(all_tukey)
        day_val = all_tukey.Day(iRow);
        % Find the day index:
        dayIdx = find(uniqueDays == day_val);
        if isempty(dayIdx)
            continue;
        end
        % Use the fields 'Platform1' and 'Platform2'
        p1 = all_tukey.Platform_1(iRow);
        p2 = all_tukey.Platform_2(iRow);
        pVal = all_tukey.pValue(iRow);
        if pVal >= 0.05
            continue;
        end
        compAge = all_tukey.Age{iRow};  % expected 'young', 'old', or 'both'
        if strcmpi(compAge, 'young')
            x1 = hBar(2*p1-1).XEndPoints(dayIdx);
            x2 = hBar(2*p2-1).XEndPoints(dayIdx);
            sigColor = clrYoungSig;
        elseif strcmpi(compAge, 'old')
            x1 = hBar(2*p1).XEndPoints(dayIdx);
            x2 = hBar(2*p2).XEndPoints(dayIdx);
            sigColor = clrOldSig;
        else
            x_candidates = [hBar(2*p1-1).XEndPoints(dayIdx), hBar(2*p1).XEndPoints(dayIdx), ...
                hBar(2*p2-1).XEndPoints(dayIdx), hBar(2*p2).XEndPoints(dayIdx)];
            x1 = min(x_candidates);
            x2 = max(x_candidates);
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
    %Legend
    % legend_labels = cell(1, 2*nPlat);
    % platform_names={'P1','P2','P3','P4'};
    % for p = 1:nPlat
    %     legend_labels{2*p-1} = sprintf('%s: Young', platform_names{p});
    %     legend_labels{2*p}   = sprintf('%s: Old', platform_names{p});
    % end
    %legend(hBar, legend_labels, 'Location', 'bestoutside');
    pubify_figure_axis_robust(14,14);
    % Save the figure
    saveas(fH, fullfile(fig_dir,'Platform', sprintf('3Platform_%s.png', stratTitle)));
    close(fH);
end



%% 8) Group Strategy x Platform Group
strategyGroups = {
    {'thigmotaxis', 'circling', 'random_path'}, 'Non-Goal Oriented';
    {'scanning', 'chaining'}, 'Procedural';
    {'directed_search', 'corrected_search', 'direct_path','perseverance'}, 'Allocentric'
    };
groupNames = strategyGroups(:, 2);  % Extract group labels
nGroups = numel(groupNames);
% Numel & unique variables
groupStrat = cell(1, nGroups);
uniqueDays=unique(data1.Day);
uniqueRats=unique(data1.x_TargetID);
uniquePlatforms=unique(data1.pd);
uniquePlatforms = uniquePlatforms(~isnan(uniquePlatforms)); % Remove NaNs
% Colors for plotting
% Base colors
clrYoungSig = [0.2196, 0.5569, 0.2353]; % green
clrOldSig   = [0.4157, 0.1059, 0.6039]; % purple

% Generate 4 shades from light to dark by blending with white and darkening
blendLevels = linspace(0.6, 0, 4)';  % from light to base color

% Blend with white
youngShades = blendLevels .* 1 + (1 - blendLevels) .* clrYoungSig;
oldShades   = blendLevels .* 1 + (1 - blendLevels) .* clrOldSig;
clrBothSig  = [0.3 0.3 0.3];% gray
% youngShades = [0.80 1.00 0.80; 0.60 0.85 0.60; 0.30 0.70 0.30; 0.00 0.60 0.00];
% oldShades   = [0.85 0.65 0.85; 0.75 0.40 0.75; 0.60 0.25 0.60; 0.50 0.00 0.50];
% clrYoungSig = [0 0.6 0];    % green
% clrOldSig   = [0.5 0 0.5];  % purple
% clrBothSig  = [0.3 0.3 0.3];% gray
% Loop over each strategy group to compute the aggregated group probability
for g = 1:nGroups
    currentStrategies = strategyGroups{g,1};
    % Start with the first strategy's probability
    groupProb = data1.(currentStrategies{1});
    % Sum probabilities across all strategies in the group
    for ii = 2:numel(currentStrategies)
        groupProb = groupProb + data1.(currentStrategies{ii});
    end
    groupStrat{g} = groupProb;
    groupTitle=groupNames{g};
    %---------------------- Build Master Table per Group ----------------------%
    % For each group, create a master table where each row = {RatID, Age, Day, P1, P2, ..., Pn}
    nPlat = numel(uniquePlatforms);
    master_data = {};
    for d = 1:numel(uniqueDays)
        dVal = uniqueDays(d);
        % Subset data1 for the current day
        dayData = data1(data1.Day == dVal, :);
        % Get the corresponding aggregated group probabilities for the current day.
        % (Assumes data1 rows remain in the same order.)
        groupProbDay = groupStrat{g}(data1.Day == dVal);
        for r = 1:numel(uniqueRats)
            ratID = uniqueRats{r};
            isCurrentRat = strcmp(dayData.x_TargetID, ratID);
            if ~any(isCurrentRat)
                continue; % This rat is not present on the current day
            end
            % Get the rat's Age (assumed constant)
            ratAge = unique(dayData.Age(isCurrentRat));
            if iscell(ratAge)
                ratAge = ratAge{1};
            end
            % For each platform, compute the mean aggregated probability for that rat on the day.
            meanUse = nan(1, nPlat);
            for p = 1:nPlat
                pVal = uniquePlatforms(p);
                idx = isCurrentRat & (dayData.pd == pVal);
                if any(idx)
                    meanUse(p) = mean(groupProbDay(idx));
                end
            end
            % Create a row: {RatID, Age, Day, P1, P2, ..., PnPlat}
            row = {ratID, ratAge, dVal};
            for p = 1:nPlat
                row = [row, {meanUse(p)}];
            end
            master_data(end+1,:) = row;
        end
    end

    % Build variable names: fixed columns and then one per platform
    measureVars = arrayfun(@(x) sprintf('P%d', x), uniquePlatforms, 'UniformOutput', false);
    measureVars = measureVars(:)';  % force row vector
    varNames = [{'RatID','Age','Day'}, measureVars];
    master_table = cell2table(master_data, 'VariableNames', varNames);


    %------------ Run ANOVA & Tukey Post Hoc for each day and combine results
    all_anova = table();
    all_tukey = table();
    for d = 1:numel(uniqueDays)
        dVal = uniqueDays(d);
        day_table = master_table(master_table.Day == dVal, :);
        % Create a reduced table with only the between-subject variable and repeated measures
        anova_tbl = day_table(:, [{'Age'}, measureVars]);

        % Run repeated measures ANOVA (within factor = Platform, between factor = Age)
        anovaResults = runMixedANOVA(anova_tbl, measureVars, 'Platform');
        postHocResults = runTukeyPostHocMixed(anova_tbl, measureVars, 'Platform');

        % Tag results with the current Day
        anovaResults.Day = repmat(dVal, height(anovaResults), 1);
         disp(anovaResults)
        postHocResults.Day = repmat(dVal, height(postHocResults), 1);
        % Remove row names to allow concatenation
        anovaResults.Properties.RowNames = {};
       
        postHocResults.Properties.RowNames = {};

        all_anova = [all_anova; anovaResults];
        all_tukey = [all_tukey; postHocResults];
    end

    % Save the combined ANOVA & Tukey results for this group.
    anovaFile = fullfile(processed_dir, sprintf('3%s_PlatformGroups_ANOVA.csv', groupTitle));
    tukeyFile = fullfile(processed_dir, sprintf('3%s_PlatformGroups_Tukey.csv', groupTitle));
    writetable(all_anova, anovaFile);
    writetable(all_tukey, tukeyFile);

    %Compute summary statistics for plotting
    % Compute means and SEMs for each Day, separately for young and old, for each platform.
    nDays = numel(uniqueDays);
    meanYoung = nan(nDays, nPlat);
    semYoung = nan(nDays, nPlat);
    meanOld = nan(nDays, nPlat);
    semOld = nan(nDays, nPlat);
    dataYoung = cell(nDays, nPlat);
    dataOld   = cell(nDays, nPlat);

    for d = 1:nDays
        dVal = uniqueDays(d);
        day_tbl = master_table(master_table.Day == dVal, :);
        for p = 1:nPlat
            colName = measureVars{p};  % e.g., 'P1', 'P2', ...
            % For young rats:
            yVals = day_tbl{strcmp(day_tbl.Age, 'young'), colName};
            meanYoung(d, p) = mean(yVals, 'omitnan');
            semYoung(d, p) = std(yVals, 'omitnan')/sqrt(sum(~isnan(yVals)));
            dataYoung{d, p} = yVals;
            % For old rats:
            oVals = day_tbl{strcmp(day_tbl.Age, 'old'), colName};
            meanOld(d, p) = mean(oVals, 'omitnan');
            semOld(d, p) = std(oVals, 'omitnan')/sqrt(sum(~isnan(oVals)));
            dataOld{d, p} = oVals;
        end
    end

    % Create matrices for grouped bars.
    % For each day, we create 2*nPlat bars (odd indices for young, even for old).
    meanMatrix = nan(nDays, 2*nPlat);
    semMatrix = nan(nDays, 2*nPlat);
    for p = 1:nPlat
        meanMatrix(:, 2*p-1) = meanYoung(:, p);
        meanMatrix(:, 2*p)   = meanOld(:, p);
        semMatrix(:, 2*p-1) = semYoung(:, p);
        semMatrix(:, 2*p)   = semOld(:, p);
    end

    %Plotting the Grouped Bar Chart with Sigstar
    fH = figure('Name', groupTitle, 'Color', [1 1 1]);
    sgtitle(groupTitle, 'FontSize', 14, 'FontWeight', 'bold');
    ylim([0 1.5])
    hold on;
    hBar = bar(meanMatrix, 'grouped');
    set(gca, 'XTick', 1:nDays, 'XTickLabel', arrayfun(@num2str, uniqueDays, 'UniformOutput', false));
    xlabel('Day'); ylabel('Mean Aggregated Strategy Usage');

    % Assign colors for each platform: young bars get youngShades; old bars get oldShades.
    for p = 1:nPlat
        colYoung = 2*p-1;
        colOld = 2*p;
        idx = min(p, size(youngShades,1));
        hBar(colYoung).FaceColor = youngShades(idx, :);
        hBar(colYoung).FaceAlpha = 0.75;
        hBar(colOld).FaceColor = oldShades(idx, :);
        hBar(colOld).FaceAlpha = 0.75;
    end

    % Add error bars.
    drawnow;
    for c = 1:size(meanMatrix,2)
        xVals = hBar(c).XEndPoints;
        yVals = meanMatrix(:, c);
        eVals = semMatrix(:, c);
        errorbar(xVals, yVals, eVals, 'k', 'LineStyle', 'none', 'LineWidth', 1.2);
    end

    % Add scatter points (jittered) for individual rat values.
    jitterAmount = 0.05;
    for d = 1:nDays
        for c = 1:2*nPlat
            if mod(c,2)==1
                pIdx = ceil(c/2);
                pts = dataYoung{d, pIdx};
            else
                pIdx = ceil(c/2);
                pts = dataOld{d, pIdx};
            end
            if isempty(pts), continue; end
            xCenter = hBar(c).XEndPoints(d);
            xJitter = xCenter + (rand(size(pts))-0.5)*jitterAmount;
            scatter(xJitter, pts, 12, 'MarkerFaceColor', hBar(c).FaceColor, ...
                'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 1);
        end
    end

    % Add significance markers using sigstar.
    sigPairs_all = {};
    sigPvals_all = [];
    sigColors_all = {};
    for iRow = 1:height(all_tukey)
        day_val = all_tukey.Day(iRow);
        % Find the day index:
        dayIdx = find(uniqueDays == day_val);
        if isempty(dayIdx)
            continue;
        end
        % Here, we assume the fields 'Day_1' and 'Day_2' actually indicate Platform numbers.
        p1 = all_tukey.Platform_1(iRow);
        p2 = all_tukey.Platform_2(iRow);
        pVal = all_tukey.pValue(iRow);
        if pVal >= 0.05
            continue;
        end
        % Determine which age group the comparison applies to.
        compAge = all_tukey.Age{iRow};  % expected 'young', 'old', or 'both'
        if strcmpi(compAge, 'young')
            x1 = hBar(2*p1-1).XEndPoints(dayIdx);
            x2 = hBar(2*p2-1).XEndPoints(dayIdx);
            sigColor = clrYoungSig;
        elseif strcmpi(compAge, 'old')
            x1 = hBar(2*p1).XEndPoints(dayIdx);
            x2 = hBar(2*p2).XEndPoints(dayIdx);
            sigColor = clrOldSig;
        else
            x_candidates = [hBar(2*p1-1).XEndPoints(dayIdx), hBar(2*p1).XEndPoints(dayIdx), ...
                hBar(2*p2-1).XEndPoints(dayIdx), hBar(2*p2).XEndPoints(dayIdx)];
            x1 = min(x_candidates);
            x2 = max(x_candidates);
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
    % %Legend
    % legend_labels = cell(1, 2*nPlat);
    % platform_names={'P1','P2','P3','P4'};
    % for p = 1:nPlat
    %     legend_labels{2*p-1} = sprintf('%s: Young', platform_names{p});
    %     legend_labels{2*p}   = sprintf('%s: Old', platform_names{p});
    % end
    % legend(hBar, legend_labels, 'Location', 'bestoutside');
    %
    pubify_figure_axis_robust(14,14);
    % Save the figure
    saveas(fH, fullfile(fig_dir, 'Platform',sprintf('3Platform_%s.png', groupTitle)));
    close(fH);
end


%% 9) Strategy Group Diff for young and old separately

% Strategy groups
strategyGroups = {
    {'thigmotaxis','circling','random_path'}, 'NonGoal';
    {'scanning','chaining'},                 'Procedural';
    {'directed_search','corrected_search','direct_path','perseverance'}, 'Allocentric'
    };
groupLabels = strategyGroups(:,2);       % {'NonGoal','Procedural','Allocentric'}
groupNames  = {'Non‑Goal Oriented','Procedural','Allocentric'}; % for display
nGroups     = numel(groupLabels);
uniqueDays  = unique(data1.Day);

% Colors for the three groups
stratColors = [

0.3961    0.2627    0.1294;   % Non‑Goal  brown
1.0000    0.7020    0.4000;   % Procedural yellow ochre
0    0.5020         0    % Allocentric green
];

rats = unique(data1.x_TargetID);
nRats = numel(rats);
nDays = numel(uniqueDays);

% Preallocate per‑rat matrix [nRats × (nGroups×nDays)]
ratMeans = nan(nRats, nGroups * nDays);
AgeCell  = cell(nRats,1);

for i = 1:nRats
    ridx = strcmp(data1.x_TargetID, rats{i});
    AgeCell{i} = data1.Age{find(ridx,1)};
    col = 1;
    for g = 1:nGroups
        % Sum the raw probabilities of the strategies in this group
        gp = zeros(sum(ridx),1);
        for s = strategyGroups{g,1}
            gp = gp + data1.(s{1})(ridx);
        end
        % Mean per day
        for d = 1:nDays
            pdx = ridx & (data1.Day==uniqueDays(d));
            ratMeans(i, col) = mean(gp(data1.Day(ridx)==uniqueDays(d))); %#ok<*FNDSB>
            col = col + 1;
        end
    end
end

% Build the table: columns = RatID | Age | <Group>_D1 ... <Group>_D4 for each group
varNames = {'RatID','Age'};
for g = 1:nGroups
    for d = 1:nDays
        varNames{end+1} = sprintf('%s_Day%d', groupLabels{g}, uniqueDays(d)); %#ok<*AGROW>
    end
end
ratTable = table(rats, AgeCell, ratMeans(:,1), ratMeans(:,2), ratMeans(:,3), ratMeans(:,4), ...
    ratMeans(:,5), ratMeans(:,6), ratMeans(:,7), ratMeans(:,8), ...
    ratMeans(:,9), ratMeans(:,10), ratMeans(:,11), ratMeans(:,12), ...
    'VariableNames', varNames);

%Repeated‑Measures ANOVA & Tukey Post Hoc, per Age Group

% Prepare WithinDesign: 12 rows = 3 groups × 4 days
% assume uniqueDays is [1;2;3;4] and groupLabels = {'NonGoal','Procedural','Allocentric'}
nDays    = numel(uniqueDays);
nGroups  = numel(groupLabels);

DaysVec  = repmat(uniqueDays, nGroups, 1);        % [1;2;3;4;1;2;3;4;1;2;3;4]
StrVec   = repelem(groupLabels(:), nDays, 1);    % {'NonGoal';'NonGoal';…4×; 'Procedural';…4×; 'Allocentric';…4×}

withinDesign = table(DaysVec, categorical(StrVec), ...
    'VariableNames', {'Day','StrategyGroup'});

measureVars = ratTable.Properties.VariableNames(3:end); % the 12 columns

for age = {'young','old'}
    ag = age{1};
    idxAge = strcmp(ratTable.Age, ag);
    subTbl  = ratTable(idxAge, :);

    % Fit RM model: no between factor, two within factors Day & StrategyGroup
    formula = sprintf('%s-%s ~ 1', measureVars{1}, measureVars{end});
    rm = fitrm(subTbl, formula, 'WithinDesign', withinDesign);

    % ANOVA: test main effects & interaction
    anovaResults = ranova(rm, 'WithinModel', 'Day*StrategyGroup');
    writetable(anovaResults, fullfile(processed_dir, sprintf('ANOVA_%s.csv', ag)));

    % Tukey post hoc for between‑strategy (by Day)
    postHocStrat = multcompare(rm, 'StrategyGroup', 'By', 'Day', 'ComparisonType', 'tukey-kramer');
    writetable(postHocStrat, fullfile(processed_dir, sprintf('PostHoc_Strat_%s.csv', ag)));

    % Tukey post hoc for within‑strategy (by StrategyGroup)
    postHocDay   = multcompare(rm, 'Day', 'By', 'StrategyGroup', 'ComparisonType', 'tukey-kramer');
    writetable(postHocDay,   fullfile(processed_dir, sprintf('PostHoc_Day_%s.csv',   ag)));

    %Plot: Grouped Bar Chart with SEM, Jittered Per‑Rat Means, sigstar
    % Compute group means & SEM across rats
    meanData = zeros(nDays, nGroups);
    semData  = zeros(nDays, nGroups);
    dataMeansPerGroup = cell(nGroups, nDays);
    for g = 1:nGroups
        for d = 1:nDays
            colName = sprintf('%s_Day%d', groupLabels{g}, uniqueDays(d));
            vals = subTbl.(colName);
            meanData(d,g) = mean(vals);
            semData(d,g)  = std(vals)/sqrt(numel(vals));
            dataMeansPerGroup{g,d} = vals;
        end
    end

    % Plot
    f = figure('Position',[95,100,1100,630]);
    hold on;
    hBar = bar(uniqueDays, meanData, 'grouped', 'BarWidth', 0.8);
    for g = 1:nGroups
        hBar(g).FaceColor = stratColors(g,:);
        hBar(g).FaceAlpha = 0.35;
    end
    drawnow;
    % Get bar centers
    barCenters = reshape([hBar.XEndPoints], nGroups, nDays)' ;
    %Re‐compute bar centers so they’re aligned
    barCenters = nan(nDays, nGroups);
    for g = 1:nGroups
        barCenters(:,g) = hBar(g).XEndPoints';
    end
    % Error bars
    for g = 1:nGroups
        errorbar(barCenters(:,g), meanData(:,g), semData(:,g), 'k', 'LineStyle','none','LineWidth',1.5);
    end

    % Jittered per-rat means
    jit = 0.08;
    for g = 1:nGroups
        for d = 1:nDays
            pts = dataMeansPerGroup{g,d};
            xj  = barCenters(d,g) + (rand(size(pts))-0.5)*jit;
            scatter(xj, pts, 20, 'MarkerFaceColor', stratColors(g,:), ...
                'MarkerEdgeColor','none','MarkerFaceAlpha',0.8);
        end
    end

    %Prepare sigstar inputs
    sigPairs = {};
    sigP     = [];
    sigCols  = {};

    % 1. Between‐strategy comparisons on each day (gray)
    cs = multcompare(rm, 'StrategyGroup', 'By', 'Day', 'ComparisonType', 'tukey-kramer');
    % cs is a table: cs.Day, cs.StrategyGroup_1, cs.StrategyGroup_2, cs.pValue, etc.
    for rr = 1:height(cs)
        pval = cs.pValue(rr);
        if pval < 0.05
            % get day index
            dVal = cs.Day(rr);
            dIdx = find(uniqueDays == dVal, 1);
            % find group indices
            lvl1 = char(cs.StrategyGroup_1(rr));
            lvl2 = char(cs.StrategyGroup_2(rr));
            g1 = find(strcmp(groupLabels, lvl1), 1);
            g2 = find(strcmp(groupLabels, lvl2), 1);
            % only plot once per unordered pair
            if ~isempty(dIdx) && ~isempty(g1) && ~isempty(g2) && (g1 < g2)
                sigPairs{end+1} = [barCenters(dIdx, g1), barCenters(dIdx, g2)];
                sigP(end+1)     = pval;
                sigCols{end+1}  = [0.5,0.5,0.5];  % gray
            end
        end
    end
    if ~isempty(sigPairs)
        hSig = sigstar(sigPairs, sigP);
        for k = 1:size(hSig,1)
            set(hSig(k,1), 'Color', sigCols{k});   % line
            set(hSig(k,2), 'Color', sigCols{k});   % stars
        end
    end


    % % 2. Within‐strategy consecutive‐day comparisons (black)
    % Recompute bar centers
    barCenters = nan(nDays, nGroups);
    for g = 1:nGroups
        barCenters(:,g) = hBar(g).XEndPoints';
    end

    % Preallocate
    sigPairs = {};
    sigP     = [];
    sigCols  = {};
    nsAnn    = struct('x',{},'y',{},'text',{},'color',{});

    % (Re)compute bar centers
    barCenters = nan(nDays, nGroups);
    for g = 1:nGroups
        barCenters(:,g) = hBar(g).XEndPoints';
    end

    % Prepare containers
    sigPairs = {};
    sigP     = [];
    sigCols  = {};

    % Gray color for 2‑out‑of‑3 lines
    darkGray = [0.2 0.2 0.2];

    % Run Tukey for Day within StrategyGroup
    cd = multcompare(rm, 'Day','By','StrategyGroup','ComparisonType','tukey-kramer');

    for d = 1:nDays-1
        % select only consecutive‐day rows
        sel = cd.Day_1==uniqueDays(d) & cd.Day_2==uniqueDays(d+1);
        sub = cd(sel,:);
        % build p‑vector
        pv = nan(1,nGroups);
        for rr = 1:height(sub)
            gidx = find(strcmp(groupLabels, char(sub.StrategyGroup(rr))),1);
            pv(gidx) = sub.pValue(rr);
        end
        sigIdx = find(pv<0.05);
        switch numel(sigIdx)
            case 0
                % none significant → do nothing
            case 3
                % all three → black line between avg centers
                x1 = mean(barCenters(d, :));
                x2 = mean(barCenters(d+1, :));
                sigPairs{end+1} = [x1, x2];
                sigP(end+1)     = min(pv(sigIdx));
                sigCols{end+1}  = [0 0 0];
            case 2
                % two significant → dark gray line
                x1 = mean(barCenters(d, sigIdx));
                x2 = mean(barCenters(d+1, sigIdx));
                sigPairs{end+1} = [x1, x2];
                sigP(end+1)     = min(pv(sigIdx));
                sigCols{end+1}  = darkGray;
                % annotate the non‑significant with 'ns' in its own color
                ns = setdiff(1:nGroups, sigIdx);
                yLine = max(max(meanData(d, sigIdx)), max(meanData(d+1, sigIdx))) + 0.02;
                x_ns  = barCenters(d, ns);
                % text(x_ns+0.5, repmat(yLine+0.55, size(x_ns)), 'n.s.', ...
                %     'Color', stratColors(ns,:), 'HorizontalAlignment','Center');
            case 1
                % single significant → colored line
                g = sigIdx;
                x1 = barCenters(d, g);
                x2 = barCenters(d+1, g);
                sigPairs{end+1} = [x1, x2];
                sigP(end+1)     = pv(g);
                sigCols{end+1}  = stratColors(g,:);
        end
    end

    % filter out invalid entries
    valid = cellfun(@(c)isnumeric(c)&&numel(c)==2, sigPairs);
    sigPairs = sigPairs(valid);
    sigP     = sigP(valid);
    sigCols  = sigCols(valid);

    % draw significance bars & stars
    if ~isempty(sigPairs)
        hSig = sigstar(sigPairs, sigP);
        for k = 1:size(hSig,1)
            set(hSig(k,1), 'Color', sigCols{k});   % line
            set(hSig(k,2), 'Color', sigCols{k});   % stars
        end
    end


    % Final formatting
    title(sprintf('%s', upper(ag)));
    xlabel('Day');
    ylabel('Probability for Strategy Group');
    xticks(uniqueDays);
    legend(groupNames, 'Location','northeastoutside');
    ylim([0 1.2]);
    pubify_figure_axis_robust(16,16);
    hold off;

    % Save figure
    saveas(f, fullfile(fig_dir, sprintf('DayStrategyUse_%s.png', ag)));
    close(f);
end

%% 10) Individual strategy differences per dya per age group
% List the eight (or however many) raw strategy probability columns
strategyList = { ...
    'thigmotaxis', 'circling',      'random_path', ...
    'scanning',    'chaining', ...
    'directed_search','corrected_search','direct_path' };

% 8 distinct strategy colours  (rows = strategies in your strategyList order)
stratColors = [
    0.850, 0.325, 0.098;   % thigmotaxis      – vermilion
    0.094, 0.509, 0.800;   % circling         – blue
    0.467, 0.675, 0.188;   % random_path      – green
    0.800, 0.475, 0.655;   % scanning         – purple-pink
    0.894, 0.766, 0.039;   % chaining         – yellow-gold
    0.436, 0.600, 0.800;   % directed_search  – sky blue
    0.905, 0.161, 0.541;   % corrected_search – magenta
    0.400, 0.400, 0.400    % direct_path      – gray
    ];


nStrats     = numel(strategyList);      % e.g. 9 strategies
uniqueDays  = unique(data1.Day);
nDays       = numel(uniqueDays);

%  display names & colours (edit as you like)

rats   = unique(data1.x_TargetID);
nRats  = numel(rats);

% ---------- Build a per-rat matrix: rows = rats, columns = (strategy × day)
ratMeans = nan(nRats, nStrats * nDays);
AgeCell  = cell(nRats,1);

for i = 1:nRats
    ridx        = strcmp(data1.x_TargetID, rats{i});
    AgeCell{i}  = data1.Age{find(ridx,1)};

    col = 1;                       % running column index
    for sIdx = 1:nStrats
        stratName = strategyList{sIdx};
        stratVec  = data1.(stratName)(ridx);    % all trials for that rat

        for d = 1:nDays
            dayMask = (data1.Day(ridx)==uniqueDays(d));
            ratMeans(i,col) = mean(stratVec(dayMask));
            col = col + 1;
        end
    end
end

% ---------- Build variable names like  <Strategy>_Day1 … _Day4
varNames = {'RatID','Age'};
for sIdx = 1:nStrats
    for d = 1:nDays
        varNames{end+1} = sprintf('%s_Day%d', strategyList{sIdx}, uniqueDays(d));
    end
end

ratTable = cell2table([rats, AgeCell, num2cell(ratMeans)], ...
    'VariableNames', varNames);

% ----------  Within-Design table for fitrm  (nStrats*nDays rows)
DaysVec  = repmat(uniqueDays, nStrats, 1);               % [1;2;3;4; … ]
StrVec   = repelem(strategyList(:), nDays, 1);           % each strategy repeated 4×
withinDesign = table( categorical(DaysVec), categorical(StrVec), ...
    'VariableNames', {'Day','Strategy'});

measureVars   = ratTable.Properties.VariableNames(3:end);  % all <Strat>_Day# cols

% Plot
for age = {'young','old'}
    ag = age{1};
    idxAge = strcmp(ratTable.Age, ag);
    subTbl  = ratTable(idxAge, :);

    % Fit RM model: no between factor, two within factors Day & StrategyGroup
    formula = sprintf('%s-%s ~ 1', measureVars{1}, measureVars{end});
    rm = fitrm(subTbl, ...
        sprintf('%s-%s ~ 1', measureVars{1}, measureVars{end}), ...
        'WithinDesign', withinDesign, ...
        'WithinModel',  'Day*Strategy');   % <-- matches the column names!


    % ANOVA: test main effects & interaction
    anovaResults = ranova(rm, 'WithinModel', 'Day*Strategy');
    writetable(anovaResults, fullfile(processed_dir, sprintf('All_ANOVA_%s.csv', ag)));
    saveTableAsFigure(anovaResults, ...
        fullfile(processed_dir, sprintf('ANOVA_%s.png', ag,...
        'Title',sprintf('ANOVA Results (%s rats)', upper(ag)) )));
    disp(anovaResults)
    % Tukey post hoc for between‑strategy (by Day)
    postHocStrat = multcompare(rm, 'Strategy', 'By', 'Day', 'ComparisonType', 'tukey-kramer');
    writetable(postHocStrat, fullfile(processed_dir, sprintf('ALL_PostHoc_Strat_%s.csv', ag)));
    % % Tukey post hoc for within‑strategy (by StrategyGroup)
    postHocDay   = multcompare(rm, 'Day', 'By', 'Strategy', 'ComparisonType', 'tukey-kramer');
    writetable(postHocDay,   fullfile(processed_dir, sprintf('ALL_PostHoc_Day_%s.csv',   ag)));
    %
    %Plot: Grouped Bar Chart with SEM, Jittered Per‑Rat Means, sigstar
    % Compute group means & SEM across rats
    meanData = zeros(nDays, nStrats);
    semData  = zeros(nDays, nStrats);
    dataMeansPerGroup = cell(nStrats, nDays);
    for g = 1:nStrats
        for d = 1:nDays
            colName = sprintf('%s_Day%d', strategyList{g}, uniqueDays(d));
            vals = subTbl.(colName);
            meanData(d,g) = mean(vals);
            semData(d,g)  = std(vals)/sqrt(numel(vals));
            dataMeansPerGroup{g,d} = vals;
        end
    end

    % Plot
    f = figure('Position',[95,100,1100,630]);
    hold on;
    hBar = bar(uniqueDays, meanData, 'grouped', 'BarWidth', 0.8);
    for g = 1:nStrats
        hBar(g).FaceColor = stratColors(g,:);
        hBar(g).FaceAlpha = 0.35;
    end
    drawnow;
    % Get bar centers
    barCenters = reshape([hBar.XEndPoints], nStrats, nDays)' ;
    %Re‐compute bar centers so they’re aligned
    barCenters = nan(nDays, nStrats);
    for g = 1:nStrats
        barCenters(:,g) = hBar(g).XEndPoints';
    end

    % Jittered per-rat means
    jit = 0.03;
    for g = 1:nStrats
        for d = 1:nDays
            pts = dataMeansPerGroup{g,d};
            xj  = barCenters(d,g) + (rand(size(pts))-0.5)*jit;
            scatter(xj, pts, 20, 'MarkerFaceColor', stratColors(g,:), ...
                'MarkerEdgeColor','none','MarkerFaceAlpha',0.8);
        end
    end
    % Error bars
    for g = 1:nStrats
        errorbar(barCenters(:,g), meanData(:,g), semData(:,g), 'k', 'LineStyle','none','LineWidth',2);
    end

    % Prepare sigstar inputs
    % sigPairs = {};
    % sigP     = [];
    % sigCols  = {};
    %
    % % 1. Between‐strategy comparisons on each day (gray)
    % cs = multcompare(rm, 'Strategy', 'By', 'Day', 'ComparisonType', 'tukey-kramer');
    % % cs is a table: cs.Day, cs.StrategyGroup_1, cs.StrategyGroup_2, cs.pValue, etc.
    % for rr = 1:height(cs)
    %     pval = cs.pValue(rr);
    %     if pval < 0.05
    %         % get day index
    %         dVal = cs.Day(rr);
    %         dIdx = find(uniqueDays == double(dVal), 1);
    %         % find group indices
    %         lvl1 = char(cs.Strategy_1(rr));
    %         lvl2 = char(cs.Strategy_2(rr));
    %         g1 = find(strcmp(strategyList, lvl1), 1);
    %         g2 = find(strcmp(strategyList, lvl2), 1);
    %         % only plot once per unordered pair
    %         if ~isempty(dIdx) && ~isempty(g1) && ~isempty(g2) && (g1 < g2)
    %             sigPairs{end+1} = [barCenters(dIdx, g1), barCenters(dIdx, g2)];
    %             sigP(end+1)     = pval;
    %             sigCols{end+1}  = [0.5,0.5,0.5];  % gray
    %         end
    %     end
    % end
    % if ~isempty(sigPairs)
    %     hSig = sigstar(sigPairs, sigP);
    %     for k = 1:size(hSig,1)
    %         set(hSig(k,1), 'Color', sigCols{k});   % line
    %         set(hSig(k,2), 'Color', sigCols{k});   % stars
    %     end
    % end
    %
    %
    % % % 2. Within‐strategy consecutive‐day comparisons (black)
    % % Recompute bar centers
    % barCenters = nan(nDays, nStrats);
    % for g = 1:nStrats
    %     barCenters(:,g) = hBar(g).XEndPoints';
    % end
    %
    % % Preallocate
    % sigPairs = {};
    % sigP     = [];
    % sigCols  = {};
    % nsAnn    = struct('x',{},'y',{},'text',{},'color',{});
    %
    % % (Re)compute bar centers
    % barCenters = nan(nDays, nStrats);
    % for g = 1:nStrats
    %     barCenters(:,g) = hBar(g).XEndPoints';
    % end
    %
    % % % Prepare containers
    % sigPairs = {};
    % sigP     = [];
    % sigCols  = {};
    %
    % % Gray color for 2‑out‑of‑3 lines
    % darkGray = [0.2 0.2 0.2];
    %
    % % Run Tukey for Day within StrategyGroup
    % cd = multcompare(rm, 'Day','By','StrategyGroup','ComparisonType','tukey-kramer');
    %
    % for d = 1:nDays-1
    %     % select only consecutive‐day rows
    %     sel = double(cd.Day_1)==uniqueDays(d) & double(cd.Day_2)==uniqueDays(d+1);
    %     sub = cd(sel,:);
    %     % build p‑vector
    %     pv = nan(1,nStrats);
    %     for rr = 1:height(sub)
    %         gidx = find(strcmp(strategyList, char(sub.StrategyGroup(rr))),1);
    %         pv(gidx) = sub.pValue(rr);
    %     end
    %     sigIdx = find(pv<0.05);
    %     switch numel(sigIdx)
    %         case 0
    %             % none significant → do nothing
    %         case 3
    %             % all three → black line between avg centers
    %             x1 = mean(barCenters(d, :));
    %             x2 = mean(barCenters(d+1, :));
    %             sigPairs{end+1} = [x1, x2];
    %             sigP(end+1)     = min(pv(sigIdx));
    %             sigCols{end+1}  = [0 0 0];
    %         case 2
    %             % two significant → dark gray line
    %             x1 = mean(barCenters(d, sigIdx));
    %             x2 = mean(barCenters(d+1, sigIdx));
    %             sigPairs{end+1} = [x1, x2];
    %             sigP(end+1)     = min(pv(sigIdx));
    %             sigCols{end+1}  = darkGray;
    %             % annotate the non‑significant with 'ns' in its own color
    %             ns = setdiff(1:nStrats, sigIdx);
    %             yLine = max(max(meanData(d, sigIdx)), max(meanData(d+1, sigIdx))) + 0.02;
    %             x_ns  = barCenters(d, ns);
    %             % text(x_ns+0.5, repmat(yLine+0.55, size(x_ns)), 'n.s.', ...
    %             %     'Color', stratColors(ns,:), 'HorizontalAlignment','Center');
    %         case 1
    %             % single significant → colored line
    %             g = sigIdx;
    %             x1 = barCenters(d, g);
    %             x2 = barCenters(d+1, g);
    %             sigPairs{end+1} = [x1, x2];
    %             sigP(end+1)     = pv(g);
    %             sigCols{end+1}  = stratColors(g,:);
    %     end
    % end
    %
    % % filter out invalid entries
    % valid = cellfun(@(c)isnumeric(c)&&numel(c)==2, sigPairs);
    % sigPairs = sigPairs(valid);
    % sigP     = sigP(valid);
    % sigCols  = sigCols(valid);
    %
    % % draw significance bars & stars
    % if ~isempty(sigPairs)
    %     hSig = sigstar(sigPairs, sigP);
    %     for k = 1:size(hSig,1)
    %         set(hSig(k,1), 'Color', sigCols{k});   % line
    %         set(hSig(k,2), 'Color', sigCols{k});   % stars
    %     end
    % end

    % Final formatting
    title(sprintf('%s', upper(ag)));
    xlabel('Day');
    ylabel('Probability for Strategy Group');
    xticks(uniqueDays);
    ylim([0 0.8]);
    legend(strategyList, 'Location','northeast');

    pubify_figure_axis_robust(16,16);
    hold off;

    % Save figure
    exportgraphics(f, fullfile(fig_dir, sprintf('Day_ALLStrategies_%s.png', ag)), ...
        'Resolution', 450);  % Set to 300 DPI

    close(f);
end
%% 11) Within strategy group ANOVA and comparisons- STATS
strategyGroups = { ...
    {'thigmotaxis','circling','random_path'},'NonGoal'; ...
    {'scanning','chaining'},  'Procedural'; ...
    {'directed_search','corrected_search','direct_path'},  'Allocentric'};

groupPrettyNames = {'Non-Goal Oriented','Procedural','Allocentric'};

% 8 strategy palette
stratColors = [
    0.850 0.325 0.098;
    0.094 0.509 0.800;
    0.467 0.675 0.188;
    0.800 0.475 0.655;
    0.894 0.766 0.039;
    0.436 0.600 0.800;
    0.905 0.161 0.541;
    0.400 0.400 0.400];

% Mkdir and initially  parameters
outDir = fullfile(processed_dir,'Strategy_Day_Comparisons');
if ~exist(outDir,'dir'); mkdir(outDir); end

uniqueDays  = unique(data1.Day);   nDays  = numel(uniqueDays);
rats        = unique(data1.x_TargetID);
nRats       = numel(rats);

%------------   Descriptive Mean ± SD over rats  --------------
descTab = table();   % we’ll concatenate rows

for g = 1:size(strategyGroups,1)
    stratList   = strategyGroups{g,1};          % cellstr of strategies in this subgroup
    subgroupTag = strategyGroups{g,2};          % 'NonGoal' / etc.

    for dIdx = 1:nDays
        dVal = uniqueDays(dIdx);
        for ageTag = ["young","old"]

            % gather all trials of that day *and* age
            dayAgeMask = (data1.Day == dVal) & strcmp(data1.Age, ageTag);

            for sIdx = 1:numel(stratList)
                stratName = stratList{sIdx};
                vals = data1.(stratName)(dayAgeMask);
                m  = mean(vals,'omitnan');
                sd = std(vals,'omitnan');

                descTab = [descTab;       %#ok<AGROW>
                    table(categorical({subgroupTag}), categorical(ageTag), dVal, ...
                    categorical({stratName}), m, sd, ...
                    'VariableNames',{'Group','Age','Day','Strategy','Mean','SD'})];
            end
        end
    end
end
writetable(descTab, fullfile(outDir,'MeanSD_StrategyUse.csv'));

%----------  Loop through each subgroup for RM-ANOVAs ---------
for g = 1:size(strategyGroups,1)
    stratList   = strategyGroups{g,1};
    subgroupTag = strategyGroups{g,2};
    prettyName  = groupPrettyNames{g};
    nStrats     = numel(stratList);

    % ---- Build per-rat matrix  (rows=rats, cols = strat×day) ----
    ratMeans = nan(nRats, nStrats*nDays);
    AgeCell  = cell(nRats,1);

    for r = 1:nRats
        rMask       = strcmp(data1.x_TargetID, rats{r});
        AgeCell{r}  = data1.Age{find(rMask,1)};
        col = 1;
        for sIdx = 1:nStrats
            vec = data1.(stratList{sIdx})(rMask);
            dayVec = data1.Day(rMask);
            for dIdx = 1:nDays
                ratMeans(r,col) = mean(vec(dayVec==uniqueDays(dIdx)));
                col = col + 1;
            end
        end
    end

    % ---- Variable names  ---------------------------
    varNames = {'RatID','Age'};
    for sIdx = 1:nStrats
        for dIdx = 1:nDays
            varNames{end+1} = sprintf('%s_Day%d', stratList{sIdx}, uniqueDays(dIdx));
        end
    end
    T = cell2table([rats, AgeCell, num2cell(ratMeans)], 'VariableNames',varNames);

    % ---- Within-design table ------------------------------------
    DaysVec = repmat(uniqueDays, nStrats,1);
    StrVec  = repelem(stratList(:), nDays,1);
    WD      = table(categorical(DaysVec), categorical(StrVec), ...
        'VariableNames',{'Day','Strategy'});
    measVars = T.Properties.VariableNames(3:end);

    % Fit the repeated‐measures model
    rm = fitrm(T, sprintf('%s-%s ~ Age', measVars{1}, measVars{end}), ...
        'WithinDesign', WD, ...
        'WithinModel',  'Day*Strategy');

    % Get the full ANOVA (including Age×Day×Strategy)
    a = ranova(rm,'WithinModel','Day*Strategy');
    % ---- compute partial η² ----
    effNames = a.Properties.RowNames;   % row names, e.g. '(Intercept):Day'
    SS       = a.SumSq;
    eta      = nan(size(SS));

    for i = 1:numel(effNames)
        eff = effNames{i};

        % Skip interceipt and error rows
        if strcmp(eff,'(Intercept)') || startsWith(eff,'Error')
            continue
        end

        % Pick the right error term by matching the effect name
        switch eff
            case 'Age'
                errName = 'Error';                    % between-subjects
            case {'(Intercept):Day','Age:Day'}
                errName = 'Error(Day)';
            case {'(Intercept):Strategy','Age:Strategy'}
                errName = 'Error(Strategy)';
            case {'(Intercept):Day:Strategy','Age:Day:Strategy'}
                errName = 'Error(Day:Strategy)';
            otherwise
                warning('No mapping for effect “%s” – η² left NaN', eff);
                continue
        end
        % find corresponding row
        idxE = find(strcmp(effNames, errName), 1);
        if isempty(idxE)
            warning('Couldn''t find error row “%s” for effect “%s”', errName, eff);
            continue
        end
        eta(i)   = SS(i) / (SS(i) + ssE);
    end
    a.etaSq = round(eta,2);

    % ---- build a table with Effect, df, SS, MS, F, p, etaSq ----
    tbl = a(:, {'DF','SumSq','MeanSq','F','pValue','etaSq'});
    tbl = addvars(tbl, effNames, 'Before','DF', 'NewVariableNames','Effect');
    tbl.Properties.VariableNames = {'Effect','df','SS','MS','F','p','etaSq'};

    % ---- write to CSV ----
    writetable(tbl, fullfile(outDir, sprintf('RM_ANOVA_%s.csv', subgroupTag)), ...
        'WriteRowNames', false);

    % ---- Post-hoc exports  ----
    phS = multcompare(rm,'Strategy','By','Day','ComparisonType','tukey-kramer');
    phS.sigFlag = phS.pValue < 0.05;
    writetable(phS, fullfile(outDir,sprintf('PostHoc_Strategy_%s.csv',subgroupTag)));

    phD = multcompare(rm,'Day','By','Strategy','ComparisonType','tukey-kramer');
    phD.sigFlag = phD.pValue < 0.05;
    writetable(phD, fullfile(outDir,sprintf('PostHoc_Day_%s.csv',subgroupTag)));

    % ---- Within-age ANOVA exports ----
    for ageTag = ["young","old"]
        idx = strcmp(T.Age, ageTag);
        subT = T(idx,:);
        rmW  = fitrm(subT, sprintf('%s-%s~1',measVars{1},measVars{end}), ...
            'WithinDesign',WD,'WithinModel','Day*Strategy');
        aW   = ranova(rmW,'WithinModel','Day*Strategy');

        % partial η² for within-age
        effNamesW = aW.Properties.RowNames;
        SSW       = aW.SumSq;
        etaW      = nan(size(SSW));
        for i = 1:numel(effNamesW)
            eff = effNamesW{i};
            if strcmp(eff,'(Intercept)') || startsWith(eff,'Error'), continue; end
            parts = regexp(eff,'Day|Strategy','match');
            if isempty(parts), errName = 'Error';
            else errName = ['Error(' strjoin(parts,':') ')']; end
            idxE = find( strcmp(effNamesW,errName) | strcmp(effNamesW,strrep(errName,':','*')),1);
            if ~isempty(idxE)
                etaW(i) = SSW(i) / (SSW(i) + SSW(idxE));
            end
        end
        aW.etaSq = round(etaW,4);

        tblW = aW(:, {'DF','SumSq','MeanSq','F','pValue','etaSq'});
        tblW.Properties.VariableNames = {'df','SS','MS','F','p','etaSq'};
        writetable(tblW, fullfile(outDir, sprintf('WithinAge_ANOVA_%s_%s.csv',subgroupTag,upper(char(ageTag)))), ...
            'WriteRowNames', true);
    end
    % --- Verify with line+errorbar plots per strategy group ---
    markers = {'o','s','^'};  % up to 9 distinct markers

    for g = 1:size(strategyGroups,1)
        stratList   = strategyGroups{g,1};     % e.g. {'thigmotaxis','circling','random_path'}
        subgroupTag = strategyGroups{g,2};     % 'NonGoal', etc.
        nStrats     = numel(stratList);

        figure('Name',subgroupTag,'Color','w'); hold on;
        legendEntries = cell(nStrats*2,1);
        ages = {'young','old'};

        for aIdx = 1:2
            ageTag = ages{aIdx};
            if strcmp(ageTag,'young')
                colC = youngColor;
            else
                colC = oldColor;
            end

            % small age-based shift so young/old lines don’t overlap exactly:
            ageShift = (aIdx - 1.5)*0.08;

            for sIdx = 1:nStrats
                stratName = stratList{sIdx};
                marker    = markers{sIdx};

                % Compute mean & SEM over rats for this age/day/strategy
                meanVals = nan(size(uniqueDays));
                semVals  = nan(size(uniqueDays));
                for d = 1:numel(uniqueDays)
                    mask = (data1.Day == uniqueDays(d)) & strcmp(data1.Age, ageTag);
                    v    = data1.(stratName)(mask);
                    meanVals(d) = mean(v,'omitnan');
                    semVals(d)  = std(v,'omitnan')/sqrt(sum(mask));
                end

                % horizontal shift by strategy so lines fan out
                stratShift = (sIdx - (nStrats+1)/2)*0.15;
                xPos = uniqueDays + ageShift + stratShift;

                % Plot errorbars + line + marker
                h = errorbar(xPos, meanVals, semVals, '-o', ...
                    'Color',       colC, ...
                    'Marker',      marker, ...
                    'MarkerFaceColor', colC, ...
                    'MarkerSize',  6, ...
                    'LineWidth',   1.5, ...
                    'DisplayName', sprintf('%s (%s)', stratName, ageTag));

                % Store legend entry automatically via DisplayName
            end
        end

        % Final formatting
        xlabel('Day');
        ylabel('Mean Strategy Probability');
        title(sprintf('%s: Strategy × Day by Age', subgroupTag), 'Interpreter','none');
        xticks(uniqueDays);
        legend('Location','bestoutside');
        hold off;
    end


end
