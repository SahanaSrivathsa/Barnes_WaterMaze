%% ------------- MUST RUN SECTION FOR ALL FURTHER ANALYSIS ------------- %%
% Dirs
base_dir='G:\WMAZE_Data\Data_Behaviour\Data_CognitiveBattery\WaterMaze\Mia_Proj';
processed_dir=fullfile(base_dir,'StrategyProcessed');
fig_dir=fullfile(base_dir,'Figures');
fig_trans_dir=fullfile(fig_dir,'Transition');
if~exist(fig_trans_dir,'dir')
    mkdir(fig_trans_dir)
end

% Define the strategy column names (adjust names if needed)
strategyNames = {'thigmotaxis','circling','random_path','scanning',...
    'chaining','directed_search','corrected_search','direct_path','perseverance'};
nStrategies = numel(strategyNames);

strategy_titles={'Thigmotaxis','Circling','Random Path','Scanning',...
    'Chaining','Directed Search','Corrected Search','Direct Path','Perseverance'};
%Load data from Excel files
data1 = readtable(fullfile(base_dir,'MWM_results.xlsx'));  % From RTrack, contains Track_ID, Strategy, Age
%% CHECK FOR CONSISTENCY in datasets (Trial no and full values)

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


%% Transition Probabilities: Define Strategies and Separate Data by Age
strategyNames = {'thigmotaxis','circling','random_path','scanning',...
    'chaining','directed_search','corrected_search','direct_path','perseverance'};
numStrategies = length(strategyNames);

% Separate data by Age group
youngData = data1(strcmpi(data1.Age, 'young'), :);
oldData   = data1(strcmpi(data1.Age, 'old'), :);

% Get unique rat IDs for each group
youngRatIDs = unique(youngData.x_TargetID);
oldRatIDs   = unique(oldData.x_TargetID);

% Assume days are numeric; determine total number of days
totalDays = max(data1.Day);

%% Compute Daily Within-Day Transitions for Young Rats
% Two cell arrays: one for discrete transitions, one for probabilistic transitions.
youngDailyDiscreteTrans = cell(length(youngRatIDs), totalDays);
youngDailyProbTrans     = cell(length(youngRatIDs), totalDays);

for i = 1:length(youngRatIDs)
    ratID = youngRatIDs(i);
    % Select and sort data for this rat by Day and Trial order
    ratData = youngData(strcmp(youngData.x_TargetID, ratID), :);
    ratData = sortrows(ratData, {'Day', 'x_Trial'});
    
    % Loop over each day (from 1 to totalDays)
    for d = 1:totalDays
        dayData = ratData(ratData.Day == d, :);
        if height(dayData) > 1
            discMat = zeros(numStrategies, numStrategies);
            probMat = zeros(numStrategies, numStrategies);
            % Loop over consecutive trials for this day
            for t = 1:height(dayData)-1
                currState = dayData.strategy(t);
                nextState = dayData.strategy(t+1);
                % DISCRETE: count the transition from currState to nextState
                discMat(currState, nextState) = discMat(currState, nextState) + 1;
                % PROBABILISTIC: use the probability vectors for a weighted update
                currProbs = table2array(dayData(t, strategyNames));
                nextProbs = table2array(dayData(t+1, strategyNames));
                probMat = probMat + (currProbs' * nextProbs);
            end
            % Normalize each row so that the sum equals 1 (if possible)
            for j = 1:numStrategies
                if sum(discMat(j,:)) > 0
                    discMat(j,:) = discMat(j,:) / sum(discMat(j,:));
                end
                if sum(probMat(j,:)) > 0
                    probMat(j,:) = probMat(j,:) / sum(probMat(j,:));
                end
            end
            youngDailyDiscreteTrans{i, d} = discMat;
            youngDailyProbTrans{i, d}     = probMat;
        else
            % If not enough trials, store NaN matrix
            youngDailyDiscreteTrans{i, d} = NaN(numStrategies, numStrategies);
            youngDailyProbTrans{i, d}     = NaN(numStrategies, numStrategies);
        end
    end
end

%Compute Daily Within-Day Transitions for Old Rats
oldDailyDiscreteTrans = cell(length(oldRatIDs), totalDays);
oldDailyProbTrans     = cell(length(oldRatIDs), totalDays);

for i = 1:length(oldRatIDs)
    ratID = oldRatIDs(i);
    ratData = oldData(strcmp(oldData.x_TargetID, ratID), :);
    ratData = sortrows(ratData, {'Day', 'x_Trial'});
    for d = 1:totalDays
        dayData = ratData(ratData.Day == d, :);
        if height(dayData) > 1
            discMat = zeros(numStrategies, numStrategies);
            probMat = zeros(numStrategies, numStrategies);
            for t = 1:height(dayData)-1
                currState = dayData.strategy(t);
                nextState = dayData.strategy(t+1);
                discMat(currState, nextState) = discMat(currState, nextState) + 1;
                currProbs = table2array(dayData(t, strategyNames));
                nextProbs = table2array(dayData(t+1, strategyNames));
                probMat = probMat + (currProbs' * nextProbs);
            end
            for j = 1:numStrategies
                if sum(discMat(j,:)) > 0
                    discMat(j,:) = discMat(j,:) / sum(discMat(j,:));
                end
                if sum(probMat(j,:)) > 0
                    probMat(j,:) = probMat(j,:) / sum(probMat(j,:));
                end
            end
            oldDailyDiscreteTrans{i, d} = discMat;
            oldDailyProbTrans{i, d}     = probMat;
        else
            oldDailyDiscreteTrans{i, d} = NaN(numStrategies, numStrategies);
            oldDailyProbTrans{i, d}     = NaN(numStrategies, numStrategies);
        end
    end
end

%% Aggregate Daily Transition Matrices (Averaging Across Rats)
% Initialize aggregated cell arrays for both groups
youngAggDiscrete = cell(totalDays, 1);
oldAggDiscrete   = cell(totalDays, 1);
youngAggProb     = cell(totalDays, 1);
oldAggProb       = cell(totalDays, 1);

for d = 1:totalDays
    %Aggregate Young Group
    count = 0;
    aggDisc = zeros(numStrategies, numStrategies);
    for i = 1:length(youngRatIDs)
         M = youngDailyDiscreteTrans{i, d};
         if ~any(isnan(M(:)))
             aggDisc = aggDisc + M;
             count = count + 1;
         end
    end
    if count > 0, aggDisc = aggDisc / count; end
    youngAggDiscrete{d} = aggDisc;
    
    count = 0;
    aggProb = zeros(numStrategies, numStrategies);
    for i = 1:length(youngRatIDs)
         M = youngDailyProbTrans{i, d};
         if ~any(isnan(M(:)))
             aggProb = aggProb + M;
             count = count + 1;
         end
    end
    if count > 0, aggProb = aggProb / count; end
    youngAggProb{d} = aggProb;
    
    %Aggregate Old Group
    count = 0;
    aggDisc = zeros(numStrategies, numStrategies);
    for i = 1:length(oldRatIDs)
         M = oldDailyDiscreteTrans{i, d};
         if ~any(isnan(M(:)))
             aggDisc = aggDisc + M;
             count = count + 1;
         end
    end
    if count > 0, aggDisc = aggDisc / count; end
    oldAggDiscrete{d} = aggDisc;
    
    count = 0;
    aggProb = zeros(numStrategies, numStrategies);
    for i = 1:length(oldRatIDs)
         M = oldDailyProbTrans{i, d};
         if ~any(isnan(M(:)))
             aggProb = aggProb + M;
             count = count + 1;
         end
    end
    if count > 0, aggProb = aggProb / count; end
    oldAggProb{d} = aggProb;
end

%% Repeated Measures ANOVA and Tukey Post Hoc for a Specific Transition
% CURRENTLY ONLY ONE STRATEGY 
strategyIdx = 1;  

% Create matrices to store daily self-transition values for each rat.
youngSelfDaily = NaN(length(youngRatIDs), totalDays);
for i = 1:length(youngRatIDs)
    for d = 1:totalDays
         M = youngDailyDiscreteTrans{i, d};
         if ~any(isnan(M(:)))
             youngSelfDaily(i, d) = M(strategyIdx, strategyIdx);
         end
    end
end

oldSelfDaily = NaN(length(oldRatIDs), totalDays);
for i = 1:length(oldRatIDs)
    for d = 1:totalDays
         M = oldDailyDiscreteTrans{i, d};
         if ~any(isnan(M(:)))
             oldSelfDaily(i, d) = M(strategyIdx, strategyIdx);
         end
    end
end

% Combine into one table for repeated measures analysis.
allSelfDaily = [youngSelfDaily; oldSelfDaily];
Age = [repmat({'young'}, size(youngSelfDaily,1), 1); repmat({'old'}, size(oldSelfDaily,1), 1)];
T = array2table(allSelfDaily, 'VariableNames', strcat('Day', string(1:totalDays)));
T.Age = Age;

% Run the mixed-design repeated measures ANOVA and Tukey post hoc tests.
anovaResults = runMixedANOVA(T, strcat('Day', string(1:totalDays)), 'Day');
postHocResults = runTukeyPostHocMixed(T, strcat('Day', string(1:totalDays)));


%% Plotting: Transition Probabilities as Heat Maps
% Example: Plot aggregated heat maps for young rats 
f1=figure('Position', [100, 100, 2200, 1200]);
for d = 1:totalDays
    % Discrete transitions heat map
    subplot(2, totalDays, d);
    h1 = heatmap(strategy_titles, strategy_titles, youngAggDiscrete{d});
    h1.Title = ['Young: Day' num2str(d)];
    h1.XLabel = 'Next Strategy';
    h1.YLabel = 'Current Strategy';
    clim([0 0.5])
    %Old
    subplot(2, totalDays, d+totalDays);
    h1 = heatmap(strategy_titles, strategy_titles, oldAggDiscrete{d});
    h1.Title = ['Old: Day' num2str(d)];
    h1.XLabel = 'Next Strategy';
    h1.YLabel = 'Current Strategy';
    clim([0 0.5])
   
end
sgtitle('Discrete Strategy Transition Probabilities')
saveas(f1,fullfile(fig_trans_dir,'DiscreteTransitionProbabilities'),'png');
close;
 % Probabilistic transitions heat map
f2=figure('Position', [100, 100, 2200, 1200]);
for d = 1:totalDays

    subplot(2, totalDays, d);
    h1 = heatmap(strategy_titles, strategy_titles, youngAggProb{d});
    h1.Title = ['Young: Day' num2str(d)];
    h1.XLabel = 'Next Trial';
    h1.YLabel = 'Current Trial';
    clim([0 0.4])
    subplot(2, totalDays, d+totalDays);
    h2 = heatmap(strategy_titles, strategy_titles, oldAggProb{d});
    h2.Title = ['Old: Day' num2str(d)];
    h2.XLabel = 'Next Trial';
    h2.YLabel = 'Current Trial';
    clim([0 0.4])
end
sgtitle('Strategy Confidence Intervals in a Transition Matrix ')
saveas(f2,fullfile(fig_trans_dir,'ConfidenceIntervalTransitionProbabilities'),'png');
close;


%% Plotting: Transitions Matrix - Network Graph
% Get unique days from the data
uniqueDays = unique(data1.Day);
nDays = length(uniqueDays);

figure('Position', [100, 100, 1400, 200*nDays]);

for d = 1:nDays
    day = uniqueDays(d);
    
    % --- Young rats for day 'day' ---
    dayYoungData = youngData(youngData.Day == day, :);
    youngDayTrans = zeros(numStrategies, numStrategies);
    for i = 1:length(youngRatIDs)
        ratData = dayYoungData(strcmp(dayYoungData.x_TargetID, youngRatIDs{i}), :);
        ratData = sortrows(ratData, 'x_Trial');
        for t = 1:height(ratData)-1
            currStrat = ratData.strategy(t);
            nextStrat = ratData.strategy(t+1);
            youngDayTrans(currStrat, nextStrat) = youngDayTrans(currStrat, nextStrat) + 1;
        end
    end
    % Normalize row-wise
    for i = 1:numStrategies
        if sum(youngDayTrans(i,:)) > 0
            youngDayTrans(i,:) = youngDayTrans(i,:) / sum(youngDayTrans(i,:));
        end
    end
    
    % --- Old rats for day 'day' ---
    dayOldData = oldData(oldData.Day == day, :);
    oldDayTrans = zeros(numStrategies, numStrategies);
    for i = 1:length(oldRatIDs)
        ratData = dayOldData(strcmp(dayOldData.x_TargetID, oldRatIDs{i}), :);
        ratData = sortrows(ratData, 'x_Trial');
        for t = 1:height(ratData)-1
            currStrat = ratData.strategy(t);
            nextStrat = ratData.strategy(t+1);
            oldDayTrans(currStrat, nextStrat) = oldDayTrans(currStrat, nextStrat) + 1;
        end
    end
    for i = 1:numStrategies
        if sum(oldDayTrans(i,:)) > 0
            oldDayTrans(i,:) = oldDayTrans(i,:) / sum(oldDayTrans(i,:));
        end
    end
    
    % --- Plotting the network graphs ---
    threshold = 0.05;
    
    % Young subplot for the day
    subplot(nDays, 2, (d-1)*2 + 1);
    G_young_day = digraph(youngDayTrans, strategy_titles);
    % Remove edges below or equal to threshold
    removeIdx = find(G_young_day.Edges.Weight <= threshold);
    G_young_day = rmedge(G_young_day, removeIdx);
    if ~isempty(G_young_day.Edges)
        LWidths = 1 + 5 * G_young_day.Edges.Weight / max(G_young_day.Edges.Weight);
        edgeLabels = round(G_young_day.Edges.Weight * 100) / 100;
        p = plot(G_young_day, 'LineWidth', LWidths);
        p.NodeColor = youngColor;
        p.NodeFontSize = 10;
        p.EdgeColor = youngColor;
        p.EdgeAlpha = 0.7;
        p.EdgeFontSize = 8;
        layout(p, 'circle');
    else
        text(0.5, 0.5, 'No transitions > threshold', 'HorizontalAlignment', 'center');
        axis([0 1 0 1]);
    end
    title(['Day ' num2str(day) ' - Young']);
    
    % Old subplot for the day
    subplot(nDays, 2, (d-1)*2 + 2);
    G_old_day = digraph(oldDayTrans, strategy_titles);
    removeIdx = find(G_old_day.Edges.Weight <= threshold);
    G_old_day = rmedge(G_old_day, removeIdx);
    if ~isempty(G_old_day.Edges)
        LWidths = 1 + 5 * G_old_day.Edges.Weight / max(G_old_day.Edges.Weight);
        edgeLabels = round(G_old_day.Edges.Weight * 100) / 100;
        p = plot(G_old_day,  'LineWidth', LWidths);
        p.NodeColor = oldColor;
        p.NodeFontSize = 10;
        p.EdgeColor = oldColor;
        p.EdgeAlpha = 0.7;
        p.EdgeFontSize = 8;
        layout(p, 'circle');
    else
        text(0.5, 0.5, 'No transitions > threshold', 'HorizontalAlignment', 'center');
        axis([0 1 0 1]);
    end
    title(['Day ' num2str(day) ' - Old']);
end

%% Self Transition Probabilities - Discrete Strategy
%Computing 
% Compute aggregated self-transition probabilities (averaged across days) per rat for each strategy.
youngSelfAgg = zeros(length(youngRatIDs), numStrategies);
for i = 1:length(youngRatIDs)
    for s = 1:numStrategies
        values = [];  % collect daily self-transition values for strategy s
        for d = 1:totalDays
            M = youngDailyDiscreteTrans{i, d};
            if ~any(isnan(M(:)))
                values = [values; M(s, s)];
            end
        end
        % Compute the mean self-transition probability for this rat and strategy.
        if ~isempty(values)
            youngSelfAgg(i, s) = mean(values);
        else
            youngSelfAgg(i, s) = NaN;
        end
    end
end

oldSelfAgg = zeros(length(oldRatIDs), numStrategies);
for i = 1:length(oldRatIDs)
    for s = 1:numStrategies
        values = [];
        for d = 1:totalDays
            M = oldDailyDiscreteTrans{i, d};
            if ~any(isnan(M(:)))
                values = [values; M(s, s)];
            end
        end
        if ~isempty(values)
            oldSelfAgg(i, s) = mean(values);
        else
            oldSelfAgg(i, s) = NaN;
        end
    end
end
% Plotting
f3=figure('Position', [100, 100, 1200, 800]);
for s = 1:numStrategies
    subplot(3,3,s);
    hold on;
    
    % Combine data for boxplot for current strategy
    data_box = [youngSelfAgg(:, s); oldSelfAgg(:, s)];
    groups = [ones(size(youngSelfAgg(:, s))); 2*ones(size(oldSelfAgg(:, s)))];
    
    % Draw the boxplot (suppressing outlier symbols)
    boxplot(data_box, groups, 'Colors', [youngColor; oldColor], 'Symbol', '');
    
    % Overlay jittered scatter points for each group
    jitter = 0.05;
    scatter(ones(size(youngSelfAgg(:, s))) + jitter*randn(size(youngSelfAgg(:, s))), youngSelfAgg(:, s), ...
        'filled', 'MarkerFaceColor', youngColor);
    scatter(2 + jitter*randn(size(oldSelfAgg(:, s))), oldSelfAgg(:, s), ...
        'filled', 'MarkerFaceColor', oldColor);
    
    % Perform a two-sample t-test between groups for this strategy
    [~, pVal] = ttest2(youngSelfAgg(:, s), oldSelfAgg(:, s));
    
    % Annotate significance using your sigstar_OnlySig function for groups [1, 2]
    sigstar_OnlySig({[1,2]}, pVal);
    
    % Customize the plot appearance
    set(gca, 'XTick', [1 2], 'XTickLabel', {'Young','Old'});
    title(strategy_titles{s});
    ylabel('Probability');
    xlim([0.5 2.5]);
    ylim([0 1.2]);
    hold off;
end
sgtitle('Self-Transition Probabilities per Strategy');
pubify_figure_axis_robust(14,14)
saveas(f3,fullfile(fig_trans_dir,'SelfTransitionProbabilities'),'png');
close;
%% Compute aggregated self-transition probabilities
% using probabilistic transitions for each rat and each strategy.
youngSelfAggProb = zeros(length(youngRatIDs), numStrategies);
for i = 1:length(youngRatIDs)
    for s = 1:numStrategies
        values = [];  % collect daily self-transition values for strategy s
        for d = 1:totalDays
            M = youngDailyProbTrans{i, d};
            if ~any(isnan(M(:)))
                values = [values; M(s, s)];
            end
        end
        if ~isempty(values)
            youngSelfAggProb(i, s) = mean(values);
        else
            youngSelfAggProb(i, s) = NaN;
        end
    end
end

oldSelfAggProb = zeros(length(oldRatIDs), numStrategies);
for i = 1:length(oldRatIDs)
    for s = 1:numStrategies
        values = [];
        for d = 1:totalDays
            M = oldDailyProbTrans{i, d};
            if ~any(isnan(M(:)))
                values = [values; M(s, s)];
            end
        end
        if ~isempty(values)
            oldSelfAggProb(i, s) = mean(values);
        else
            oldSelfAggProb(i, s) = NaN;
        end
    end
end

% Plot the aggregated probabilistic self-transition probabilities for each strategy
fProb = figure('Position', [100, 100, 1200, 800]);
for s = 1:numStrategies
    subplot(3,3,s);
    hold on;
    
    % Combine data for boxplot for current strategy
    data_box = [youngSelfAggProb(:, s); oldSelfAggProb(:, s)];
    groups   = [ones(size(youngSelfAggProb(:, s))); 2*ones(size(oldSelfAggProb(:, s)))];
    
    % Draw the boxplot (suppressing outlier symbols)
    boxplot(data_box, groups, 'Colors', [youngColor; oldColor], 'Symbol', '');
    
    % Overlay jittered scatter points for each group
    jitter = 0.05;
    scatter(ones(size(youngSelfAggProb(:, s))) + jitter*randn(size(youngSelfAggProb(:, s))), youngSelfAggProb(:, s), ...
        'filled', 'MarkerFaceColor', youngColor);
    scatter(2 + jitter*randn(size(oldSelfAggProb(:, s))), oldSelfAggProb(:, s), ...
        'filled', 'MarkerFaceColor', oldColor);
    
    % Perform a two-sample t-test between groups for the current strategy
    [~, pVal] = ttest2(youngSelfAggProb(:, s), oldSelfAggProb(:, s));
    
    % Annotate significance using sigstar_OnlySig for groups [1,2]
    sigstar_OnlySig({[1,2]}, pVal);
    
    % Customize plot appearance
    set(gca, 'XTick', [1 2], 'XTickLabel', {'Young','Old'});
    title(strategy_titles{s});
    ylabel('Self-Transition Probability');
    xlim([0.5 2.5]);
    ylim([0 1.2]);
    hold off;
end
sgtitle(' Self-Transition Probabilities: Confidence Intervals');
pubify_figure_axis_robust(14,14)
saveas(fProb, fullfile(fig_trans_dir, 'SelfTransitionProbabilities_CI '), 'png');
close;


%% Plotting: Time Series Change in Strategy Use
% Here we assume that you want to plot, for each age group, how often each strategy is used 
% on each day. We can compute the average use of a strategy per day as follows.
% (For demonstration, we'll compute the average probability that the strategy was "active"
% on a trial using the aggregated discrete matrices.)

% Initialize matrices to store daily average strategy use (frequency)
avgStrategyYoung = zeros(totalDays, numStrategies);
avgStrategyOld   = zeros(totalDays, numStrategies);

for d = 1:totalDays
    % For young, average across rats (assume valid matrices)
    values = [];
    for i = 1:length(youngRatIDs)
        M = youngDailyDiscreteTrans{i, d};
        if ~any(isnan(M(:)))
            % For each rat, use the row sum (or alternatively, the diagonal or another measure)
            % Here we compute the average usage probability per strategy as the mean across rows.
            values = [values; M];
        end
    end
    if ~isempty(values)
        % Average across all rat matrices: here, compute mean along first dimension
        avgStrategyYoung(d, :) = mean(values, 1);
    else
        avgStrategyYoung(d, :) = NaN(1, numStrategies);
    end
    
    % For old group
    values = [];
    for i = 1:length(oldRatIDs)
        M = oldDailyDiscreteTrans{i, d};
        if ~any(isnan(M(:)))
            values = [values; M];
        end
    end
    if ~isempty(values)
        avgStrategyOld(d, :) = mean(values, 1);
    else
        avgStrategyOld(d, :) = NaN(1, numStrategies);
    end
end

%Plot of strategy use for each age group across days - EXTRA Average 
figure('Position', [100, 100, 1200, 600]);
% Young Group
subplot(2,1,1);
for s = 1:numStrategies
    plot(1:totalDays, avgStrategyYoung(:, s), '-o', 'LineWidth', 2);
    hold on;
end
xlabel('Day');
ylabel('Average Transition Probability');
title('Young: Strategy Use across days');
legend(strategy_titles, 'Location', 'best');
grid on;

% Old Group
subplot(2,1,2);
for s = 1:numStrategies
    plot(1:totalDays, avgStrategyOld(:, s), '-o', 'LineWidth', 2);
    hold on;
end
xlabel('Day');
ylabel('Average Transition Probability');
title('Old: Time Series of Strategy Use');
legend(strategy_titles, 'Location', 'best');
grid on;
