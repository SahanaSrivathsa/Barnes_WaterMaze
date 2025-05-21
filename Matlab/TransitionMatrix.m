%% ------------- MUST RUN SECTION FOR ALL FURTHER ANALYSIS ------------- %%
% Dirs

base_dir='C:\DATA\WaterMaze\Mia_Proj';
processed_dir=fullfile(base_dir,'StrategyProcessed');
fig_dir=fullfile(base_dir,'Figures');
fig_trans_dir=fullfile(fig_dir,'Transition');
if~exist(fig_trans_dir,'dir')
    mkdir(fig_trans_dir)
end

% Define the strategy column names (adjust names if needed)
strategyNames = {'thigmotaxis','circling','random_path','scanning',...
    'chaining','directed_search','corrected_search','direct_path'};
nStrategies = numel(strategyNames);

strategy_titles={'Thigmotaxis','Circling','Random Path','Scanning',...
    'Chaining','Directed Search','Corrected Search','Direct Path'};
%Load data from Excel files
data1 = readtable(fullfile(base_dir,'MWM_results.xlsx'));  % From RTrack, contains Track_ID, Strategy, Age
%CHECK FOR CONSISTENCY in datasets (Trial no and full values)

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


%% COLORS
oldColor=[0.4157,0.1059,0.6039]; %green
 youngColor=[0.2196,0.5569,0.2353];% purple
%% Transition Probabilities: Define Strategies and Separate Data by Age
numStrategies = length(strategyNames);

% Separate data by Age group
youngData = data1(strcmpi(data1.Age, 'young'), :);
oldData   = data1(strcmpi(data1.Age, 'old'), :);

% Get unique rat IDs 
youngRatIDs = unique(youngData.x_TargetID);
oldRatIDs   = unique(oldData.x_TargetID);
totalDays = max(data1.Day);

%% Compute Daily Within-Day Transitions for  Rats
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

%% Avg Transition matrices by rat across days
% Initialize aggregated cell arrays for both groups
youngAggDiscrete = cell(totalDays, 1);
oldAggDiscrete   = cell(totalDays, 1);
youngAggProb     = cell(totalDays, 1);
oldAggProb       = cell(totalDays, 1);
% Old method - did not normalize all rows to 1
% for d = 1:totalDays
%     %Aggregate Young Group
%     count = 0;
%     aggDisc = zeros(numStrategies, numStrategies);
%     for i = 1:length(youngRatIDs)
%          M = youngDailyDiscreteTrans{i, d};
%          if ~any(isnan(M(:)))
%              aggDisc = aggDisc + M;
%              count = count + 1;
%          end
%     end
%     if count > 0, aggDisc = aggDisc / count; end
%     youngAggDiscrete{d} = aggDisc;
% 
%     count = 0;
%     aggProb = zeros(numStrategies, numStrategies);
%     for i = 1:length(youngRatIDs)
%          M = youngDailyProbTrans{i, d};
%          if ~any(isnan(M(:)))
%              aggProb = aggProb + M;
%              count = count + 1;
%          end
%     end
%     if count > 0, aggProb = aggProb / count; end
%     youngAggProb{d} = aggProb;
% 
%     %Aggregate Old Group
%     count = 0;
%     aggDisc = zeros(numStrategies, numStrategies);
%     for i = 1:length(oldRatIDs)
%          M = oldDailyDiscreteTrans{i, d};
%          if ~any(isnan(M(:)))
%              aggDisc = aggDisc + M;
%              count = count + 1;
%          end
%     end
%     if count > 0, aggDisc = aggDisc / count; end
%     oldAggDiscrete{d} = aggDisc;
% 
%     count = 0;
%     aggProb = zeros(numStrategies, numStrategies);
%     for i = 1:length(oldRatIDs)
%          M = oldDailyProbTrans{i, d};
%          if ~any(isnan(M(:)))
%              aggProb = aggProb + M;
%              count = count + 1;
%          end
%     end
%     if count > 0, aggProb = aggProb / count; end
%     oldAggProb{d} = aggProb;
% end

% New method only for discrete strategies

for d = 1:totalDays
    % --- Young aggregate ---
    % Build a 3D array [strategy x strategy x rat]
    YM = nan(numStrategies, numStrategies, nYoung);
    for i = 1:nYoung
        YM(:,:,i) = youngDailyDiscreteTrans{i, d};
    end
    % Element‐wise mean across rats, ignoring NaNs
    A = nanmean(YM, 3);
    % Optional: re‐normalize each row to sum to 1
    rowSums = sum(A, 2);
    A = bsxfun(@rdivide, A, rowSums + eps);
    youngAggDiscrete{d} = A;

    % --- Old aggregate ---
    OM = nan(numStrategies, numStrategies, nOld);
    for i = 1:nOld
        OM(:,:,i) = oldDailyDiscreteTrans{i, d};
    end
    B = nanmean(OM, 3);
    rowSums = sum(B, 2);
    B = bsxfun(@rdivide, B, rowSums + eps);
    oldAggDiscrete{d} = B;
end
%% KL Divergence calculation
% Compute KL on discrete transitions ——— %%
epsKL = 1e-12;
klDiv = @(P,Q) nansum( P(:) .* log2( (P(:)+epsKL) ./ (Q(:)+epsKL) ) );

nInt      = totalDays - 1;
nYoung    = numel(youngRatIDs);
nOld      = numel(oldRatIDs);

KL_youngD = NaN(nYoung, nInt);
KL_oldD   = NaN(nOld,   nInt);

% Young discrete
for i = 1:nYoung
    for d = 1:nInt
        P = youngDailyDiscreteTrans{i, d};
        Q = youngDailyDiscreteTrans{i, d+1};
        if all(~isnan(P(:))) && all(~isnan(Q(:)))
            KL_youngD(i, d) = klDiv(P, Q);
        end
    end
end

% Old discrete
for i = 1:nOld
    for d = 1:nInt
        P = oldDailyDiscreteTrans{i, d};
        Q = oldDailyDiscreteTrans{i, d+1};
        if all(~isnan(P(:))) && all(~isnan(Q(:)))
            KL_oldD(i, d) = klDiv(P, Q);
        end
    end
end

% Build table for mixed‐design ANOVA ——— %%
% Combine into one matrix and add Age as between‐subject factor
allKL_D = [KL_youngD; KL_oldD];
Age      = [ repmat({'young'}, nYoung, 1)
             repmat({'old'},   nOld,   1) ];

varNames = strcat('Int', string(1:nInt));  % {'Int1','Int2',…}

TklD = array2table(allKL_D, 'VariableNames', varNames);
TklD.Age = Age;

%Run the mixed‐design ANOVA and Tukey post‐hoc ——— %%
anovaKL_D   = runMixedANOVA(   TklD, varNames, 'Interval');
disp(anovaKL_D)
postHocKL_D = runTukeyPostHocMixed(TklD, varNames);



%% Shannon Entropy
uniqueDays = sort( unique(data1.Day) );
nDays      = numel(uniqueDays);
nYoung     = numel(youngRatIDs);
nOld       = numel(oldRatIDs);

H_young = NaN(nYoung, nDays);
H_old   = NaN(nOld,   nDays);

% ---- Young
for r = 1:nYoung
    for dIdx = 1:nDays
        day = uniqueDays(dIdx);
        P   = youngDailyDiscreteTrans{r, day};
        if all(~isnan(P(:)))
            H_young(r,dIdx) = matrixEntropy(P);
        end
    end
end

% ---- Old
for r = 1:nOld
    for dIdx = 1:nDays
        day = uniqueDays(dIdx);
        P   = oldDailyDiscreteTrans{r, day};
        if all(~isnan(P(:)))
            H_old(r,dIdx) = matrixEntropy(P);
        end
    end
end

%-Bonferroni Day-to-day comparisons (all rats pooled) 
allH   = [H_young; H_old];
pairs  = nchoosek(1:nDays,2);
rawP   = zeros(size(pairs,1),1);
for k = 1:numel(rawP)
    [~, rawP(k)] = ttest( allH(:,pairs(k,1)), allH(:,pairs(k,2)) );
end
bonfP = min(rawP * numel(rawP), 1);

fprintf('\nPairwise day comparisons (Bonferroni corrected)\n');
for k = 1:numel(rawP)
    fprintf(' Day %d vs %d : p = %.4f\n', ...
        uniqueDays(pairs(k,1)), uniqueDays(pairs(k,2)), bonfP(k));
end

% Bonferroni Young vs Old for each day
pAge_raw  = nan(1,nDays);
for dIdx = 1:nDays
    [~, pAge_raw(dIdx)] = ttest2( H_young(:,dIdx), H_old(:,dIdx), 'Vartype','unequal' );
end
pAge_bonf = min(pAge_raw * nDays, 1);

fprintf('\nYoung vs Old per day (Bonferroni corrected)\n');
for dIdx = 1:nDays
    fprintf(' Day %d : p = %.4f\n', uniqueDays(dIdx), pAge_bonf(dIdx));
end

%-Mixed-design ANOVA
varNames = strcat('Day', string(uniqueDays));
T_H      = array2table([H_young; H_old], 'VariableNames', varNames);
T_H.Age  = [repmat({'young'},nYoung,1); repmat({'old'},nOld,1)];

anovaEntropy   = runMixedANOVA(T_H, varNames, 'Day');
postHocEntropy = runTukeyPostHocMixed(T_H, varNames);

disp(anovaEntropy)

%--------------4) Plot means ± SEM
meanY = nanmean(H_young,1);                          % dim=1
semY  = nanstd (H_young,0,1) ./ sqrt(sum(~isnan(H_young),1));
meanO = nanmean(H_old,1);
semO  = nanstd (H_old,0,1) ./ sqrt(sum(~isnan(H_old),1));

f1=figure;
hold on
errorbar(1:nDays, meanY, semY, 'o-','LineWidth',4,'Color',youngColor);
errorbar(1:nDays, meanO, semO, 'o-','LineWidth',4,'Color',oldColor);
xticks(1:nDays); xticklabels(string(uniqueDays));
xlabel('Day'); ylabel('Shannon entropy');
legend('Young','Old','Location','northeast');
xlim([0.9 4.1])
ylim([1.5 2.6])
%title('Average Entropy of Discrete-Transition Matrices');
pubify_figure_axis_robust(16,16)
exportgraphics(f1, fullfile(fig_trans_dir, 'Per_RatTransitionEntropy.png'), ...
     'Resolution', 450); 
%----------------- Plot all rats 
% f2=figure('Position',[100 100 800 400]); hold on
% x = 1:nDays;
% 
% % --- Young rats (faint lines) ---
% for r = 1:nYoung
%     plot( x, H_young(r,:), '-', 'Color',[youngColor 0.25], 'LineWidth',1 );
% end
% % mean line
% plot( x, meanY, '-o', 'Color',youngColor, 'LineWidth',2, 'MarkerFaceColor',youngColor );
% 
% % --- Old rats ---
% for r = 1:nOld
%     plot( x, H_old(r,:), '-', 'Color',[oldColor 0.25], 'LineWidth',1 );
% end
% plot( x, meanO, '-s', 'Color',oldColor, 'LineWidth',2, 'MarkerFaceColor',oldColor );
% 
% xticks(x); xticklabels(string(uniqueDays));
% xlabel('Day'); ylabel('Shannon entropy (bits)');
% legend({'Young rats','Mean Young','Old rats','Mean Old'},'Location','best');
% title('Per-rat entropy trajectories with group means');
% pubify_figure_axis_robust(14,14);
% exportgraphics(f2, fullfile(fig_trans_dir, sprintf('RatTransitionEntropy_%s.png', ag)), ...
%     'Resolution', 450);  % Set to 300 DPI
% %------------------- Bar Plots
% f3=figure('Position',[150 150 800 400]); hold on
% 
% barWidth = 0.35;
% offset   = 0.17;         % horizontal offset for young vs old points
% for dIdx = 1:nDays
%     % ---- bars ----
%     bar(dIdx-offset, meanY(dIdx), barWidth, 'FaceColor', youngColor, 'FaceAlpha',0.3);
%     bar(dIdx+offset, meanO(dIdx), barWidth, 'FaceColor', oldColor,   'FaceAlpha',0.3);
% 
%     % error bars (SEM)
%     errorbar(dIdx-offset, meanY(dIdx), semY(dIdx), 'k', 'LineStyle','none', 'LineWidth',1);
%     errorbar(dIdx+offset, meanO(dIdx), semO(dIdx), 'k', 'LineStyle','none', 'LineWidth',1);
% 
%     % ---- scatter every rat ----
%     jitter = 0.04;
%     scatter( dIdx-offset + jitter*randn(nYoung,1), H_young(:,dIdx), 30, ...
%              'MarkerEdgeColor','none','MarkerFaceColor',youngColor);
%     scatter( dIdx+offset + jitter*randn(nOld,1),   H_old(:,dIdx),   30, ...
%              'MarkerEdgeColor','none','MarkerFaceColor',oldColor);
% end
% 
% xticks(x); xticklabels(string(uniqueDays));
% xlabel('Day'); ylabel('Shannon entropy (bits)');
% legend({'Mean Young','Mean Old'},'Location','best');
% title('Transition Matrix Entropy Calculation');
% pubify_figure_axis_robust(14,14);
% exportgraphics(f3, fullfile(fig_trans_dir, 'Bar_RatTransitionEntropy.png'), ...
%     'Resolution', 450);  % Set to 300 DPI

%% Calculate Entropy on Aggregate Matrix
H_young_all = nan(1,nDays);
H_old_all   = nan(1,nDays);

for dIdx = 1:nDays
    % pull out the aggregated discrete‐transition matrix for that day
    P_y = youngAggDiscrete{dIdx};
    P_o = oldAggDiscrete{dIdx};

    % now simply call the fixed matrixEntropy
    H_young_all(dIdx) = matrixEntropy(P_y);
    H_old_all(dIdx)   = matrixEntropy(P_o);
end

% plot
f4=figure; hold on
errorbar(1:nDays, H_young_all, zeros(1,nDays), 'o-','LineWidth',2,'Color',youngColor);
errorbar(1:nDays, H_old_all,   zeros(1,nDays),   'o-','LineWidth',2,'Color',oldColor);
xticks(1:nDays); xticklabels(string(uniqueDays));
xlabel('Day'); ylabel('Entropy (bits)');
legend('Young','Old','Location','northeast');
xlim([0.9 4.1])
%ylim([1 2.4])
title('Entropy of The Aggregated Discrete Transition Matrix');
pubify_figure_axis_robust(15,15)
exportgraphics(f4, fullfile(fig_trans_dir, 'Agg_RatTransitionEntropy.png'), ...
     'Resolution', 450); 

%% ——— Compute Spectral Gap & Stationary‐Distribution Entropy per Rat per Day ——— %%
epsEnt = 1e-12;    % to guard against log2(0)

nYoung = numel(youngRatIDs);
nOld   = numel(oldRatIDs);
nDays  = totalDays;

gap_young = NaN(nYoung, nDays);
gap_old   = NaN(nOld,   nDays);
H_young   = NaN(nYoung, nDays);
H_old     = NaN(nOld,   nDays);

for i = 1:nYoung
    for d = 1:nDays
        P = youngDailyDiscreteTrans{i,d};
        if all(~isnan(P(:)))
            % Spectral gap
            ev   = eig(P);
            lam  = sort(abs(ev), 'descend');
            gap_young(i,d) = 1 - lam(2);

            %Stationary distribution entropy
            [V,D] = eig(P.');
            % find eigenvalue closest to 1
            [~, idx] = min(abs(diag(D) - 1));
            pi = real(V(:,idx));
            pi = pi / sum(pi);              % normalize
            % clip any small negatives and renormalize
            pi(pi<0) = 0;
            pi = pi / sum(pi);
            H_young(i,d) = -sum(pi .* log2(pi + epsEnt));
        end
    end
end

for i = 1:nOld
    for d = 1:nDays
        P = oldDailyDiscreteTrans{i,d};
        if all(~isnan(P(:)))
            ev   = eig(P);
            lam  = sort(abs(ev), 'descend');
            gap_old(i,d) = 1 - lam(2);

            [V,D] = eig(P.');
            [~, idx] = min(abs(diag(D) - 1));
            pi = real(V(:,idx));
            pi = pi / sum(pi);
            pi(pi<0) = 0;
            pi = pi / sum(pi);
            H_old(i,d) = -sum(pi .* log2(pi + epsEnt));
        end
    end
end

%——— Assemble Tables for ANOVA ——— %%
allGap = [gap_young; gap_old];
allH   = [H_young;   H_old];
Age    = [ repmat({'young'},nYoung,1); repmat({'old'},nOld,1) ];
dayVars = strcat('Day', string(1:nDays));

T_gap = array2table(allGap, 'VariableNames', dayVars);
T_gap.Age = Age;

T_ent = array2table(allH, 'VariableNames', dayVars);
T_ent.Age = Age;

%——— Run Mixed‐Design ANOVAs ——— %%
anovaGap = runMixedANOVA(   T_gap, dayVars, 'Day');

anovaEnt = runMixedANOVA(   T_ent, dayVars, 'Day');
disp(anovaEnt)

% Optional Tukey post‐hoc
postHocGap = runTukeyPostHocMixed(T_gap, dayVars);
postHocEnt = runTukeyPostHocMixed(T_ent, dayVars);

% ——— Between‐Age T‐Tests on Stationary‐Distribution Entropy ——— %%
pEntDay    = nan(1, nDays);
for dIdx = 1:nDays
    y = H_young(:,dIdx);
    o = H_old(:,  dIdx);
    % unequal‐variance two‐sample t‐test
    [~, pEntDay(dIdx)] = ttest2(y, o, 'Vartype','unequal');
end

% Bonferroni correction across nDays tests
pEntDay_bonf = min(pEntDay * nDays, 1);

fprintf('\nStationary‐Distribution Entropy: Young vs Old per day\n');
for dIdx = 1:nDays
    day = uniqueDays(dIdx);
    fprintf(' Day %d: p = %.4f   (Bonferroni p = %.4f)\n', ...
            day, pEntDay(dIdx), pEntDay_bonf(dIdx));
end


figure; hold on
errorbar(1:nDays, nanmean(H_young), nanstd(H_young)./sqrt(nYoung), 'o-','LineWidth',2,'Color',youngColor);
errorbar(1:nDays, nanmean(H_old),   nanstd(H_old)./sqrt(nOld),   'o-','LineWidth',2,'Color',oldColor);
xticks(1:nDays);
xticklabels(string(uniqueDays));
xlabel('Day');
ylabel('Stationary Distribution');
legend('Young','Old','Location','best');
grid on;
title('Stationary Distribtuion of Discrete Transition Matrices');

%% Entropy on Probability Transition  Matrix
% Get your actual days
uniqueDays = sort(unique(data1.Day));
nDays      = numel(uniqueDays);

% Compute means ± SEM
meanGapY = nanmean(gap_young,1);
semGapY  = nanstd(gap_young,[],1)./sqrt(sum(~isnan(gap_young),1));
meanGapO = nanmean(gap_old,  [],1);
semGapO  = nanstd(gap_old,  [],1)./sqrt(sum(~isnan(gap_old),1));

meanHY   = nanmean(H_young, 1);
semHY    = nanstd(H_young,  [],1)./sqrt(sum(~isnan(H_young),1));
meanHO   = nanmean(H_old,    1);
semHO    = nanstd(H_old,    [],1)./sqrt(sum(~isnan(H_old),1));

% ---- Figure 1: Spectral Gap ----
figure; hold on
errorbar(1:nDays, meanGapY, semGapY, 'o-','LineWidth',1.5);
errorbar(1:nDays, meanGapO, semGapO, 's--','LineWidth',1.5);
xticks(1:nDays); xticklabels(string(uniqueDays));
xlabel('Day'); ylabel('Spectral gap (1 - |\lambda_2|)');
legend('Young','Old','Location','best');
title('Spectral gap by Day and Age');
grid on

% ---- Figure 2: Stationary‐Distribution Entropy ----
figure; hold on
errorbar(1:nDays, meanHY, semHY, 'o-','LineWidth',1.5);
errorbar(1:nDays, meanHO, semHO, 's--','LineWidth',1.5);
xticks(1:nDays); xticklabels(string(uniqueDays));
xlabel('Day'); ylabel('Entropy of \pi (bits)');
legend('Young','Old','Location','best');
title('Stationary‐Distribution Entropy by Day and Age');
grid on

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
 oldColor=[0.4157,0.1059,0.6039]; %purple
 youngColor=[0.2196,0.5569,0.2353];% green
% Generate colormaps with increasing intensity (light to dark)
nShades = 256;
youngCMap = [linspace(1, youngColor(1), nShades)', ...
             linspace(1, youngColor(2), nShades)', ...
             linspace(1, youngColor(3), nShades)'];

oldCMap = [linspace(1, oldColor(1), nShades)', ...
           linspace(1, oldColor(2), nShades)', ...
           linspace(1, oldColor(3), nShades)'];
f1=figure('Position',[-2505,-898,2063,1103]);
for d = 1:totalDays
    % --- Young Group Heatmap ---
    subplot(2, totalDays, d);
    h1 = heatmap(strategy_titles, strategy_titles, youngAggDiscrete{d});
    h1.Title = ['Day ' num2str(d)];
    h1.XLabel = '\bfNext Strategy';
    h1.FontSize = 14;
    h1.ColorLimits = [0 0.5];
    h1.GridVisible = 'on';
    h1.CellLabelColor = 'black';
    h1.Colormap = youngCMap;

    if d == 1
        h1.YLabel = '\bfCurrent Strategy';
    else
        h1.YLabel = '';
    end

    % Only show colorbar for last day
    h1.ColorbarVisible = (d == totalDays);

    % --- Old Group Heatmap ---
    subplot(2, totalDays, d + totalDays);
    h2 = heatmap(strategy_titles, strategy_titles, oldAggDiscrete{d});
    h2.Title = [' Day ' num2str(d)];
    h2.XLabel = '\bfNext Strategy';
    h2.FontSize = 14;
    h2.ColorLimits = [0 0.5];
    h2.GridVisible = 'on';
    
    h2.CellLabelColor = 'black';
    h2.Colormap = oldCMap;

    if d == 1
        h2.YLabel = '\bfCurrent Strategy';
    else
        h2.YLabel = '';
    end

    h2.ColorbarVisible = (d == totalDays);
end

sgtitle('Strategy Transition Probabilities', 'FontSize', 18, 'FontWeight', 'bold');
saveas(f1, fullfile(fig_trans_dir, 'DiscreteTransitionProbabilities'), 'png');
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

%% HEAT MAP Prob Transition Matrix
f1=figure('Position',[-2505,-898,2063,1103]);
for d = 1:totalDays
    % --- Young Group Heatmap ---
    subplot(2, totalDays, d);
    h1 = heatmap(strategy_titles, strategy_titles, youngAggProb{d});
    h1.Title = ['Day ' num2str(d)];
    h1.XLabel = '\bfNext Strategy';
    h1.FontSize = 14;
    h1.ColorLimits = [0 0.5];
    h1.GridVisible = 'on';
    h1.CellLabelColor = 'black';
    h1.Colormap = youngCMap;

    if d == 1
        h1.YLabel = '\bfCurrent Strategy';
    else
        h1.YLabel = '';
    end

    % Only show colorbar for last day
    h1.ColorbarVisible = (d == totalDays);

    % --- Old Group Heatmap ---
    subplot(2, totalDays, d + totalDays);
    h2 = heatmap(strategy_titles, strategy_titles, oldAggProb{d});
    h2.Title = [' Day ' num2str(d)];
    h2.XLabel = '\bfNext Strategy';
    h2.FontSize = 14;
    h2.ColorLimits = [0 0.5];
    h2.GridVisible = 'on';
    
    h2.CellLabelColor = 'black';
    h2.Colormap = oldCMap;

    if d == 1
        h2.YLabel = '\bfCurrent Strategy';
    else
        h2.YLabel = '';
    end

    h2.ColorbarVisible = (d == totalDays);
end

sgtitle('Strategy Transition Probabilities', 'FontSize', 18, 'FontWeight', 'bold');
saveas(f1, fullfile(fig_trans_dir, 'ProbTransitionProbabilities'), 'png');
close;
 % Probabilistic transitions heat map
% f3=figure('Position', [100, 100, 2200, 1200]);
% for d = 1:totalDays
% 
%     subplot(2, totalDays, d);
%     h1 = heatmap(strategy_titles, strategy_titles, youngAggProb{d});
%     h1.Title = ['Young: Day' num2str(d)];
%     h1.XLabel = 'Next Trial';
%     h1.YLabel = 'Current Trial';
%     clim([0 0.4])
%     subplot(2, totalDays, d+totalDays);
%     h2 = heatmap(strategy_titles, strategy_titles, oldAggProb{d});
%     h2.Title = ['Old: Day' num2str(d)];
%     h2.XLabel = 'Next Trial';
%     h2.YLabel = 'Current Trial';
%     clim([0 0.4])
% end
% sgtitle('Strategy Confidence Intervals in a Transition Matrix ')
% saveas(f3,fullfile(fig_trans_dir,'Prob_TrnsitionMatrix'),'png');
% close;

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


