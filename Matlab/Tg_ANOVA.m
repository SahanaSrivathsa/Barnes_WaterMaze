%% RUN FIRST
% Rtrack strategy sheet
strategy_data = readtable('/Users/miasponseller/Desktop/Lab/Rtrack/Tg/Tg_MWM_results_06-17-2025.xlsx', VariableNamingRule='preserve');

% Age column to numeric
strategy_data.Age = str2double(strategy_data.Age);

% Age groups
age = strategy_data.Age;
age_group = NaN(height(strategy_data), 1);
age_group(age >= 4 & age <= 6) = 5;
age_group(age >= 7 & age <= 9) = 8;
age_group(age >= 10 & age <= 12) = 11;
age_group(age >= 13 & age <= 14) = 13.5;
age_group(age >= 15 & age <= 16) = 15.5;
age_group(age >= 17 & age <= 19) = 18;
age_group(age >= 20 & age <= 21) = 20; % only have age 20, so calling it 20 for now
strategy_data.age_group = categorical(age_group);

% Convert to categorical
strategy_data.Sex = categorical(strategy_data.Sex);
strategy_data.APP = categorical(strategy_data.APP);

% Strategy categories
strat_cats = {
    {'thigmotaxis', 'circling', 'random path'}, 'NonGoal';
    {'scanning', 'chaining'}, 'Procedural';
    {'directed search', 'corrected path', 'direct path', 'perseverance'} 'Allocentric'
    };

cat_labels = strat_cats(:, 2);
cat_names = {'Non-Goal Oriented', 'Procedural', 'Allocentric'};
n_cats = numel(cat_labels);

rats = unique(strategy_data.("_TargetID"));
n_rats = numel(rats);
unique_days = unique(strategy_data.("_Day"));
n_days = numel(unique_days);

%% Strategy Category Average Table

rat_means = nan(n_rats, n_cats * n_days);
AgeGroup = categorical(zeros(n_rats, 1));
Sex = categorical(zeros(n_rats, 1));
Genotype = categorical(zeros(n_rats, 1));

for i = 1:n_rats
    ridx = strcmp(strategy_data.("_TargetID"), rats{i});
    
    % Assign Age, Sex, and APP for this rat
    AgeGroup(i) = strategy_data.age_group(find(ridx, 1));
    Sex(i) = strategy_data.Sex(find(ridx, 1));
    Genotype(i) = strategy_data.APP(find(ridx, 1));

    col = 1;
    for c = 1:n_cats
        strat_prob = zeros(sum(ridx), 1);
        
        % Extract strategy names for category c
        strategy_list = strat_cats{c, 1};
        
        for s = 1:numel(strategy_list)
            strat_prob = strat_prob + strategy_data.(strategy_list{s})(ridx);
        end
        
        for d = 1:n_days
            pdx = ridx & (strategy_data.("_Day") == unique_days(d));
            
            % Mean of summed strategy probabilities for this day
            rat_means(i, col) = mean(strat_prob(strategy_data.("_Day")(ridx) == unique_days(d)));
            col = col + 1;
        end
    end
end

% Build table
var_names = {'RatID', 'AgeGroup', 'Sex', 'Genotype'};

for g = 1:n_cats
    for d = 1:n_days
        var_names{end + 1} = sprintf('%s_Day%d', cat_labels{g}, unique_days(d));
    end
end

% Convert rat_means into a cell array (one variable per column)
rat_means_cell = num2cell(rat_means, 1);  

% Build final table
rat_table = table(rats, AgeGroup, Sex, Genotype, ...
                 rat_means_cell{:}, ...
                 'VariableNames', var_names);



%% Run ANOVA and post-hoc

withinDesign = table((1:n_days)', 'VariableNames', {'Day'});
withinDesign.Day = categorical(withinDesign.Day);

% Get unique combinations of between-subject factors
[grpIDs, sexGrp, genoGrp, ageGrp] = findgroups(rat_table.Sex, rat_table.Genotype, rat_table.AgeGroup);

% Count rats per group
counts = splitapply(@numel, grpIDs, grpIDs);

% Store unique group combinations as a table
groupTable = table(sexGrp, genoGrp, ageGrp, counts, 'VariableNames', {'Sex', 'Genotype', 'AgeGroup', 'Count'});

% Filter groups with at least 3 rats
validGroups = groupTable(groupTable.Count >= 3, :);

for c = 1:numel(cat_labels)
    measureVars = cell(1, n_days);
    for d = 1:n_days
        measureVars{d} = sprintf('%s_Day%d', cat_labels{c}, unique_days(d));
    end

    for g = 1:height(validGroups)
        % Extract current group's factor values
        curSex = validGroups.Sex(g);
        curGenotype = validGroups.Genotype(g);
        curAgeGroup = validGroups.AgeGroup(g);

        % Filter rat_table to only rats in this group
        rowsInGroup = rat_table.Sex == curSex & ...
                      rat_table.Genotype == curGenotype & ...
                      rat_table.AgeGroup == curAgeGroup;

        subTable = rat_table(rowsInGroup, :);

        % Skip if less than 3 rats in the group
        if height(subTable) < 3
            continue
        end

        % Build formula for within-subjects factor (Day)
        formula = sprintf('%s-%s ~ 1', measureVars{1}, measureVars{end});

        try
            rm = fitrm(subTable, formula, 'WithinDesign', withinDesign);
            anovaResults = ranova(rm, 'WithinModel', 'Day');

            fprintf('ANOVA results for category: %s, Group: Sex=%s, Genotype=%s, AgeGroup=%s\n', ...
                    cat_labels{c}, string(curSex), string(curGenotype), string(curAgeGroup));
            disp(anovaResults);

            % Post-hoc pairwise comparisons for Day
            posthocResults = multcompare(rm, 'Day');
            fprintf('Post-hoc results:\n');
            disp(posthocResults);

        catch ME
            warning('ANOVA failed for group %d: %s', g, ME.message);
        end
    end
end
