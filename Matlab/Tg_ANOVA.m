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

% convert between-subject variables to categorical
strategy_data.Sex = categorical(strategy_data.Sex);
strategy_data.APP = categorical(strategy_data.APP);

% Strategy categories
strat_cats = struct();
strat_cats.Allocentric = {'direct path', 'corrected path', 'directed search'};
strat_cats.Procedural = {'scanning', 'chaining'};
strat_cats.NonGoalOriented = {'thigmotaxis', 'circling', 'random path'};

strat_cat_names = fieldnames(strat_cats);
%% Strategy Category Average Table

% Make column names valid
all_strategies = [strat_cats.Allocentric, strat_cats.Procedural, strat_cats.NonGoalOriented];
valid_strategy_vars = matlab.lang.makeValidName(all_strategies); 
strategy_data.Properties.VariableNames = matlab.lang.makeValidName(strategy_data.Properties.VariableNames);

% Convert x_TargetID and x_Day to categorical
strategy_data.x_TargetID = categorical(strategy_data.x_TargetID);
strategy_data.x_Day = categorical(strategy_data.x_Day);

% Unique rats and days
rats = unique(strategy_data.x_TargetID);
days = unique(strategy_data.x_Day);

% Output table
category_avg_table = table();

for i = 1:length(rats)
    for j = 1:length(days)
        rat_id = rats(i);
        day = days(j);

        % Get trials for this rat and day
        mask = strategy_data.x_TargetID == rat_id & strategy_data.x_Day == day;
        subset = strategy_data(mask, :);

        % Row with ID and between-subject variables
        row = table();
        row.Animal = rat_id;
        row.Day = day;
        row.Sex = subset.Sex(1);
        row.APP = subset.APP(1);
        row.AgeGroup = subset.age_group(1);

        % Calculate averages per strategy category
        for c = 1:length(strat_cat_names)
            cat = strat_cat_names{c};
            strat_names = matlab.lang.makeValidName(strat_cats.(cat));
            values = table2array(subset(:, strat_names));
            row.(cat) = mean(values(:), 'omitnan');
        end

        category_avg_table = [category_avg_table; row];
    end
end

%% Run ANOVA and post-hoc for each category
