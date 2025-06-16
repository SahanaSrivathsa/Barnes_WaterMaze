library(tidyr)
library(dplyr)
library(ggplot2)
library(readxl)
library(scales)
library(Rtrack)

# Load data
strat_sheet <- read_excel('/Users/miasponseller/Desktop/Lab/Rtrack/Tg/Tg_MWM_results_06-16-2025.xlsx')
all_rats_spatial <- read.csv('/Users/miasponseller/Desktop/Lab/Rtrack/Tg/Tg_AllRats_Spatial_cleaned.csv')
description_file <- '/Users/miasponseller/Desktop/Lab/Rtrack/Tg/Tg_exp_desc.xlsx'

# Add if running Rtrack plots
all_trials <- '/Users/miasponseller/Desktop/Lab/Rtrack/Tg/All Tg Tracks'

# Figures Folder
fig_folder = '/Users/miasponseller/Desktop/Lab/Rtrack/Tg/Figures'
dir.create(fig_folder, recursive = TRUE, showWarnings = FALSE)

sex_counts <- strat_sheet %>% 
  distinct(`_TargetID`, .keep_all = TRUE) %>% 
  count(Sex)
#print(sex_counts)

n_male = sex_counts$n[sex_counts$Sex == 'M']
n_female = sex_counts$n[sex_counts$Sex == 'F']

# Strategy categories
non_goal_oriented <- c('thigmotaxis', 'circling', 'random path')
procedural <- c('scanning', 'chaining')
allocentric <- c('directed search', 'corrected path', 'direct path')

list_of_strats = c('direct path', 'corrected path', 'directed search', 'chaining', 
                   'scanning', 'random path', 'circling', 'thigmotaxis')

list_of_strat_cats = c('Allocentric', 'Procedural', 'NonGoalOriented')

# Add StratCat and AgeCat columns
strat_sheet_1 <- strat_sheet %>% 
  mutate(AgeCat = case_when(
    Age <= 13 ~ 'Young',
    Age >13 & Age <= 17 ~ 'Middle',
    Age > 17 ~ 'Old'
  )) %>% 
  mutate(
    StratCat = case_when(
      name %in% non_goal_oriented ~ 'NonGoalOriented',
      name %in% procedural ~ 'Procedural',
      name %in% allocentric ~ 'Allocentric'
    )
  )

# Add Group column
strat_sheet_1 <- strat_sheet_1 %>% 
  mutate(Group = paste0(Sex, '-', AgeCat, '-', APP))

# StratCat levels
strat_sheet_1$StratCat <- factor(
  strat_sheet_1$StratCat, levels = c('NonGoalOriented', 'Procedural', 'Allocentric')
)

# Add Allocentric, Procedural, and Non-Goal Oriented columns with avg prob for trial
strat_sheet_1$NonGoalOrientedProb <- rowSums(strat_sheet_1[,non_goal_oriented])
strat_sheet_1$ProceduralProb <- rowSums(strat_sheet_1[,procedural])
strat_sheet_1$AllocentricProb <- rowSums(strat_sheet_1[,allocentric])

# colors
strat_colors <- c(
  'direct path' = '#2B6519',
  'corrected path' = '#5F972B',
  'directed search' = '#A3CA3F',
  'chaining' = '#FCF050',
  'scanning' = '#F7D148',
  'random path' = '#F3B341',
  'circling' = '#AE7A38',
  'thigmotaxis' = '#69403F'
)

strat_cat_colors <- c(
  'Allocentric' = '#1b7837',
  'Procedural' = '#ffd700',
  'NonGoalOriented' = '#654321'
)


# Counts ------------------------------------------------------------------

sex_counts <- strat_sheet_1 %>% 
  filter(Sex %in% c('M', 'F')) %>% 
  group_by(Sex) %>% 
  summarize(UniqueCount = n_distinct(`_TargetID`))
print(sex_counts)

group_counts <- strat_sheet_1 %>% 
  group_by(Group) %>% 
  summarize(n_animals = n_distinct(`_TargetID`), .groups = 'drop')
print(group_counts)

# Rtrack Strategy Plots ---------------------------------------------------

bulk_strategy_calling <- function() {
  experiment_file <- description_file
  
  experiment <<- Rtrack::read_experiment(experiment_file, data.dir = all_trials)
  strategies <<- Rtrack::call_strategy(experiment)
  list(experiment = experiment, strategies = strategies)
}

strategy_plots <- function() {
  # if bulk_strategy_calling() has not alrady been run, run it
  if (!exists("experiment", envir = .GlobalEnv) || !exists("strategies", envir = .GlobalEnv)) {
    bulk_strategy_calling()
  }
  
  par(cex.lab = 1.4,   # Adjust axis label size
      cex.axis = 1.3,  # Adjust tick mark labels
      cex.main = 1.6,  # Main title size
      cex.sub = 1.2)   # Subtitle size
  
  # Strategy plot, across all rats
  Rtrack::plot_strategies(strategies, experiment = experiment)
  
  # Strategy plots, by Cohort
  #Rtrack::plot_strategies(strategies, experiment = experiment, factor = "Cohort")
  
  # Strategy plots, by Age
  Rtrack::plot_strategies(strategies, experiment = experiment, factor = "Sex")
  
  print("Plots successfully created")
}

#strategy_plots()


# Strategy Use Proportions by Rat -----------------------------------------

# List of unique groups
groups <- unique(strat_sheet_1$Group)

# Loop through groups and generate plots for each
for (g in groups) {
  group_data <- strat_sheet_1 %>% 
    filter(Group == g) %>% 
    group_by(`_Day`, name) %>% 
    summarize(count = n(), .groups = 'drop') %>% 
    group_by(`_Day`) %>% 
    mutate(prop = count/sum(count)) %>% 
    ungroup()
  
  plt <- ggplot(group_data, aes(x = factor(`_Day`), y = prop, 
                                fill = factor(name, levels = list_of_strats))) +
    geom_bar(stat = 'identity', position = 'stack') +
    scale_fill_manual(values = strat_colors) + 
    scale_y_continuous() +
    labs(
      title = paste('Dominant Strategy Usage by Day:', g),
      x = 'Day', y = 'Proportion of Trials',
      fill = 'Strategy'
    ) + 
    theme_minimal(base_size = 14) +
    theme(
      axis.title = element_text(size = 16, face = 'bold'),
      axis.text = element_text(size = 14, face = 'bold'),
      plot.title = element_text(size = 18, face = 'bold', hjust = .5),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12)
    )
  
  # Save plot to Figures folder
  # safe_g <- gsub("[^A-Za-z0-9_]", "_", g)
  # file_path <- file.path(fig_folder, paste0("Dominant_Strategy_Rats_Days_", safe_g, ".jpeg"))
  # ggsave(file_path, plt, width = 8, height = 6, dpi = 300)
  
  print(plt)
}








