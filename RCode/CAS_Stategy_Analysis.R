library(tidyr)
library(dplyr)
library(ggplot2)
library(readxl)
library(scales)

# Load data
strat_sheet <- read_excel('/Users/miasponseller/Desktop/Lab/Rtrack/CAS/CAS_MWM_results_05-28-2025.xlsx')
all_rats_spatial <- read.csv('/Users/miasponseller/Desktop/Lab/Rtrack/CAS/CAS_AllRats_Spatial.csv')

# Figures Folder
fig_folder = '/Users/miasponseller/Desktop/Lab/Rtrack/CAS/Figures'
dir.create(fig_folder, recursive = TRUE, showWarnings = FALSE)

age_counts <- all_rats_spatial %>% 
  distinct(Animal, .keep_all = TRUE) %>% 
  count(Age)
print(age_counts)

n_young = age_counts$n[age_counts$Age == '6mo']
n_middle = age_counts$n[age_counts$Age == '15mo']
n_old = age_counts$n[age_counts$Age == '23mo']

# Age groups
age_groups <- c('young', 'middle', 'old')

# Strategy categories
non_goal_oriented <- c('thigmotaxis', 'circling', 'random path')
procedural <- c('scanning', 'chaining')
allocentric <- c('directed search', 'corrected path', 'direct path')

list_of_strats = c('direct path', 'corrected path', 'directed search', 'chaining', 
                   'scanning', 'random path', 'circling', 'thigmotaxis')
list_of_strat_cats = c('Allocentric', 'Procedural', 'NonGoalOriented')

# Add Age/Performance column
strat_sheet_1 <- strat_sheet %>% 
  mutate(Group = paste0(Age, '/', Performance))

# Add Age/Performance and StratCat columns
strat_sheet_1 <- strat_sheet %>% 
  mutate(Group = paste0(Age, "/", Performance)) %>% 
  mutate(
    StratCat = case_when(
      name %in% non_goal_oriented ~ 'NonGoalOriented',
      name %in% procedural ~ 'Procedural',
      name %in% allocentric ~ 'Allocentric'
    )
  )

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

# Dominant Strat Across Days for Age/Perf Groups --------------------------

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
      x = 'Day', y = 'Proportion of Rats',
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
  safe_g <- gsub("[^A-Za-z0-9_]", "_", g)
  file_path <- file.path(fig_folder, paste0("Dominant_Strategy_Days_", safe_g, ".jpeg"))
  ggsave(file_path, plt, width = 8, height = 6, dpi = 300)
  
  print(plt)
}


# Dominant Strat Type Across Days for Age/Perf Groups ---------------------

# List of unique groups
groups <- unique(strat_sheet_1$Group)

# Loop through groups and generate plots for each
for (g in groups) {
  group_data <- strat_sheet_1 %>% 
    filter(Group == g) %>% 
    group_by(`_Day`, StratCat) %>% 
    summarize(count = n(), .groups = 'drop') %>% 
    group_by(`_Day`) %>% 
    mutate(prop = count/sum(count)) %>% 
    ungroup()
  
  plt <- ggplot(group_data, aes(x = factor(`_Day`), y = prop, 
                                fill = factor(StratCat, levels = list_of_strat_cats))) +
    geom_bar(stat = 'identity', position = 'stack') +
    scale_fill_manual(values = strat_cat_colors) + 
    scale_y_continuous() +
    labs(
      title = paste('Dominant Strategy Category Usage by Day:', g),
      x = 'Day', y = 'Proportion of Trials',
      fill = 'Category'
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
  safe_g <- gsub("[^A-Za-z0-9_]", "_", g)
  file_path <- file.path(fig_folder, paste0("Dominant_Strat_Cat_Days_", safe_g, ".jpeg"))
  ggsave(file_path, plt, width = 8, height = 6, dpi = 300)
  
  print(plt)
}


# Dominant Strat by Performance Groups Each Day ---------------------------

# List of unique age groups
perf_groups <- unique(strat_sheet_1$Performance)

for (p in perf_groups) {
  
  perf_data <- strat_sheet_1 %>% 
    filter(Performance == p)
  
  # Mean and SE
  summary_df <- perf_data %>% 
    group_by(`_Day`, StratCat) %>% 
    summarize(
      mean_prob = mean(confidence, na.rm = TRUE),
      se = sd(confidence, na.rm = TRUE)/sqrt(n()),
      .groups = 'drop'
    )
  
  plt <- ggplot(summary_df, aes(x = factor(`_Day`), y = mean_prob, fill = StratCat)) +
    geom_bar(stat = 'identity', position = position_dodge(width = .9), alpha = .5) +
    geom_jitter(data = perf_data,
                aes(x = factor(`_Day`), y = confidence, color = StratCat),
                position = position_jitterdodge(jitter.width = .2, dodge.width = .9),
                size = 2, alpha = .7,
                inherit.aes = FALSE,
                show.legend = FALSE) +
    geom_errorbar(aes(ymin = mean_prob - se, ymax = mean_prob + se), 
                  position = position_dodge(width = .9),
                  width = .2) +
    labs(
      title = paste('Strategy Category Usage by Day -', tools::toTitleCase(p), 'Performance'),
      x = 'Day', y = 'Probability of Strategy Category',
      fill = 'Category'
    ) +
    scale_fill_manual(values = strat_cat_colors) +
    scale_color_manual(values = strat_cat_colors) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_minimal(base_size = 14) +
    theme(
      axis.title = element_text(size = 16, face = 'bold'),
      axis.text = element_text(size = 14, face = 'bold'),
      plot.title = element_text(size = 18, face = 'bold', hjust = .5),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12))

  print(plt)
  
  # Save plot to Figures folder
  safe_g <- gsub("[^A-Za-z0-9_]", "_", g)
  file_path <- file.path(fig_folder, paste0("Dominant_Strat_Cat_Days_Performance", safe_g, ".jpeg"))
  ggsave(file_path, plt, width = 8, height = 6, dpi = 300)
  
}


# Prob Strat Cat Age/Perf -------------------------------------------------

# Average strategy category by rat, day, and group
rat_day_avg <- strat_sheet_1 %>% 
  group_by(Group, `_TargetID`, `_Day`) %>% 
  summarize(
    Allocentric = mean(AllocentricProb, na.rm = TRUE),
    Procedural = mean(ProceduralProb, na.rm = TRUE),
    NonGoalOriented = mean(NonGoalOrientedProb, na.rm = TRUE),
    .groups = 'drop'
  )

# Pivot long
rat_day_long <- rat_day_avg %>% 
  pivot_longer(cols = c(Allocentric, Procedural, NonGoalOriented),
               names_to = 'StrategyCategory',
               values_to = 'Probability')

group_day_summary <- rat_day_long %>% 
  group_by(Group, `_Day`, StrategyCategory) %>% 
  summarize(MeanProb = mean(Probability, na.rm = TRUE),
            SD = sd(Probability, na.rm = TRUE),
            SEM = SD/sqrt(n()),
            .groups = 'drop')

for (g in unique(group_day_summary$Group)) {
  
  # Filter group
  summary_data <- group_day_summary %>% filter(Group == g)
  rat_data <- rat_day_long %>% filter(Group == g)
  
  plt <- ggplot(summary_data, aes(x = as.factor(`_Day`), y = MeanProb, fill = StrategyCategory)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = .8, alpha = .5) +
    geom_jitter(data = rat_data,
                aes(x = as.factor(`_Day`), y = Probability, color = StrategyCategory),
                position = position_jitterdodge(jitter.width = .2, dodge.width = .9),
                alpha = .8, size = 2, inherit.aes = FALSE) +
    geom_errorbar(aes(ymin = MeanProb - SEM, ymax = MeanProb + SEM),
                  width = .2, position = position_dodge(width = .9)) +
    labs(
      title = paste("Strategy Use:", g),
      x = "Day",
      y = "Probability of Strategy Use",
      fill = "Strategy Category",
      color = 'Strategy Category'
    ) +
    scale_fill_manual(values = strat_cat_colors) +
    scale_color_manual(values = strat_cat_colors) +
    coord_cartesian(ylim = c(0,1)) +
    theme_minimal(base_size = 14)
  
  # Print the plot
  print(plt)
  
  # Save plot to Figures folder
  # safe_g <- gsub("[^A-Za-z0-9_]", "_", g)
  # file_path <- file.path(fig_folder, paste0("AvgStratCatUse_DayxPerfGroup", safe_g, ".jpeg"))
  # ggsave(file_path, p, width = 8, height = 6, dpi = 300)
  
}


  

# Group Counts ------------------------------------------------------------

group_counts <- strat_sheet_1 %>% 
  distinct(`_TargetID`, Group, .keep_all = TRUE) %>% 
  count(Group, sort = TRUE)

#print(group_counts)















