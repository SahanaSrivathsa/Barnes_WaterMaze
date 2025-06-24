library(tidyr)
library(dplyr)
library(ggplot2)
library(readxl)
library(scales)
library(Rtrack)
library(afex)
library(emmeans)
library(colorspace)
library(ggpubr)


# MUST RUN ----------------------------------------------------------------
# Load data
strat_sheet <- read_excel('/Users/miasponseller/Desktop/Lab/Rtrack/Tg/Tg_MWM_results_06-17-2025.xlsx')
all_rats_spatial <- read.csv('/Users/miasponseller/Desktop/Lab/Rtrack/Tg/Tg_AllRats_Spatial_cleaned.csv')
description_file <- '/Users/miasponseller/Desktop/Lab/Rtrack/Tg/Tg_exp_desc.xlsx'

# Add if running Rtrack plots
all_trials <- '/Users/miasponseller/Desktop/Lab/Rtrack/Tg/All Tg Tracks'

# Fix floating-point precision issue
# strat_sheet$Age <- as.numeric(strat_sheet$Age)
# strat_sheet$Age <- round(strat_sheet$Age, 1)

# Figures Folder
fig_folder = '/Users/miasponseller/Desktop/Lab/Rtrack/Tg/Figures'
dir.create(fig_folder, recursive = TRUE, showWarnings = FALSE)

sex_counts <- strat_sheet %>% 
  group_by(Sex) %>% 
  summarize(UniqueCount = n_distinct(`_TargetID`))

n_male = sex_counts$n[sex_counts$Sex == 'M']
n_female = sex_counts$n[sex_counts$Sex == 'F']

# Strategy categories
non_goal_oriented <- c('thigmotaxis', 'circling', 'random path')
procedural <- c('scanning', 'chaining')
allocentric <- c('directed search', 'corrected path', 'direct path')

list_of_strats = c('direct path', 'corrected path', 'directed search', 'chaining', 
                   'scanning', 'random path', 'circling', 'thigmotaxis')

list_of_strat_cats = c('Allocentric', 'Procedural', 'NonGoalOriented')

# Add StratCat, AgeCat, and Group columns
strat_sheet <- strat_sheet %>% 
  mutate(
    StratCat = case_when(
      name %in% non_goal_oriented ~ 'NonGoalOriented',
      name %in% procedural ~ 'Procedural',
      name %in% allocentric ~ 'Allocentric'
    )) %>% 
  mutate(
    AgeCat = case_when(
      Age %in% c(4, 5, 6) ~ '5',
      Age %in% c(7, 8, 9) ~ '8', 
      Age %in% c(10, 11, 12) ~ '11',
      Age %in% c(13, 14) ~ '13.5',
      Age %in% c(15, 16) ~ '15.5',
      Age %in% c(20) ~ '20'
    )) %>% 
  mutate(Group = paste0(Sex, '-', APP))

# StratCat levels
strat_sheet$StratCat <- factor(
  strat_sheet$StratCat, levels = c('NonGoalOriented', 'Procedural', 'Allocentric')
)

# Group levels
group_levels <- c(
  "F-5", "M-5",
  "F-8", "M-8",
  "F-11", "M-11",
  "F-13.5", "M-13.5",
  "F-15.5", "M-15.5",
  "F-20", "M-20"
)

age_cats <- c("5", "8", "11", "13.5", "15.5", "20")
sexes <- c("F", "M")
sex_age_groups <- as.vector(outer(sexes, age_cats, paste, sep = "-"))

# Add Allocentric, Procedural, and Non-Goal Oriented columns with avg prob for trial
strat_sheet$NonGoalOrientedProb <- rowSums(strat_sheet[,non_goal_oriented])
strat_sheet$ProceduralProb <- rowSums(strat_sheet[,procedural])
strat_sheet$AllocentricProb <- rowSums(strat_sheet[,allocentric])

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

sex_age_colors <- c(
  "F-5"     = "#e6194b",  # Vivid Red
  "M-5"     = "#f58231",  # Vivid Orange
  "F-8"     = "#3cb44b",  # Vivid Green
  "M-8"     = "#0082c8",  # Vivid Blue
  "F-11"    = "#911eb4",  # Vivid Purple
  "M-11"    = "#46f0f0",  # Bright Cyan
  "F-15.5"  = "#000075",  # Dark Blue (adds contrast)
  "M-15.5"  = "#f032e6",  # Yellow *included optionally* — change if unwanted
  "F-20"    = "#bfef45",  # Lime Green (also bright)
  "M-20"    = "#7f7fff"   # Light Magenta (swap if pastel-like)
)


# Counts ------------------------------------------------------------------

sex_counts <- strat_sheet %>% 
  group_by(Sex) %>% 
  summarize(n_animals = n_distinct(`_TargetID`), .groups = 'drop')
print(sex_counts)

age_counts <- strat_sheet %>% 
  group_by(Age) %>% 
  summarize(n_animals = n_distinct(`_TargetID`))
print(age_counts)

group_counts <- strat_sheet %>% 
  group_by(Group) %>% 
  summarize(n_animals = n_distinct(`_TargetID`), .groups = 'drop')
print(group_counts)

# Group summary, cols Group, n_animals, and one for each AgeCat
group_summary <- strat_sheet %>%
  distinct(Group, AgeCat, `_TargetID`) %>%
  count(Group, AgeCat) %>%
  pivot_wider(names_from = AgeCat, values_from = n, values_fill = 0) %>%
  mutate(n_animals = rowSums(across(-Group))) %>%
  relocate(Group, n_animals, '5', '8', '11', '15.5', '20')
totals_row <- group_summary %>%
  summarise(across(where(is.numeric), sum)) %>%
  mutate(Group = "TOTAL") %>%
  relocate(Group)
final_summary <- bind_rows(group_summary, totals_row)
print(final_summary)

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
groups <- unique(strat_sheet$Group)

# Loop through groups and generate plots for each
for (g in groups) {
  group_data <- strat_sheet %>% 
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
  # file_path <- file.path(fig_folder, paste0("Dominant_Strategy_PropOfTrials", safe_g, ".jpeg"))
  # ggsave(file_path, plt, width = 8, height = 6, dpi = 300)
  
  print(plt)
}

# Prob Strat Cat Sex/Genotype ---------------------------------------------

# Average strategy category by rat, day, and group
rat_day_avg <- strat_sheet %>% 
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
               names_to = 'StratCat',
               values_to = 'Probability') %>% 
  mutate(StratCat = factor(
    StratCat,
    levels = c('NonGoalOriented', 'Procedural', 'Allocentric')
  ))

group_day_summary <- rat_day_long %>% 
  group_by(Group, `_Day`, StratCat) %>% 
  summarize(MeanProb = mean(Probability, na.rm = TRUE),
            SD = sd(Probability, na.rm = TRUE),
            SEM = SD/sqrt(n()),
            .groups = 'drop')

for (g in unique(group_day_summary$Group)) {
  
  # Filter group
  summary_data <- group_day_summary %>% filter(Group == g)
  rat_data <- rat_day_long %>% filter(Group == g)
  
  plt <- ggplot(summary_data, aes(x = as.factor(`_Day`), y = MeanProb, fill = StratCat)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = .8, alpha = .5) +
    geom_jitter(data = rat_data,
                aes(x = as.factor(`_Day`), y = Probability, color = StratCat),
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
  safe_g <- gsub("[^A-Za-z0-9_]", "_", g)
  file_path <- file.path(fig_folder, paste0("AvgStratCatUse_Day_Sex_Genotype_", safe_g, ".jpeg"))
  ggsave(file_path, plt, width = 8, height = 6, dpi = 300)
  
}


# Prob StratCat Sex/Genotype/AgeCat ---------------------------------------
# Average strategy category by rat, day, sex, genotype, and age category
rat_day_avg2 <- strat_sheet %>%
  group_by(Sex, APP, AgeCat, `_TargetID`, `_Day`) %>%
  summarize(
    Allocentric = mean(AllocentricProb, na.rm = TRUE),
    Procedural = mean(ProceduralProb, na.rm = TRUE),
    NonGoalOriented = mean(NonGoalOrientedProb, na.rm = TRUE),
    .groups = 'drop'
  )

# Pivot long
rat_day_long2 <- rat_day_avg2 %>%
  pivot_longer(
    cols = c(Allocentric, Procedural, NonGoalOriented),
    names_to = 'StratCat',
    values_to = 'Probability'
  ) %>%
  mutate(StratCat = factor(
    StratCat,
    levels = c('NonGoalOriented', 'Procedural', 'Allocentric')
  ))

# Summary for plotting
group_day_summary2 <- rat_day_long2 %>%
  group_by(Sex, APP, AgeCat, `_Day`, StratCat) %>%
  summarize(
    MeanProb = mean(Probability, na.rm = TRUE),
    SD = sd(Probability, na.rm = TRUE),
    SEM = SD / sqrt(n()),
    .groups = 'drop'
  )

groups_to_plot <- unique(group_day_summary2 %>% select(Sex, APP, AgeCat))

for(i in seq_len(nrow(groups_to_plot))) {
  current_group <- groups_to_plot[i, ]
  
  summary_data <- group_day_summary2 %>%
    filter(
      Sex == current_group$Sex,
      APP == current_group$APP,
      AgeCat == current_group$AgeCat
    )
  
  rat_data <- rat_day_long2 %>%
    filter(
      Sex == current_group$Sex,
      APP == current_group$APP,
      AgeCat == current_group$AgeCat
    )
  
  plot_title <- paste("Strategy Use:", current_group$Sex, current_group$APP, "Age:", current_group$AgeCat)
  
  plt <- ggplot(summary_data, aes(x = as.factor(`_Day`), y = MeanProb, fill = StratCat)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8, alpha = 0.5) +
    geom_jitter(data = rat_data,
                aes(x = as.factor(`_Day`), y = Probability, color = StratCat),
                position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9),
                alpha = 0.8, size = 2, inherit.aes = FALSE) +
    geom_errorbar(aes(ymin = MeanProb - SEM, ymax = MeanProb + SEM),
                  width = 0.2, position = position_dodge(width = 0.9)) +
    labs(
      title = plot_title,
      x = "Day",
      y = "Probability of Strategy Use",
      fill = "Strategy Category",
      color = "Strategy Category"
    ) +
    scale_fill_manual(values = strat_cat_colors) +
    scale_color_manual(values = strat_cat_colors) +
    coord_cartesian(ylim = c(0, 1)) +
    theme_minimal(base_size = 14)
  
  print(plt)
  
  # Save plot to Figures folder
  # safe_title <- gsub("[^A-Za-z0-9_]", "_", plot_title)
  # file_path <- file.path(fig_folder, paste0("AvgStratCatUse_Day_Sex_Genotype_Age", safe_title, ".jpeg"))
  # ggsave(file_path, plt, width = 8, height = 6, dpi = 300)
}



# Across Days StratCat Line Plots -----------------------------------------

rat_day_avg2 <- strat_sheet %>%
  group_by(Sex, APP, AgeCat, `_TargetID`, `_Day`) %>%
  summarize(
    Allocentric = mean(AllocentricProb, na.rm = TRUE),
    Procedural = mean(ProceduralProb, na.rm = TRUE),
    NonGoalOriented = mean(NonGoalOrientedProb, na.rm = TRUE),
    .groups = 'drop'
  )

# Pivot long
rat_day_long2 <- rat_day_avg2 %>%
  pivot_longer(
    cols = c(Allocentric, Procedural, NonGoalOriented),
    names_to = 'StratCat',
    values_to = 'Probability'
  ) %>%
  mutate(StratCat = factor(
    StratCat,
    levels = c('NonGoalOriented', 'Procedural', 'Allocentric')
  ))

# Summary for plotting
group_day_summary2 <- rat_day_long2 %>%
  group_by(Sex, APP, AgeCat, `_Day`, StratCat) %>%
  summarize(
    MeanProb = mean(Probability, na.rm = TRUE),
    SD = sd(Probability, na.rm = TRUE),
    SEM = SD / sqrt(n()),
    .groups = 'drop'
  )

plot_strategy_lines <- function(data, sex_filter = NULL, title_prefix = "", save_folder = NULL) {
  if (!is.null(sex_filter)) {
    data <- data %>% filter(Sex == sex_filter)
  }
  
  data <- data %>% 
    mutate(
      ColorGroup = paste(Sex, AgeCat, sep = '-'),
      LineType = factor(APP, levels = c('WT', 'APP/+')) # WT solid, APP dashed
    )
  
  for (strategy in levels(data$StratCat)) {
    strategy_data <- data %>% filter(StratCat == strategy)
    
    plt <- ggplot(strategy_data, aes(x = as.factor(`_Day`), y = MeanProb, 
                                     group = interaction(ColorGroup, LineType), 
                                     color = ColorGroup, linetype = LineType)) +
      geom_line(linewidth = 1) +
      geom_point(size = 2) +
      geom_errorbar(aes(ymin = MeanProb - SEM, ymax = MeanProb + SEM), width = 0.1) +
      labs(
        title = paste0(gsub("_", "", title_prefix), " Strategy Use Over Days: ", strategy),
        x = "Day",
        y = "Average Probability of Strategy Use",
        color = "Sex-Age Group",
        linetype = 'Genotype'
      ) +
      scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
      theme_minimal(base_size = 14)
    
    print(plt)
    
    # Save plot to Figures folder
    if (!is.null(save_folder)) {
      safe_g <- paste0(gsub("[^A-Za-z0-9_]", "_", paste0(title_prefix, strategy)), ".jpeg")
      ggsave(file.path(save_folder, safe_g), plt, width = 10, height = 6, dpi = 300)}
    
  }
}

# --- Generate Plots ---

# All groups
plot_strategy_lines(group_day_summary2, sex_filter = NULL, title_prefix = "All_", save_folder = fig_folder)

# Females only
plot_strategy_lines(group_day_summary2, sex_filter = "F", title_prefix = "Female_", save_folder = fig_folder)

# Males only
plot_strategy_lines(group_day_summary2, sex_filter = "M", title_prefix = "Male_", save_folder = fig_folder)



# Across Days Allocentric Line Plots --------------------------------------

# Allocentric strategies to long
allocentric_long <- strat_sheet %>%
  select(Group, `_Day`, Sex, APP, AgeCat, all_of(allocentric)) %>%
  pivot_longer(cols = all_of(allocentric), 
               names_to = "Strategy", 
               values_to = "Probability")

# Summarize
allocentric_summary <- allocentric_long %>%
  group_by(Group, `_Day`, Strategy, Sex, APP, AgeCat) %>%
  summarize(
    MeanProb = mean(Probability, na.rm = TRUE),
    SD = sd(Probability, na.rm = TRUE),
    SEM = SD / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(GroupLabel = paste(Sex, AgeCat, sep = "-"))

plot_allocentric_strategies <- function(data, sex_filter = NULL, title_prefix = "", save_folder = NULL) {
  if (!is.null(sex_filter)) {
    data <- data %>% filter(Sex == sex_filter)
  }
  
  data <- data %>%
    mutate(
      GroupLabel = factor(paste(Sex, AgeCat, sep = "-"), levels = group_levels),
      ColorGroup = GroupLabel,
      LineType = factor(APP, levels = c("WT", "APP/+"))
    )
  
  for (strategy in unique(data$Strategy)) {
    strategy_data <- data %>% filter(Strategy == strategy)
    
    plt <- ggplot(strategy_data, aes(x = as.factor(`_Day`), y = MeanProb,
                                     group = interaction(GroupLabel, APP), color = GroupLabel, linetype = APP)) +
      geom_line(linewidth = 1) +
      geom_point(size = 2) +
      geom_errorbar(aes(ymin = MeanProb - SEM, ymax = MeanProb + SEM), width = 0.1) +
      labs(
        title = paste0(
          if (is.null(sex_filter)) {
            "All"
          } else if (sex_filter == "M") {
            "Male"
          } else if (sex_filter == "F") {
            "Female"
          } else {
            sex_filter
          },
          " ", tools::toTitleCase(strategy), " Use Over Days"
        ),
        x = "Day",
        y = "Average Probability of Strategy Use",
        color = "Sex-Age Group",
        linetype = "Genotype"
      ) +
      scale_y_continuous(limits = c(0, .5), expand = c(0, 0)) +
      scale_color_manual(values = sex_age_colors) +
      scale_linetype_manual(values = c("WT" = "solid", "APP/+" = "22")) +
      guides(
        linetype = guide_legend(order = 1),
        color = guide_legend(order = 2)
      ) +
      theme_minimal(base_size = 14) +
      theme(
        legend.box = "vertical"
      )
    
    
    print(plt)
    
    # Save to Figures folder
    if (!is.null(save_folder)) {
      safe_g <- paste0(title_prefix, strategy) %>% gsub("[^A-Za-z0-9]", "_", .)
      ggsave(file.path(save_folder, paste0(safe_g, ".jpeg")), plt, width = 10, height = 6, dpi = 300)}
  }
}

# All groups
plot_allocentric_strategies(allocentric_summary, sex_filter = NULL, title_prefix = "All_", save_folder = fig_folder)

# Females only
plot_allocentric_strategies(allocentric_summary, sex_filter = "F", title_prefix = "Female_", save_folder = fig_folder)

# Males only
plot_allocentric_strategies(allocentric_summary, sex_filter = "M", title_prefix = "Male_", save_folder = fig_folder)








# Across Days Procedural Line Plots ---------------------------------------

# Procedural strategies to long
procedural_long <- strat_sheet %>%
  select(Group, `_Day`, Sex, APP, AgeCat, all_of(procedural)) %>%
  pivot_longer(cols = all_of(procedural), 
               names_to = "Strategy", 
               values_to = "Probability")

# Summarize
procedural_summary <- procedural_long %>%
  group_by(Group, `_Day`, Strategy, Sex, APP, AgeCat) %>%
  summarize(
    MeanProb = mean(Probability, na.rm = TRUE),
    SD = sd(Probability, na.rm = TRUE),
    SEM = SD / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(GroupLabel = paste(Sex, AgeCat, sep = "-"))

plot_procedural_strategies <- function(data, sex_filter = NULL, title_prefix = "", save_folder = NULL) {
  if (!is.null(sex_filter)) {
    data <- data %>% filter(Sex == sex_filter)
  }
  
  data <- data %>%
    mutate(
      GroupLabel = factor(paste(Sex, AgeCat, sep = "-"), levels = group_levels),
      ColorGroup = GroupLabel,
      LineType = factor(APP, levels = c("WT", "APP/+"))
    )
  
  for (strategy in unique(data$Strategy)) {
    strategy_data <- data %>% filter(Strategy == strategy)
    
    plt <- ggplot(strategy_data, aes(x = as.factor(`_Day`), y = MeanProb,
                                     group = interaction(GroupLabel, APP), color = GroupLabel, linetype = APP)) +
      geom_line(linewidth = 1) +
      geom_point(size = 2) +
      geom_errorbar(aes(ymin = MeanProb - SEM, ymax = MeanProb + SEM), width = 0.1) +
      labs(
        title = paste0(
          if (is.null(sex_filter)) {
            "All"
          } else if (sex_filter == "M") {
            "Male"
          } else if (sex_filter == "F") {
            "Female"
          } else {
            sex_filter
          },
          " ", tools::toTitleCase(strategy), " Use Over Days"
        ),
        x = "Day",
        y = "Average Probability of Strategy Use",
        color = "Sex-Age Group",
        linetype = "Genotype"
      ) +
      scale_y_continuous(limits = c(0, .31), expand = c(0, 0)) +
      scale_color_manual(values = sex_age_colors) +
      scale_linetype_manual(values = c("WT" = "solid", "APP/+" = "22")) +
      guides(
        linetype = guide_legend(order = 1),
        color = guide_legend(order = 2)
      ) +
      theme_minimal(base_size = 14) +
      theme(
        legend.box = "vertical"
      )
    
    
    print(plt)
    
    # Save to Figures folder
    if (!is.null(save_folder)) {
      safe_g <- paste0(title_prefix, strategy) %>% gsub("[^A-Za-z0-9]", "_", .)
      ggsave(file.path(save_folder, paste0(safe_g, ".jpeg")), plt, width = 10, height = 6, dpi = 300)}
  }
}

# All groups
plot_procedural_strategies(procedural_summary, sex_filter = NULL, title_prefix = "All_", save_folder = fig_folder)

# Females only
plot_procedural_strategies(procedural_summary, sex_filter = "F", title_prefix = "Female_", save_folder = fig_folder)

# Males only
plot_procedural_strategies(procedural_summary, sex_filter = "M", title_prefix = "Male_", save_folder = fig_folder)






