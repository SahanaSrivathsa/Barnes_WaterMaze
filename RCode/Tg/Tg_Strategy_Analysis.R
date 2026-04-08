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
library(rstatix)
library(car)

# MUST RUN ----------------------------------------------------------------

# Load data
strat_sheet <- read_excel('/Users/miasponseller/Desktop/Lab/Rtrack/Tg/Tg_MWM_results_04-04-2026.xlsx')
all_rats_spatial <- read.csv('/Users/miasponseller/Desktop/Lab/Rtrack/Tg/Tg_AllRats_Spatial_cleaned.csv')
description_file <- '/Users/miasponseller/Desktop/Lab/Rtrack/Tg/Tg_exp_desc.xlsx'

strat_sheet$Age = as.numeric(strat_sheet$Age)

# Add if running Rtrack plots
all_trials <- '/Users/miasponseller/Desktop/Lab/Rtrack/Tg/All Tg Tracks'

# Fix floating-point precision issue for Age column
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
platform_independent <- c('thigmotaxis', 'circling', 'random path')
procedural <- c('scanning', 'chaining')
allocentric <- c('directed search', 'corrected path', 'direct path')

list_of_strats = c('direct path', 'corrected path', 'directed search', 'chaining', 
                   'scanning', 'random path', 'circling', 'thigmotaxis')

list_of_strat_cats = c('Allocentric', 'Procedural', 'PlatformIndependent')

# Add StratCat, AgeCat, AgeGroup, and Group columns
strat_sheet <- strat_sheet %>% 
  mutate(
    StratCat = case_when(
      name %in% platform_independent ~ 'PlatformIndependent',
      name %in% procedural ~ 'Procedural',
      name %in% allocentric ~ 'Allocentric'
    )) %>% 
  mutate(
    AgeCat = case_when(
      Age < 7  ~ '5',
      Age >= 7 & Age < 10 ~ '8', 
      Age >= 10 & Age < 13 ~ '11',
      Age >= 13 & Age < 15 ~ '13.5',
      Age>= 15 & Age < 16 ~ '15.5',
      Age >= 20 ~ '20'
    )) %>% 
  mutate(Group = paste0(Sex, '-', APP)) %>% 
  mutate(
    AgeGroup = case_when(
      Age < 9 ~ 'Young',
      Age >= 9 & Age < 15 ~ 'Middle',
      Age >= 15 ~ 'Old'
    )
  )

# StratCat levels
strat_sheet$StratCat <- factor(
  strat_sheet$StratCat, levels = c('PlatformIndependent', 'Procedural', 'Allocentric')
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

# Add Allocentric, Procedural, and PlatformIndependent columns with avg prob for trial
strat_sheet$PlatformIndependentProb <- rowSums(strat_sheet[,platform_independent])
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
  'PlatformIndependent' = '#654321'
)

sex_age_colors <- c(
  "F-5"     = "#e6194b",
  "M-5"     = "#f58231",
  "F-8"     = "#3cb44b",
  "M-8"     = "#0082c8",
  "F-11"    = "#911eb4",
  "M-11"    = "#46f0f0",
  "F-15.5"  = "#000075",
  "M-15.5"  = "#f032e6",
  "F-20"    = "#bfef45",
  "M-20"    = "#7f7fff"
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

r_track_strat_plots <- function() {
  bulk_strategy_calling <- function() {
    experiment_file <- description_file
    
    experiment <<- Rtrack::read_experiment(experiment_file, data.dir = all_trials)
    strategies <<- Rtrack::call_strategy(experiment)
    list(experiment = experiment, strategies = strategies)
  }
  
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
  Rtrack::plot_strategies(strategies, experiment = experiment, factor = "Genotype")
  
  print("Plots successfully created")
}
bulk_strategy_calling <- function() {
  experiment_file <- description_file
  
  experiment <<- Rtrack::read_experiment(experiment_file, data.dir = all_trials)
  strategies <<- Rtrack::call_strategy(experiment)
  list(experiment = experiment, strategies = strategies)
}

r_track_strat_plots()

# Strategy Use Proportions by Rat -----------------------------------------

strategy_proportions_graphs <- function() {
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
}

strategy_proportions_graphs()

# Prob Strat Cat Sex/Genotype ---------------------------------------------

sex_genotype <- function() {
  # Average strategy category by rat, day, and group
  rat_day_avg <- strat_sheet %>% 
    group_by(Group, `_TargetID`, `_Day`) %>% 
    summarize(
      Allocentric = mean(AllocentricProb, na.rm = TRUE),
      Procedural = mean(ProceduralProb, na.rm = TRUE),
      PlatformIndependent = mean(PlatformIndependentProb, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Pivot long
  rat_day_long <- rat_day_avg %>% 
    pivot_longer(cols = c(Allocentric, Procedural, PlatformIndependent),
                 names_to = 'StratCat',
                 values_to = 'Probability') %>% 
    mutate(StratCat = factor(
      StratCat,
      levels = c('PlatformIndependent', 'Procedural', 'Allocentric')
    ))
  
  group_day_summary <- rat_day_long %>% 
    group_by(Group, `_Day`, StratCat) %>% 
    summarize(MeanProb = mean(Probability, na.rm = TRUE),
              SD = sd(Probability, na.rm = TRUE),
              SEM = SD/sqrt(n()),
              .groups = 'drop')
  
  for (g in unique(group_day_summary$Group)) {
    
    # Filter data for group g
    summary_data <- group_day_summary %>% filter(Group == g)
    rat_data <- rat_day_long %>% filter(Group == g)
    
    # Make _Day a factor with all 4 days explicitly
    rat_data <- rat_data %>%
      mutate(
        DayFactor = factor(`_Day`, levels = c("1", "2", "3", "4"))
      )
    
    # Print counts of observations per StratCat and Day
    cat("\n=== Counts of observations by StratCat and Day for Group:", g, "===\n")
    print(
      rat_data %>% 
        group_by(StratCat, DayFactor) %>% 
        tally() %>% 
        arrange(StratCat, DayFactor)
    )
    
    # Run pairwise t-tests within each strategy category (all pairs possible)
    pvals_df_all <- rat_data %>%
      group_by(StratCat) %>%
      pairwise_t_test(
        Probability ~ DayFactor,
        p.adjust.method = "bonferroni"
      ) %>%
      rename(
        DayNum1 = group1,
        DayNum2 = group2
      )
    
    # Print all day pairs tested
    cat("\n=== All day pairs being compared for Group:", g, "===\n")
    pvals_df_all %>%
      select(StratCat, DayNum1, DayNum2) %>%
      distinct() %>%
      arrange(StratCat, DayNum1, DayNum2) %>%
      print(n = Inf)
    
    # Filter significant pairs and add labels
    pvals_df <- pvals_df_all %>%
      mutate(
        p.adj.label = case_when(
          p.adj < 0.001 ~ "***",
          p.adj < 0.01  ~ "**",
          p.adj < 0.05  ~ "*",
          TRUE ~ NA_character_
        )
      ) %>%
      filter(!is.na(p.adj.label)) %>%
      mutate(
        y.position = case_when(
          StratCat == "PlatformIndependent" ~ 0.8,
          StratCat == "Procedural" ~ 0.9,
          StratCat == "Allocentric" ~ 1.0
        ) + as.numeric(factor(paste(DayNum1, DayNum2, sep = "_"))) * 0.025
      )
    
    # Calculate dodge offsets
    dodge_width <- 0.9
    strat_levels <- c('PlatformIndependent', 'Procedural', 'Allocentric')
    n_strat <- length(strat_levels)
    bar_width <- dodge_width / n_strat
    dodge_offsets <- setNames(
      seq(0, by = bar_width, length.out = n_strat) - dodge_width/2 + bar_width/2,
      strat_levels
    )
    
    # Ensure DayNum1 and DayNum2 factors have correct levels
    levels_day <- levels(rat_data$DayFactor)
    
    pvals_df <- pvals_df %>%
      mutate(
        DayNum1 = factor(DayNum1, levels = levels_day),
        DayNum2 = factor(DayNum2, levels = levels_day),
        DayNum1_num = as.numeric(DayNum1),
        DayNum2_num = as.numeric(DayNum2),
        xmin = DayNum1_num + dodge_offsets[StratCat],
        xmax = DayNum2_num + dodge_offsets[StratCat]
      )
    
    # Print significant p-values
    cat("\n=== Significant pairwise p-values for Group:", g, "===\n")
    print(pvals_df %>% select(StratCat, DayNum1, DayNum2, p, p.adj, p.adj.label))
    
    # Plot
    plt <- ggplot(summary_data, aes(x = as.factor(`_Day`), y = MeanProb, fill = StratCat)) +
      geom_bar(stat = "identity", position = position_dodge(width = dodge_width), width = .8, alpha = .5) +
      geom_jitter(data = rat_data,
                  aes(x = as.factor(`_Day`), y = Probability, color = StratCat),
                  position = position_jitterdodge(jitter.width = .2, dodge.width = dodge_width),
                  alpha = .8, size = 2, inherit.aes = FALSE) +
      geom_errorbar(aes(ymin = MeanProb - SEM, ymax = MeanProb + SEM),
                    width = .2, position = position_dodge(width = dodge_width)) +
      stat_pvalue_manual(
        data = pvals_df %>%
          rename(
            group1 = DayNum1,
            group2 = DayNum2
          ),
        label = "p.adj.label",
        y.position = "y.position",
        xmin = "xmin",
        xmax = "xmax",
        tip.length = 0.01,
        step.increase = 0.02,
        bracket.size = 0.3
      ) +
      labs(
        title = paste("Strategy Use:", g),
        x = "Day",
        y = "Probability of Strategy Use",
        fill = "Strategy Category",
        color = 'Strategy Category'
      ) +
      scale_fill_manual(values = strat_cat_colors) +
      scale_color_manual(values = strat_cat_colors) +
      coord_cartesian(ylim = c(0, 1.1)) +
      theme_minimal(base_size = 14)
    
    print(plt)
    
    # Save plot to Figures folder
    # safe_g <- gsub("[^A-Za-z0-9_]", "_", g)
    # file_path <- file.path(fig_folder, paste0("AvgStratCatUse_Day_Sex_Genotype_", safe_g, ".jpeg"))
    # ggsave(file_path, plt, width = 8, height = 6, dpi = 300)
    
  }
}

sex_genotype()

# Prob StratCat Sex/Genotype/AgeCat ---------------------------------------

sex_genoetype_agecat <- function() {
  # Average strategy category by rat, day, sex, genotype, and age category
  rat_day_avg2 <- strat_sheet %>%
    group_by(Sex, APP, AgeCat, `_TargetID`, `_Day`) %>%
    summarize(
      Allocentric = mean(AllocentricProb, na.rm = TRUE),
      Procedural = mean(ProceduralProb, na.rm = TRUE),
      PlatformIndependent = mean(PlatformIndependentProb, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Pivot long
  rat_day_long2 <- rat_day_avg2 %>%
    pivot_longer(
      cols = c(Allocentric, Procedural, PlatformIndependent),
      names_to = 'StratCat',
      values_to = 'Probability'
    ) %>%
    mutate(StratCat = factor(
      StratCat,
      levels = c('PlatformIndependent', 'Procedural', 'Allocentric')
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
}

sex_genoetype_agecat()

# Across Days StratCat Line Plots -----------------------------------------

rat_day_avg2 <- strat_sheet %>%
  group_by(Sex, APP, AgeCat, `_TargetID`, `_Day`) %>%
  summarize(
    Allocentric = mean(AllocentricProb, na.rm = TRUE),
    Procedural = mean(ProceduralProb, na.rm = TRUE),
    PlatformIndependent = mean(PlatformIndependentProb, na.rm = TRUE),
    .groups = 'drop'
  )

# Pivot long
rat_day_long2 <- rat_day_avg2 %>%
  pivot_longer(
    cols = c(Allocentric, Procedural, PlatformIndependent),
    names_to = 'StratCat',
    values_to = 'Probability'
  ) %>%
  mutate(StratCat = factor(
    StratCat,
    levels = c('PlatformIndependent', 'Procedural', 'Allocentric')
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


# ANOVA: Sex x AgeGroup x Genotype (excluding old) ------------------------

sex_agegroup_genotype_anova <- function() {
  
  fp <- '/Users/miasponseller/Desktop/Lab/Rtrack/Tg/tg_anova_no_old_day4.txt'
  # Open sink to save all console output to a file
  sink(fp)
  
  # Filter and prepare data as before
  rat_strategy_avg <- strat_sheet %>%
    group_by(`_TargetID`, StratCat, Sex, AgeGroup, APP) %>%
    filter(AgeGroup %in% c('Young', 'Middle')) %>% # Exclude old age group
    filter(`_Day` == 4) %>% # Filter to a specific day
    summarize(
      MeanProb = mean(case_when(
        StratCat == "Allocentric" ~ AllocentricProb,
        StratCat == "Procedural" ~ ProceduralProb,
        StratCat == "PlatformIndependent" ~ PlatformIndependentProb,
        TRUE ~ NA_real_
      ), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(TargetID = `_TargetID`) %>%
    mutate(
      AgeGroup = factor(AgeGroup),
      Sex = factor(Sex),
      APP = factor(APP),
      StratCat = factor(StratCat, levels = c("PlatformIndependent", "Procedural", "Allocentric"))
    )
  
  # Run ANOVA loop, output goes to file
  for (strat in levels(rat_strategy_avg$StratCat)) {
    cat("\n\n### ANOVA for strategy category:", strat, "###\n")
    
    df <- rat_strategy_avg %>% filter(StratCat == strat)
    
    anova_res <- aov_ez(
      id = "TargetID",
      dv = "MeanProb",
      data = df,
      between = c("Sex", "AgeGroup", "APP")
    )
    
    print(anova_res$anova_table)
    
    emm <- emmeans(anova_res, ~ Sex * AgeGroup * APP)
    print(emm)
    
    cat("\nPairwise comparisons by AgeGroup within Sex × APP:\n")
    print(pairs(emm, by = c("Sex", "APP")))
    
    cat("\nPairwise comparisons by Sex within AgeGroup × APP:\n")
    print(pairs(emm, by = c("AgeGroup", "APP")))
    
    cat("\nPairwise comparisons by APP within Sex × AgeGroup:\n")
    print(pairs(emm, by = c("Sex", "AgeGroup")))
  }
  
  sink()
  message('ANOVA results saved to: ', fp)
}

sex_agegroup_genotype_anova()


# ANOVA: AgeGroup x Genotype ----------------------------------------------

agegroup_genotype_anova <- function() {
  fp <- "/Users/miasponseller/Desktop/Lab/Rtrack/tg_anova_no_sex.txt"
  
  sink(fp)
  
  rat_strategy_avg <- strat_sheet %>%
    group_by(`_TargetID`, StratCat, Sex, AgeGroup, APP) %>%
    summarize(
      MeanProb = mean(case_when(
        StratCat == "Allocentric" ~ AllocentricProb,
        StratCat == "Procedural" ~ ProceduralProb,
        StratCat == "PlatformIndependent" ~ PlatformIndependentProb,
        TRUE ~ NA_real_
      ), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(TargetID = `_TargetID`) %>%
    mutate(
      AgeGroup = factor(AgeGroup),
      APP = factor(APP),
      StratCat = factor(StratCat, levels = c("PlatformIndependent", "Procedural", "Allocentric"))
    )
  
  # Run ANOVA loop
  for (strat in levels(rat_strategy_avg$StratCat)) {
    cat("\n\n### ANOVA for strategy category:", strat, " (excluding Sex) ###\n")
    
    df <- rat_strategy_avg %>% filter(StratCat == strat)
    
    anova_res <- aov_ez(
      id = "TargetID",
      dv = "MeanProb",
      data = df,
      between = c("AgeGroup", "APP")
    )
    
    print(anova_res$anova_table)
    
    emm <- emmeans(anova_res, ~ AgeGroup * APP)
    print(emm)
    
    cat("\nPairwise comparisons by AgeGroup within APP:\n")
    print(pairs(emm, by = "APP"))
    
    cat("\nPairwise comparisons by APP within AgeGroup:\n")
    print(pairs(emm, by = "AgeGroup"))
  }
  
  sink()
  message("ANOVA results saved to: ", fp)
}

agegroup_genotype_anova()


# Day vs. Prob StratCat Use - Sex/Genotype --------------------------------

plot_strategy_probabilities <- function() {
  long_data <- strat_sheet %>% 
    filter(as.numeric(strat_sheet$Age) < 19) %>% 
    pivot_longer(cols = c(PlatformIndependentProb, ProceduralProb, AllocentricProb),
                 names_to = "StrategyCategory",
                 values_to = "Probability")
  
  # average across trials for each rat
  rat_avg <- long_data %>%
    group_by(`_TargetID`, Group, `_Day`, StrategyCategory) %>%
    summarize(rat_mean_prob = mean(Probability, na.rm = TRUE), .groups = "drop")
  
  # group means and SE (across rats)
  plot_data <- rat_avg %>%
    group_by(Group, `_Day`, StrategyCategory) %>%
    summarize(mean_prob = mean(rat_mean_prob, na.rm = TRUE),
              se = sd(rat_mean_prob, na.rm = TRUE) / sqrt(n()),
              n_rats = n(),
              .groups = "drop")
  
  ggplot(plot_data, aes(x = factor(`_Day`),
                        y = mean_prob,
                        fill = Group)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(aes(ymin = mean_prob - se, ymax = mean_prob + se),
                  position = position_dodge(width = 0.8), width = 0.2) +
    facet_wrap(~ StrategyCategory, ncol = 1) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(x = "Day",
         y = "Average Probability of Strategy Use (per rat)",
         fill = "Group (Sex × Genotype)") +
    theme_bw(base_size = 14) +
    theme(strip.background = element_rect(fill = "gray90", color = "gray60"),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank())
}

plot_strategy_probabilities()


# Age vs. Prob StratCat Use Line Plots ------------------------------------

# Day 4 mean PlatformIndependentProb by Age / Sex / APP -------------------

platformindependent_prob_day4 <- strat_sheet %>% 
  filter(`_Day` == 4) %>% 
  group_by(Age, Sex, APP) %>% 
  summarize(
    prob = mean(PlatformIndependentProb, na.rm = TRUE),
    sd = sd(PlatformIndependentProb, na.rm = TRUE),
    n = sum(!is.na(PlatformIndependentProb)),
    se = sd/sqrt(n),
    .groups = "drop"
  ) %>% 
  mutate(Group = paste(Sex, APP, sep = "_"))

ggplot(platformindependent_prob_day4,
       aes(x = Age, y = prob,
           color = Group,
           group = Group)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = prob - se, ymax = prob + se), width = .2) + 
  labs(
    title = "Day 4 Mean Platform-Independent Probability by Age",
    x = "Age",
    y = "Mean Probability",
    color = "Sex / Genotype"
  ) +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 1))

# Day 4 mean ProceduralProb by Age / Sex / APP -------------------

procedural_prob_day4 <- strat_sheet %>% 
  filter(`_Day` == 4) %>% 
  group_by(Age, Sex, APP) %>% 
  summarize(
    prob = mean(ProceduralProb, na.rm = TRUE),
    sd = sd(ProceduralProb, na.rm = TRUE),
    n = sum(!is.na(ProceduralProb)),
    se = sd/sqrt(n),
    .groups = "drop"
  ) %>% 
  mutate(Group = paste(Sex, APP, sep = "_"))

ggplot(procedural_prob_day4,
       aes(x = Age, y = prob,
           color = Group,
           group = Group)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = prob - se, ymax = prob + se), width = .2) +
  labs(
    title = "Day 4 Mean Procedural Probability by Age",
    x = "Age",
    y = "Mean Probability",
    color = "Sex / Genotype"
  ) +
  theme_() +
  scale_y_continuous(limits = c(0, 1))

# Day 4 mean AllocentricProb by Age / Sex / APP -------------------

allocentric_prob_day4 <- strat_sheet %>% 
  filter(`_Day` == 4) %>% 
  group_by(Age, Sex, APP) %>% 
  summarize(
    prob = mean(AllocentricProb, na.rm = TRUE),
    sd = sd(AllocentricProb, na.rm = TRUE),
    n = sum(!is.na(AllocentricProb)),
    se = sd/sqrt(n),
    .groups = "drop"
  ) %>% 
  mutate(Group = paste(Sex, APP, sep = "_"))

ggplot(allocentric_prob_day4,
       aes(x = Age, y = prob,
           color = Group,
           group = Group)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = prob - se, ymax = prob + se), width = .2) +
  labs(
    title = "Day 4 Mean Allocentric Probability by Age",
    x = "Age",
    y = "Mean Probability",
    color = "Sex / Genotype"
  ) +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 1))


# Single Trial Plot -------------------------------------------------------

arena_single = Rtrack::read_arena('/Users/miasponseller/Desktop/Lab/Rtrack/Tg/Arena Files/Cohort1MArena.txt')
path_single = Rtrack::read_path('/Users/miasponseller/Desktop/Lab/Rtrack/Tg/All Tg Tracks/Coh1M_Trial205.csv', arena_single, id = 'test205', track.format = 'anymaze.csv')
metrics_single = Rtrack::calculate_metrics(path_single, arena_single)
Rtrack::plot_path(metrics_single)


# Plot All Paths (Save as PDF) --------------------------------------------

plot_all_paths <- function(output_folder = "/Users/miasponseller/Desktop/Lab/Rtrack/Tg/Plot_PDFs") {
  
  # Create output folder if it doesn't exist
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
    message("Created directory: ", output_folder)
  }
  
  desc_file <- "/Users/miasponseller/Desktop/Lab/Rtrack/Tg/Tg_exp_desc.xlsx"
  data_dir <- "/Users/miasponseller/Desktop/Lab/Rtrack/Tg/All Tg Tracks"
  
  # Read experiment
  experiment <- tryCatch({
    Rtrack::read_experiment(desc_file, data.dir = data_dir)
  }, error = function(e) {
    stop("Error reading experiment: ", e$message)
  })
  
  # Get strategies
  strategies <- tryCatch({
    call_strategy(experiment)
  }, error = function(e) {
    warning("Error in call_strategy: ", e$message)
    return(NULL)
  })
  
  if (is.null(strategies)) {
    message("No strategies found. Exiting.")
    return(NULL)
  }
  
  pdf_file <- file.path(output_folder, "Tg_all_paths.pdf")
  
  # Open PDF
  pdf(file = pdf_file)
  message("Saving plots to: ", pdf_file)
  
  # Plot each path
  for (i in seq_along(experiment$metrics)) {
    tryCatch({
      plot_path(
        experiment$metrics[[i]],
        title = paste0(
          experiment$metrics[[i]]$id,
          " - ",
          strategies$calls[i, "name"]
        )
      )
    }, error = function(e) {
      warning(paste("Error plotting metric", i, ":", e$message))
    })
  }
  
  message("Finished saving plots.")
}

plot_all_paths()




# for age 5mo (day4), pool plot w/ trajectory 
#  start w male, then females

# CIPL vs prob scatterplots (MATLAB)
# modify function to also add CIPL to plots




