library(tidyr)
library(dplyr)
library(ggplot2)
library(readxl)


fp <- read_excel("/Users/miasponseller/Desktop/Lab/Rtrack/MWM_results.xlsx")
n_young <- 28
n_old <- 57


# Strategy counts table, young and old total and across days --------------

strat_counts <- function(fp) {
  if ("_Day" %in% colnames(fp)) {
    fp <- fp %>% rename(Day = "_Day")
  }
  
  summary_df <- fp %>% 
    group_by(Age, Day, name) %>% 
    summarise(count = n(), .groups = "drop")
  
  result <- summary_df %>% 
    pivot_wider(names_from = c(Day, Age),
                values_from = count,
                names_glue = "day{Day}{Age}",
                values_fill = list(count = 0)) %>% 
    
    group_by(name) %>% 
    summarise(
      young_total = sum(c_across(starts_with("day1young"):starts_with("day4young"))),
      old_total = sum(c_across(starts_with("day1old"):starts_with("day4old"))),
      day1young = sum(day1young, na.rm = TRUE),
      day2young = sum(day2young, na.rm = TRUE),
      day3young = sum(day3young, na.rm = TRUE),
      day4young = sum(day4young, na.rm = TRUE),
      day1old = sum(day1old, na.rm = TRUE),
      day2old = sum(day2old, na.rm = TRUE),
      day3old = sum(day3old, na.rm = TRUE),
      day4old = sum(day4old, na.rm = TRUE)
    )
  return(result)
}


# Strategies across days histogram ----------------------------------------

strat_hist <- function(result) {
  long_data <- result %>% 
    pivot_longer(cols = starts_with("day"),
                 names_to = c("Day", "Age"),
                 names_pattern = "day(\\d)(\\w+)",
                 values_to = "Count") %>% 
    mutate(Age = ifelse(Age == "young", "young", "old"))
  
  young_data <- long_data %>% filter(Age == "young")
  old_data <- long_data %>% filter(Age == "old")
  
  # young and old plots for strats across days
  young_plot <- ggplot(young_data, aes(x = Day, y = Count/25, fill = name)) +
    geom_bar(stat = "identity", position = "stack") + 
    labs(title = "Young Rats Strategies Across Days",
         x = "Day",
         y = "Frequency",
         fill = "Strategy") +
    scale_fill_brewer(palette = "Set3") +
    theme_minimal()
  
  old_plot <- ggplot(old_data, aes(x = Day, y = Count/54, fill = name)) +
    geom_bar(stat = "identity", position = "stack") + 
    labs(title = "Old Rats Strategies Across Days",
         x = "Day",
         y = "Frequency",
         fill = "Strategy") +
    scale_fill_brewer(palette = "Set3") +
    theme_minimal()
  
  return(list(young_plot = young_plot, old_plot = old_plot))
}



# Strategies by Age each day ----------------------------------------------

strat_bar_chart <- function(result) {
  
  long_data <- result %>% 
    pivot_longer(
      cols = starts_with("day"),
      names_to = c("Day", "Age"),
      names_pattern = "day(\\d)(\\w+)",
      values_to = "Count"
    ) %>% 
    mutate(Day = as.factor(Day),
           Age = factor(Age, levels = c("young", "old")))
  
  # Add Count_by_num column
  long_data <- long_data %>%
    mutate(
      Count_by_num = ifelse(
        Age == "young",
        Count / n_young,
        Count / n_old
      )
    )
  
  # mean and SEM for each group
  data_summary <- long_data %>% 
    group_by(name, Day, Age) %>% 
    summarise(
      Mean = mean(Count_by_num),
      SEM = sd(Count_by_num),
      .groups = "drop"
    )
  
  # Plot each day in a separate chart
  bar_chart <- ggplot(data_summary, aes(x = Day, y = Mean, fill = Age)) +
    geom_bar(stat = "identity", position = "dodge", color = "black") +
    geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM),
                  position = position_dodge(width = .8), width = .2) +  # Fix here
    labs(title = "Strategy Frequencies by Day and Age",
         x = "Day",
         y = "Frequency",
         fill = "Age Group") +
    scale_fill_brewer(palette = "Set2") +
    theme_minimal() + 
    facet_wrap(~name, scales = "free_y")
  
  return(bar_chart)
}



# Combined chart ----------------------------------------------------------

strats_all_days_ages <- function(result) {
  long_data <- result %>%
    pivot_longer(
      cols = starts_with("day"),
      names_to = c("Day", "Age"),
      names_pattern = "day(\\d)(\\w+)",
      values_to = "Count"
    ) %>%
    mutate(
      Day = as.factor(Day), 
      Age = factor(Age, levels = c("young", "old"))
    )
  
  # Plot all data in a single chart
  combined_chart <- ggplot(long_data, aes(x = Day, y = Count, fill = name)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black") +
    facet_wrap(~ Age, nrow = 1) +
    labs(
      title = "Combined Strategy Frequencies by Day and Age Group",
      x = "Day",
      y = "Frequency",
      fill = "Strategy"
    ) +
    scale_fill_brewer(palette = "Set3") +
    theme_minimal() +
    theme(
      legend.position = "right",
      strip.background = element_rect(fill = "lightgrey", color = "black"),
      strip.text = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    )
  
  return(combined_chart)
}



# By types -----------------------------------------------------------------

by_type <- function(result) {
  
  non_goal_oriented <- c("thigmotaxis", "circling", "random path")
  procedural <- c("scanning", "chaining")
  allocentric <- c("directed search", "corrected path", "direct path")
  
  # stacked histogram by type for young and old
  
  long_data <- result %>% 
    pivot_longer(
      cols = starts_with("day"),
      names_to = c("Day", "Age"),
      names_pattern = "day(\\d)(\\w+)",
      values_to = "Count"
    ) %>% 
    mutate(
      Day = as.factor(Day),
      Age = factor(Age, levels = c("young", "old")),
      type = case_when(
        name %in% non_goal_oriented ~ "Non-goal Oriented",
        name %in% procedural ~ "Procedural",
        name %in% allocentric ~ "Allocentric"
      )
    )
  
  young_data <- long_data %>% filter(Age == "young")
  old_data <- long_data %>% filter(Age == "old")
  
  young_type_plot <- ggplot(young_data, aes(x = Day, y = Count/25, fill = type)) + 
    geom_bar(stat = "identity", position = "stack") +
    labs(title = "Young Rats Strategies by Type Across Days",
         x = "Day",
         y = "Frequency") +
    scale_fill_brewer(palette = "Set1") +
    theme_minimal()
  
  old_type_plot <- ggplot(old_data, aes(x = Day, y = Count/54, fill = type)) + 
    geom_bar(stat = "identity", position = "stack") +
    labs(title = "Old Rats Strategies by Type Across Days",
         x = "Day",
         y = "Frequency") +
    scale_fill_brewer(palette = "Set1") +
    theme_minimal()
  
  return(list(young_type_plot = young_type_plot, old_type_plot = old_type_plot))
}  



# prob_means ------------------------------------------------------------
prob_means <- function(fp) {
  
  # df with count of each strat for each day for each rat
  strat_count_df <<- fp %>%
    group_by(`_TargetID`, `_Day`, Age, name) %>%
    summarise(count = n(), .groups = "drop") %>%
    pivot_wider(names_from = name, values_from = count, values_fill = list(count = 0))
  
  # df with mean prob of each strat for each day for each rat
  avg_df <<- fp %>% 
    group_by(`_TargetID`, `_Day`, Age) %>%
    summarise(
      thigmotaxis_avg = mean(thigmotaxis, na.rm = TRUE),
      circling_avg = mean(circling, na.rm = TRUE),
      random_path_avg = mean(`random path`, na.rm = TRUE),
      chaining_avg = mean(chaining, na.rm = TRUE),
      directed_search_avg = mean(`directed search`, na.rm = TRUE),
      corrected_serach_avg = mean(`corrected search`, na.rm = TRUE),
      direct_path_avg = mean(`direct path`, na.rm = TRUE),
      perseverance_avg = mean(perseverance, na.rm = TRUE),
      .groups = "drop"
    )
  write_xlsx(list("avg_df" = avg_df), "/Users/miasponseller/Desktop/Lab/Rtrack/avg_df.xlsx")
  write_xlsx(list("strat_count_df" = strat_count_df), "/Users/miasponseller/Desktop/Lab/Rtrack/strat_count_df.xlsx")
  
}

#prob_means(fp)

# Run Functions -----------------------------------------------------------
#strat_counts(fp)
#strat_hist(strat_counts(fp))
#strat_bar_chart(strat_counts(fp))
#strats_all_days_ages(strat_counts(fp))
#by_type(strat_counts(fp))


