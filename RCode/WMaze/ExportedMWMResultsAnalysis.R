library(tidyr)
library(dplyr)
library(ggplot2)
library(readxl)

fp <- read_excel("/Users/miasponseller/Desktop/Lab/Rtrack/MWM_results_04-06-2025.xlsx")
n_young <- 37
n_old <- 66


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

# bar chart showing the strategy types avg probability each day for young and old
by_type <- function(fp) {
  
  # Strategy categories
  non_goal_oriented <- c("thigmotaxis", "circling", "random_path")
  procedural <- c("scanning", "chaining")
  allocentric <- c("directed_search", "corrected_search", "direct_path")
  
  # Convert to long format and categorize
  long_data <- fp %>% 
    pivot_longer(cols = all_of(c(non_goal_oriented, procedural, allocentric)),
                 names_to = "strat_type", values_to = "probability") %>%
    mutate(strategy_type = case_when(
      strat_type %in% non_goal_oriented ~ "Non-Goal Oriented",
      strat_type %in% procedural ~ "Procedural",
      strat_type %in% allocentric ~ "Allocentric"
    ))
  
  # Summarize mean and standard error
  df_summary <- long_data %>% 
    group_by(Age, `_Day`, strategy_type) %>%
    summarize(mean_prob = mean(probability, na.rm = TRUE), 
              se = sd(probability, na.rm = TRUE)/sqrt(n()),
              .groups = "drop") %>%
    mutate(age_strategy = paste(Age, strategy_type, sep = "_"))
  
  # Factor levels
  df_summary$Age <- factor(df_summary$Age, levels = c("young", "old"))
  df_summary$age_strategy <- factor(df_summary$age_strategy,
                                    levels = c(
                                      "young_Allocentric", "young_Procedural", "young_Non-Goal Oriented",
                                      "old_Allocentric", "old_Procedural", "old_Non-Goal Oriented"
                                    ))

  # custom colors
  color_map <- c(
    "young_Allocentric" = "#98FB98",  
    "young_Procedural" = "#32CD32",   
    "young_Non-Goal Oriented" = "#228B22", 
    "old_Allocentric" = "#b19cd9",     
    "old_Procedural" = "#9370DB",     
    "old_Non-Goal Oriented" = "#6A5ACD" 
  )
  
  # Plot
  ggplot(df_summary, aes(x = factor(`_Day`), y = mean_prob, fill = age_strategy)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_errorbar(aes(ymin = mean_prob - se, ymax = mean_prob + se),
                  position = position_dodge(width = 0.9), width = 0.2) +
    scale_fill_manual(
      values = color_map,
      labels = c(
        "young_Allocentric" = "Young Allocentric",
        "young_Procedural" = "Young Procedural",
        "young_Non-Goal Oriented" = "Young Non-Goal Oriented",
        "old_Allocentric" = "Old Allocentric",
        "old_Procedural" = "Old Procedural",
        "old_Non-Goal Oriented" = "Old Non-Goal Oriented"
      )
    ) +
    labs(title = "Average Strategy Type Probability by Day and Age Group", 
         x = "Day", y = "Mean Probability", fill = "Strategy Type") +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16),
      plot.title = element_text(size = 18, face = "bold"),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12)
    )
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
      random_path_avg = mean(random_path, na.rm = TRUE),
      chaining_avg = mean(chaining, na.rm = TRUE),
      directed_search_avg = mean(directed_search, na.rm = TRUE),
      corrected_serach_avg = mean(corrected_search, na.rm = TRUE),
      direct_path_avg = mean(direct_path, na.rm = TRUE),
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
by_type(fp)


