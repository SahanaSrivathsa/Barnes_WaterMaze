library(ggplot2)
library(tidyr)
library(dplyr)
library(car)
library(readxl)


avg_df <- read_xlsx("/Users/miasponseller/Desktop/Lab/Rtrack/avg_df.xlsx")
strat_count_df <- read_xlsx("/Users/miasponseller/Desktop/Lab/Rtrack/strat_count_df.xlsx")
mwm_results <- read_xlsx("/Users/miasponseller/Desktop/Lab/Rtrack/MWM_results.xlsx")

# Age Avg Prob of Strats Across Days --------------------------------------
avg_age_prob_strat_across_days <- function() {
  avg_long <- avg_df %>% 
    pivot_longer(cols = ends_with("_avg"), names_to = "strategy", values_to = "avg_probability")
  
  summary_df <- avg_long %>% 
    group_by(`_Day`, Age, strategy) %>%
    summarize(
      mean_avg = mean(avg_probability, na.rm = TRUE),
      sem = sd(avg_probability, na.rm = TRUE) / sqrt(sum(!is.na(avg_probability))),
      .groups = "drop"
    )
  
  ggplot(summary_df, aes(x = `_Day`, y = mean_avg, color = Age, group = Age)) +
    geom_point(size = 3, position = position_dodge(width = 0.2)) +
    geom_errorbar(aes(ymin = mean_avg - sem, ymax = mean_avg + sem),
                  width = .2, position = position_dodge(width = 0.2)) +
    scale_color_manual(values = c("young" = "green", "old" = "purple")) +
    facet_wrap(~ strategy, scales = "fixed") + # scales = free_y 
    labs(title = "Mean Strategy Probability by Day and Age",
         x = "Day",
         y = "Mean Probability",
         color = "Age") +
    theme_minimal()
}

avg_age_prob_strat_across_days()


# ANOVA -------------------------------------------------------------------

anova_on_prob_df <- function(avg_df) {
  strategy_cols <- grep("_avg$", names(avg_df), value = TRUE)
  
  anova_results <- lapply(strategy_cols, function(strategy) {
    formula <- as.formula(paste(strategy, "~ Age * `_Day`"))
    model <- aov(formula, data = avg_df)
    summary(model)
  })
  
  names(anova_results) <- strategy_cols
  anova_results
}

#anova_on_prob_df(avg_df)


# T-test ------------------------------------------------------------------


anova_on_prob_df <- function(avg_df) {
  strategy_cols <- grep("_avg$", names(avg_df), value = TRUE)
  
  results_list <- lapply(strategy_cols, function(strategy) {
    formula <- as.formula(paste(strategy, "~ Age * `_Day`"))
    model <- aov(formula, data = avg_df)
    
    # Extract ANOVA p-values
    anova_summary <- summary(model)[[1]]
    anova_pvalues <- anova_summary$`Pr(>F)`
    
    # Post-hoc Bonferroni-corrected t-tests
    pairwise_tests <- pairwise.t.test(avg_df[[strategy]], 
                                      interaction(avg_df$Age, avg_df$`_Day`), 
                                      p.adjust.method = "bonferroni")
    
    # Convert pairwise results into a data frame
    pairwise_df <- as.data.frame(as.table(pairwise_tests$p.value))
    colnames(pairwise_df) <- c("Group1", "Group2", "p_value")
    pairwise_df$Strategy <- strategy  # Add strategy name
    
    # Return a list of results
    list(
      ANOVA_Pvalues = data.frame(
        Strategy = strategy,
        Age_p = anova_pvalues[1], 
        Day_p = anova_pvalues[2], 
        Interaction_p = anova_pvalues[3]
      ),
      Pairwise_Comparisons = pairwise_df
    )
  })
  
  # Combine all ANOVA p-values into one data frame
  anova_results_df <<- do.call(rbind, lapply(results_list, `[[`, "ANOVA_Pvalues"))
  
  # Combine all pairwise comparisons into one data frame
  pairwise_results_df <<- do.call(rbind, lapply(results_list, `[[`, "Pairwise_Comparisons"))
  
  # Store results in a global variable for later access
  assign("anova_results", list(ANOVA_Results = anova_results_df, Pairwise_Comparisons = pairwise_results_df), envir = .GlobalEnv)
}


#anova_on_prob_df(avg_df)
#anova_results$Pairwise_Comparisons  # pairwise comparison










