library(tidyverse)

df <- read_excel("/Users/miasponseller/Desktop/Lab/Rtrack/MWM_results_Probe_04-06-2025.xlsx")

# Reshape data to long format
df_long <- df %>%
  pivot_longer(
    cols = c(thigmotaxis, circling, random_path, scanning,
             chaining, directed_search, corrected_search, direct_path),
    names_to = "strat",
    values_to = "probability"
  )

# Scatterplot with strategy probability vs. goal crossings
ggplot(df_long, aes(x = goal.crossings, y = probability, color = age)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ strat, scales = "free_y") +
  labs(
    title = "Strategy Probability vs. Goal Crossings by Age Group",
    x = "Goal Crossings",
    y = "Strategy Probability")


# boxplot
ggplot(df, aes(x = age, y = goal.crossings)) +
  geom_boxplot() +
  geom_jitter(aes(color = age), width = 0.2, alpha = 0.5) +
  scale_fill_manual(values = c(
    "old" = rgb(0.5, 0, 0.5),
    "young" = rgb(0, 0.6, 0)
  )) +
  scale_color_manual(values = c(
    "old" = rgb(0.5, 0, 0.5),
    "young" = rgb(0, 0.6, 0)
  )) +
  labs(title = "Goal Crossings by Age Group",
       x = "Age Group", 
       y = "Goal Crossings") +
  theme_minimal()



# boxplot with mean and error bars
df %>%
  group_by(age) %>%
  summarise(mean_gc = mean(goal.crossings),
            se = sd(goal.crossings) / sqrt(n())) %>%
  ggplot(aes(x = age, y = mean_gc, fill = age)) +
  geom_col() +
  geom_errorbar(aes(ymin = mean_gc - se, ymax = mean_gc + se), width = 0.2) +
  scale_fill_manual(values = c(
    "old" = rgb(0.5, 0, 0.5),
    "young" = rgb(0, 0.6, 0)
  )) +
  labs(title = "Average Goal Crossings by Age Group",
       x = "Age Group", 
       y = "Mean Goal Crossings") +
  theme_minimal()


# violin plot with boxplot overlay
ggplot(df, aes(x = age, y = goal.crossings, fill = age)) +
  geom_violin(trim = FALSE, alpha = 0.3) + 
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  scale_fill_manual(values = c(
    "old" = rgb(0.5, 0, 0.5),
    "young" = rgb(0, 0.6, 0)
  )) +
  scale_color_manual(values = c(
    "old" = rgb(0.5, 0, 0.5),
    "young" = rgb(0, 0.6, 0)
  )) +
  labs(title = "Goal Crossings by Age Group",
       x = "Age Group", 
       y = "Goal Crossings") +
  theme_minimal()




