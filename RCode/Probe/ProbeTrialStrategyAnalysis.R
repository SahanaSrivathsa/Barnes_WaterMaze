library(readxl)
library(tidyverse)
library(ggplot2)
library(ggsignif)

prob_sheet <- read.xlsx("/Users/miasponseller/Desktop/Lab/Rtrack/MWM_results_Probe_04-03-2025.xlsx")
all_probe <- read.csv("/Users/miasponseller/Desktop/Lab/Rtrack/AllMorrisWaterMazeData_Probe.csv")

prob_sheet$age <- factor(prob_sheet$age, levels = c("young", "old"))

strategies <- c("perseverance", "circling", "perseverance", "scanning", "perseverance", "perseverance", "perseverance", "perseverance")

# Box plot of strat prob by age group -------------------------------------

ggplot(prob_sheet, aes(x = age, y = perseverance)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = age), width = 0.2, size = 2) +
  scale_fill_manual(values = c("young" = rgb(0, 0.6, 0), "old" = rgb(0.5, 0, 0.5))) +
  scale_color_manual(values = c("young" = rgb(0, 0.6, 0), "old" = rgb(0.5, 0, 0.5))) +
  theme_minimal() +
  labs(title = "Perseverance Probability on Probe Trial",
       x = "Age Group",  y = "Perseverance Probability") +
  geom_signif(comparisons = list(c("young", "old")), 
              map_signif_level = TRUE)


# strat prob vs. CIPL -----------------------------------------------------

# ggplot(prob_sheet, aes(x = all_probe$Platform_CIPL[match(`_TargetID`, all_probe$Animal)],
#                        y = perseverance,
#                        color = age)) +
#   geom_point(size = 3, alpha = .7) +
#   geom_smooth(method = "lm", se = FALSE) +
#   scale_color_manual(values = c("young" = rgb(0, 0.6, 0), "old" = rgb(0.5, 0, 0.5))) +  # Custom colors
#   theme_minimal() +
#   labs(title = "Perseverance Probability vs. Platform_CIPL",
#        x = "CIPL",
#        y = "Perseverance Probability",
#        color = "Age Group")

# t-test
# t_test_results <- list()
# 
# for (strat in strategies) {
#   t_test_results[[strat]] <- t.test(prob_sheet[[strat]] ~ prob_sheet$age)
#   print(paste("T-test for", strat))
#   print(t_test_results[[strat]])
# }
