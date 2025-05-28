library(tidyr)
library(dplyr)
library(ggplot2)
library(readxl)

# Load data
strat_sheet <- read_excel('/Users/miasponseller/Desktop/Lab/Rtrack/CAS/CAS_MWM_results_05-24-2025.xlsx')
all_rats_spatial <- read.csv('/Users/miasponseller/Desktop/Lab/Rtrack/CAS/CAS_AllRats_Spatial.csv')

age_counts <- all_rats_spatial %>% 
  distinct(Animal, .keep_all = TRUE) %>% 
  count(Age)
print(age_counts)

n_young = age_counts$n[age_counts$Age == '6mo']
n_middle = age_counts$n[age_counts$Age == '15mo']
n_old = age_counts$n[age_counts$Age == '23mo']

# Strategy categories
non_goal_oriented <- c("thigmotaxis", "circling", "random_path")
procedural <- c("scanning", "chaining")
allocentric <- c("directed_search", "corrected_search", "direct_path")

# Dominant Strat Across Days for Age/Perf Groups --------------------------

strat_sheet <- strat_sheet %>%
  mutate(
    `_TargetID` = as.character(`_TargetID`),
    `_Trial` = as.integer(`_Trial`)
  )

all_rats_spatial <- all_rats_spatial %>%
  mutate(
    Animal = as.character(Animal),
    Trial = as.integer(Trial)
  )

merged <- strat_sheet %>%
  left_join(
    all_rats_spatial %>%
      select(Animal, Trial, Age, Performance),
    by = c(`_TargetID` = "Animal", `_Trial` = "Trial")
  ) %>%
  mutate(AgePerformance = paste(Age, Performance))

# calculate proportions
strat_props <- merged %>% 
  group_by(AgePerformance, Day, name) %>% 
  summarise(rats_with_strat = n_distinct(Animal), .groups = 'drop') %>% 
  group_by(AgePerformance, Day) %>% 
  mutate(tot_rats = sum(rats_with_strat), 
         prop = rats_with_strat / tot_rats) %>% 
  ungroup()

# Plot
ggplot(strat_props, aes(x = name, y = prop, fill = name)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  facet_grid(AgePerformance ~ Day) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Proportion of Rats Using Each Strategy",
    x = "Strategy",
    y = "Proportion",
    fill = "Strategy")
