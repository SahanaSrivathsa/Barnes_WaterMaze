library(tidyr)
library(dplyr)
library(ggplot2)
library(readxl)
library(scales)

# Load data
strat_sheet <- read_excel('/Users/miasponseller/Desktop/Lab/Rtrack/Tg/rats_removed_Tg_MWM_results_06-11-2025.xlsx')
all_rats_spatial <- read.csv('/Users/miasponseller/Desktop/Lab/Rtrack/Tg/Tg_AllRats_Spatial_test.csv')

# Figures Folder
fig_folder = '/Users/miasponseller/Desktop/Lab/Rtrack/Tg/Figures'
dir.create(fig_folder, recursive = TRUE, showWarnings = FALSE)

sex_counts <- all_rats_spatial %>% 
  distinct(Animal, .keep_all = TRUE) %>% 
  count(Sex)
#print(sex_counts)

n_male = sex_counts$n[sex_counts$Sex == 'M']
n_female = sex_counts$n[sex_counts$Sex == 'F']

# Sex groups
sex_groups <- c('male', 'female')

# Strategy categories
non_goal_oriented <- c('thigmotaxis', 'circling', 'random path')
procedural <- c('scanning', 'chaining')
allocentric <- c('directed search', 'corrected path', 'direct path')

list_of_strats = c('direct path', 'corrected path', 'directed search', 'chaining', 
                   'scanning', 'random path', 'circling', 'thigmotaxis')

list_of_strat_cats = c('Allocentric', 'Procedural', 'NonGoalOriented')

# Add StratCat column
strat_sheet_1 <- strat_sheet %>% 
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


# Dominant Strat Across Days M/F ------------------------------------------

# List of unique sexes
sexes <- unique(strat_sheet_1$Sex)

for (s in sexes) {
  
  sex_data <- strat_sheet_1 %>% 
    filter(Sex == s)
  
  # Mean and SE
  summary_df <- sex_data %>% 
    group_by(`_Day`, StratCat) %>% 
    summarize(
      mean_prob = mean(confidence, na.rm = TRUE),
      se = sd(confidence, na.rm = TRUE)/sqrt(n()),
      .groups = 'drop'
    )
  
  plt <- ggplot(summary_df, aes(x = factor(`_Day`), y = mean_prob, fill = StratCat)) +
    geom_bar(stat = 'identity', position = position_dodge(width = .9), alpha = .5) +
    geom_jitter(data = sex_data,
                aes(x = factor(`_Day`), y = confidence, color = StratCat),
                position = position_jitterdodge(jitter.width = .2, dodge.width = .9),
                size = 2, alpha = .7,
                inherit.aes = FALSE,
                show.legend = FALSE) +
    geom_errorbar(aes(ymin = mean_prob - se, ymax = mean_prob + se), 
                  position = position_dodge(width = .9),
                  width = .2) +
    labs(
      title = paste('Strategy Category Usage by Day -', tools::toTitleCase(s)),
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
  # safe_g <- gsub("[^A-Za-z0-9_]", "_", g)
  # file_path <- file.path(fig_folder, paste0("Dominant_Strat_Cat_Days_Sex", safe_g, ".jpeg"))
  # ggsave(file_path, plt, width = 8, height = 6, dpi = 300)
  
}

