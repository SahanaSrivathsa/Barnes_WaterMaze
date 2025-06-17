library(readr)
library(dplyr)
library(writexl)

# Load CSV file (From AllRatsCSV.R)
input_file <- '/Users/miasponseller/Desktop/Lab/Rtrack/Tg/Tg_AllRats_Spatial_cleaned.csv'
output_file <- '/Users/miasponseller/Desktop/Lab/Rtrack/Tg/Tg_exp_desc.xlsx'


df <- read_csv(input_file)

# Create new columns
df <- df %>% 
  mutate(
    `_TrackID` = paste0('Coh', Cohort, '_test', Test),
    `_TargetID` = Animal,
    `_Trial` = Trial,
    `_Day` = pmin(((Trial - 1) %/% 6 + 1), 4),
    `_TrackFileFormat` = 'anymaze.csv',
    `_Arena` = paste0('Arena Files/Cohort', Cohort, 'Arena.txt'),
    `_TrackFile` = paste0('Coh', Cohort, '_Trial', Test, '.csv')
  )

expected_trials <- 1:24

# Identify animals with 24 trials
trial_sets <- df %>% 
  group_by(`_TargetID`) %>% 
  summarize(trials = list(unique(na.omit(`_Trial`)))) %>% 
  mutate(is_valid = sapply(trials, function(x) setequal(sort(x), expected_trials)))

valid_animals <- trial_sets %>% filter(is_valid) %>% pull(`_TargetID`)
invalid_animals <- trial_sets %>% filter(!is_valid)

# Filter to only animals with 24 trials (valid)
filtered_df <- df %>% filter(`_TargetID` %in% valid_animals)

# Select output columns
output_columns <- c(
  '_TrackID', '_TargetID', '_Trial', '_Day', '_TrackFileFormat', 'Sex', 'Cohort', 
  '_Arena', '_TrackFile', 'APP', 'CIPL', 'Age'
)

output_df <- filtered_df %>% select(all_of(output_columns))

# Save to Excel file
write_xlsx(output_df, output_file)
cat(paste0('Experiment description file saved to ', output_file, '\n'))

# Print invalid animals
output <- character(nrow(invalid_animals))
for (i in seq_len(nrow(invalid_animals))) {
  rat <- invalid_animals$`_TargetID`[i]
  trials <- sort(unlist(invalid_animals$trials[i]))
  output[i] <- paste0('Rat ', rat, ' Trials: ', paste(trials, collapse = ', '))
}
cat('\nAnimals excluded (missing trials 1-24 or trial labeled incorrectly):\n')
cat(paste(output, collapse = '\n'), '\n')






