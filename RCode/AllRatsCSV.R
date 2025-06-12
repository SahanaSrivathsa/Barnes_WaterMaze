library(readxl)
library(dplyr)
library(stringr)
library(purrr)
library(writexl)

# Folder with all ANY-Maze spreadsheets, which should include Key and Spatial sheets.
# Key columns: BarnesID or CowenID
# not necessary, but can include for Key: Age, Sex, treatment type
# Spatial columns: Test, Animal, Trial
# not necessary, but can include for Spatial: Duration, Distance, Mean speed, Path efficiency, CIPL
spatial_sheets_folder <- '/Users/miasponseller/Desktop/Lab/Rtrack/Tg/Spatial Sheets'
output_file <- '/Users/miasponseller/Desktop/Lab/Rtrack/Tg/Tg_AllRats_Spatial.csv'

# CAS ONLY - Performance sheets and ID performance dictionaries
# young = read_excel('/Users/miasponseller/Desktop/Lab/Rtrack/CAS/Performance Sheets/2023-01-18_Rats_New_age_6.xlsx')
# middle = read_excel('/Users/miasponseller/Desktop/Lab/Rtrack/CAS/Performance Sheets/2023-01-18_Rats_New_age_15.xlsx')
# old = read_excel('/Users/miasponseller/Desktop/Lab/Rtrack/CAS/Performance Sheets/2023-01-18_Rats_New_age_23.xlsx')

# young_perf_dict <- setNames(young$performance, young$ID)
# middle_perf_dict <- setNames(middle$performance, middle$ID)
# old_perf_dict <- setNames(old$performance, old$ID)
# all_perf_dict <- c(young_perf_dict, middle_perf_dict, old_perf_dict)

# List all .xlsx files in the folder
spatial_files <- list.files(spatial_sheets_folder, pattern = '\\.xlsx$', full.names = TRUE)

# Empty list to store dataframes
all_data <- list()

for (file_path in spatial_files) {
  file_name <- basename(file_path)
  
  # Get cohort number
  cohort_num <- str_match(file_name, 'Coh([^-]+)-')[,2]
  if (is.na(cohort_num)) {
    message(paste('No cohort found for', file_name))
    cohort_num <- 'Unknown'
  }
  
  # Read sheets
  tryCatch({
    key_df <- read_excel(file_path, sheet = 'Key')
    spatial_df <- read_excel(file_path, sheet = 'Spatial') %>%
      select(any_of(c('Test', 'Animal', 'Trial', 'CIPL'))) %>%
      mutate(Animal = as.character(Animal))
    
    id_cols <- intersect(c('BarnesID', 'CowenID'), colnames(key_df))
    if (length(id_cols) == 0) {
      stop(paste('No valid ID columns (BarnesID or CowenID) found in', file_name))
    }
    
    key_df <- key_df %>% select(all_of(c(id_cols, 'Sex', 'APP')))
    
    
    # Pivot to long format
    key_long <- key_df %>% 
      pivot_longer(cols = all_of(id_cols), names_to = 'IDtype', values_to = 'Animal') %>% 
      mutate(Animal = as.character(Animal))
    
    # Merge and add cohort number
    merged_df <- spatial_df %>% 
      left_join(key_long, by = 'Animal') %>% 
      mutate(Cohort = cohort_num) %>% 
      select(-IDtype) 
    
    # Get performance function (CAS)
    # get_performance <- function(animal_id, age) {
    #   age_str <- as.character(age)
    #   if (str_detect(age_str, "6")) {
    #     return(young_perf_dict[[animal_id]] %||% all_perf_dict[[animal_id]])
    #   } else if (str_detect(age_str, "15")) {
    #     return(middle_perf_dict[[animal_id]] %||% all_perf_dict[[animal_id]])
    #   } else if (str_detect(age_str, "23")) {
    #     return(old_perf_dict[[animal_id]] %||% all_perf_dict[[animal_id]])
    #   } else {
    #     return(all_perf_dict[[animal_id]])
    #   }
    # }
    
    #merged_df <- merged_df %>% 
      # rowise() %>% 
      # mutate(Performannce == get_performance(Animal, Age)) %>% 
      # ungroup()
    
  all_data[[length(all_data) + 1]] <- merged_df
  }, error = function(e) {
    message(paste('Error reading', file_name, ':', e$message))
  })
}

# Combine all data and write to file
final_df <- bind_rows(all_data)
write.csv(final_df, output_file, row.names = FALSE)

cat("All data combined and saved to", output_file, "\n")

