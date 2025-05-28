library(openxlsx)
library(readxl)
library(Rtrack)

rtrack_folder <- "/Users/miasponseller/Desktop/Lab/Rtrack/CAS" # Rtrack folder file path (main directory)
all_rats_desc_fp <- file.path('/Users/miasponseller/Desktop/Lab/Rtrack/CAS/CAS_exp_desc.xlsx') # all rats experiment description file, must be an .xlsx
trials_dir <- '/Users/miasponseller/Desktop/Lab/Rtrack/CAS/All CAS Tracks' # folder with all cohorts' trials path data
strat_export_fp <- paste0("/Users/miasponseller/Desktop/Lab/Rtrack/CAS/CAS_MWM_results_", format(Sys.Date(), "%m-%d-%Y"), ".xlsx") # output file path

cohort_list <- c('1O', '1YM', '2M', '2O', '2Y', '3M', '3O', '3Y', '4M','4O', '4Y', 
                 '5M', '5O', '5Y', '6M', '6O', '6Y', '7M', '7O', '7Y', '8M', '8O', 
                 '8Y', '9M', '9O', '10O', '10YM', '11M', '11O', '11Y', '12O', 
                 '12YM', '13OP', '13OR', '13YM', '14BR', '14P', '14RP', '14Bk', 
                 '15BR', '15G', '16B', '16G', '16R')

bulk_strategy_calling <- function() {
  experiment <<- Rtrack::read_experiment(all_rats_desc_fp, data.dir = trials_dir)
  strategies <<- Rtrack::call_strategy(experiment)
  list(experiment = experiment, strategies = strategies)
}

strat_export <- function() {
  if (!exists("experiment", envir = .GlobalEnv) || !exists("strategies", envir = .GlobalEnv)) {
    bulk_strategy_calling()
  }

  Rtrack::export_results(experiment, strategies, file = strat_export_fp)
  
  # to export a subset, specify the indices or names of tracks you want to export
  #old = experiment$factors$Age == "old"
  #Rtrack::export_results(experiment, tracks = old, file = strat_export_fp)
  
  # export tracks above threshold strategy call
  #thresholded = Rtrack::threshold_strategies(strategies, threshold = 0.2)
  #Rtrack::export_results(experiment, strategies, tracks = rownames(thresholded$calls), file = strat_export_fp)
  
  
  # rename numbered strategy columns to match their strategy name
  mwm_results <- read_xlsx(strat_export_fp)  # Read the MWM results file to confirm
  mwm_results <- mwm_results %>%
    rename(
      thigmotaxis = `1`,
      circling = `2`,
      random_path = `3`,
      scanning = `4`,
      chaining = `5`,
      directed_search = `6`,
      corrected_search = `7`,
      direct_path = `8`,
      perseverance = `9`
    )
  
  # pd column
  if (!"pd" %in% colnames(mwm_results)) {
    mwm_results$pd <- NA 
  }
  
  # Loop through cohorts and update pd based on spatial data
  for (cohort in cohort_list) {
    cohort_path <- file.path(rtrack_folder, "Cohorts", paste0("Cohort", cohort))
    cohort_files <- list.files(cohort_path, pattern = paste0("^Coh", cohort, "_.*\\.xlsx$"), full.names = TRUE)
    
    if (length(cohort_files) == 0) {
      message(paste("No file found for Cohort", cohort))
      next
    }
    
    cohort_file <- cohort_files[1]  # Use the first matched file
    
    spatial_data <- tryCatch({
      read_excel(cohort_file, sheet = "Spatial") %>%
        mutate(pd = case_when(
          Drop %in% c(2, 8) ~ 1,
          Drop %in% c(3, 7) ~ 2,
          Drop %in% c(4, 6) ~ 3,
          Drop == 5 ~ 4,
          is.na(Drop) ~ NA_real_, 
        )) %>%
        select(Animal, Trial, pd)  
    }, error = function(e) {
      message(paste("Error reading 'Spatial' sheet for Cohort", cohort, ":", e$message))
      return(NULL)
    })
    
    if (is.null(spatial_data)) {
      next
    }
    
    # Convert Animal column to character to match _TargetID
    spatial_data$Animal <- as.character(spatial_data$Animal)
    mwm_results$`_TargetID` <- as.character(mwm_results$`_TargetID`)
    
    # Merge spatial data with mwm_results based on Animal and Trial columns
    mwm_results <- mwm_results %>%
      left_join(spatial_data, by = c("_TargetID" = "Animal", "_Trial" = "Trial")) %>%
      mutate(pd = coalesce(pd.y, pd.x)) %>%  
      select(-pd.x, -pd.y)  
  }
  
  write.xlsx(mwm_results, strat_export_fp, sheetName = "Spatial", overwrite = TRUE)
  
  message("\n")
  message(paste("File has been created successfully:", strat_export_fp))
}

strat_export()