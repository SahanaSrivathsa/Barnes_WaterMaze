library(openxlsx)
library(readxl)

all_rats_probe_desc <- file.path("/Users/miasponseller/Desktop/Lab/Rtrack/All_Rats_Probe.xlsx") # all rats experiment description file, must be an .xlsx
trials_dir <- "/Users/miasponseller/Desktop/Lab/Rtrack/All Probe Trials" # folder with all trials path data
strat_export_fp <- paste0("/Users/miasponseller/Desktop/Lab/Rtrack/MWM_results_Probe_", format(Sys.Date(), "%m-%d-%Y"), ".xlsx") # output file path


bulk_strategy_calling <- function() {
  experiment <<- Rtrack::read_experiment(all_rats_probe_desc, data.dir = trials_dir)
  strategies <<- Rtrack::call_strategy(experiment)
  list(experiment = experiment, strategies = strategies)
}


probe_strat_export <- function() {
  if (!exists("experiment", envir = .GlobalEnv) || !exists("strategies", envir = .GlobalEnv)) {
    bulk_strategy_calling()
  }
  
  Rtrack::export_results(experiment, strategies, file = strat_export_fp)
  
  # rename numbered strat columns to match their strat name
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
  
  write.xlsx(mwm_results, strat_export_fp, overwrite = TRUE)
  
  message("\n")
  message(paste("File has been created successfully:", strat_export_fp))
}

probe_strat_export()




