library(Rtrack)
library(dplyr)
library(openxlsx)
library(readxl)
library(ggplot2)
library(data.table)

# Rtrack folder file path (main directory)
rtrack_folder <- "/Users/miasponseller/Desktop/Lab/Rtrack"

# Rat list file path (contains cohort, rat number, and age group)
rat_list <- read.csv("/Users/miasponseller/Desktop/Lab/Rtrack/Rat_List.csv")

# Cohort list and their corresponding file paths
cohort_list <- c(1, 2, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21) # Cohort numbers
cohort_file_paths <- sapply(cohort_list, function(cohort) {
  file.path(rtrack_folder, "/Cohorts/", paste0("Cohort", cohort))},
  USE.NAMES = FALSE)

# all rats experiment description file (if exists, must be an .xlsx)
all_rats_desc_fp <- file.path(rtrack_folder, "All_Rats.xlsx")


# Add to or create AllRats.xlsx (all cohorts combined description file) ------------------------

all_rats_experiment_desc_file <- function(cohort_list, cohort_file_paths, rat_list, rtrack_folder) {
  all_rats_fp <- file.path(rtrack_folder, "All_Rats.xlsx") # output file path
  all_rats <- data.frame()
  
  # Ensure cohort list and cohort file paths are the same length
  if (length(cohort_list) != length(cohort_file_paths)) {
    stop("The cohort list and cohort file paths list are not the same length.")
  }
  
  for (i in seq_along(cohort_list)) {
    cohort <- cohort_list[i]
    cohort_fp <- cohort_file_paths[i]
    exp_desc_file <- file.path(cohort_fp, paste0("cohort", cohort, "_exp_desc.xlsx"))
    
    if (file.exists(exp_desc_file)) {
      cohort_data <- tryCatch({
        read_excel(exp_desc_file)
      }, error = function(e) {
        warning(paste("Error reading", exp_desc_file, ":", e$message))
        return(NULL)
      })
      
      if (!is.null(cohort_data)) {
        required_cols <- c("_TrackID", "_TargetID", "_Day", "_Trial", "_Arena", "_TrackFile", "_TrackFileFormat")
        missing_cols <- setdiff(required_cols, colnames(cohort_data))
        if (length(missing_cols) > 0) {
          warning(paste("Missing columns in", exp_desc_file, ":", paste(missing_cols, collapse = ", ")))
          next
        }
        
        rat_list_unique <- rat_list %>% distinct(Cohort, Number, Age)
        
        # Add Cohort and Age from rat_list
        cohort_data <- cohort_data %>% 
          left_join(rat_list_unique %>% filter(Cohort == cohort), by = c("_TargetID" = "Number"), relationship = "many-to-many") %>% 
          mutate(Cohort = cohort)
        
        all_rats <- bind_rows(all_rats, cohort_data)
      }
    } else {
      warning(paste("Experiment description file not found for cohort", cohort))
    }
  }
  
  # Combine with existing data if output file already exists
  if (file.exists(all_rats_fp)) {
    existing_data <- read.csv(all_rats_fp, stringsAsFactors = FALSE)
    all_rats <- anti_join(all_rats, existing_data, by = names(all_rats))
    all_rats <- bind_rows(existing_data, all_rats)
  }
  
  # Write combined data to an Excel file
  if (nrow(all_rats) > 0) {
    dir.create(rtrack_folder, showWarnings = FALSE, recursive = TRUE)
    openxlsx::write.xlsx(all_rats, all_rats_fp, overwrite = TRUE)
    message(paste("File", basename(all_rats_fp), "has been created/updated successfully"))
  } else {
    message("No new data to append.")
  }
}


# Single Trial Processing and Plots ---------------------------------------

process_single_path <- function(cohort_num, trial_num, output_folder = "/Users/miasponseller/Desktop/Lab/Rtrack/Plot_PDFs") {
  if (!(cohort_num %in% cohort_list)) {
    stop("Cohort number not found.")
  }
  
  cohort_index <- which(cohort_list == cohort_num)
  cohort_path <- cohort_file_paths[cohort_index]
  track_folder <- file.path(cohort_path, paste0("Cohort", cohort_num, " Trials"))
  print(track_folder)
  trial_file <- file.path(track_folder, paste0("Coh", cohort_num, "_Trial", trial_num, ".csv"))
  print(trial_file)
  
  if (!file.exists(trial_file)) {
    stop(paste("Trial", trial_num, "file not found for Cohort", cohort_num))
  }
  
  # Read in data and calculate metrics
  arena = Rtrack::read_arena(file.path(cohort_path, paste0("Cohort", cohort_num, "Arena.txt")))
  path = Rtrack::read_path(trial_file, arena, id = paste0("trial", trial_num), track.format = "anymaze.csv")
  metrics = Rtrack::calculate_metrics(path, arena)
  
  # Create output folder if it doesn't exist
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
    message("Created output directory: ", output_folder)
  }
  
  plot_path(metrics)
  #plot_density(metrics)
  
  # # Save the plots to a PDF file
  # pdf_file <- file.path(output_folder, paste0("cohort", cohort_num, "_trial", trial_num, "_plots.pdf")) # output file name
  # 
  # tryCatch({
  #   pdf(file = pdf_file)
  #   message("Saving path and density plot for cohort ", cohort_num, ", trial ", trial_num, " to ", pdf_file)
  #   plot_path(metrics) # path trajectory
  #   plot_density(metrics) # density heat map
  #   dev.off()
  #   message("Finished saving plots for cohort ", cohort_num, ", trial ", trial_num)
  # }, error = function(e) {
  #   warning(paste("Error generating PDF for cohort", cohort_num, "trial", trial_num, ": ", e$message))
  # })
  
}


# Bulk Experiment Processing ----------------------------------------------

bulk_process_experiment <- function() {
  for (cohort in cohort_list) {
    desc_file <- file.path(cohort_file_paths[which(cohort_list == cohort)], paste0("cohort", cohort, "_exp_desc.xlsx"))
    experiment <- Rtrack::read_experiment(desc_file, data.dir = file.path(cohort_file_paths[which(cohort_list == cohort)], "Trials"))
  }
}



# Plot All Paths and Save -------------------------------------------------
plot_all_paths <- function(cohort_list, output_folder = "/Users/miasponseller/Desktop/Lab/Rtrack/Plot_PDFs") {
  # Create output folder if it doesn't exist
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
    message("Created directory: ", output_folder)
  }
  
  # Iterate through each cohort in the list
  for (cohort in cohort_list) {
    # Define the path to the experiment description file
    desc_file <- file.path(cohort_file_paths[which(cohort_list == cohort)], paste0("cohort", cohort, "_exp_desc.xlsx"))
    
    # Attempt to read the experiment description file and handle errors
    experiment <- tryCatch({
      Rtrack::read_experiment(desc_file, data.dir = file.path(cohort_file_paths[which(cohort_list == cohort)], paste0("Cohort", cohort, " Trials")))
    }, error = function(e) {
      warning(paste("Error reading experiment for cohort", cohort, ":", e$message))
      return(NULL)
    })
    
    # Skip processing if experiment data couldn't be read
    if (is.null(experiment)) {
      next
    }
    
    # Extract strategies for the experiment
    strategies <- tryCatch({
      call_strategy(experiment)
    }, error = function(e) {
      warning(paste("Error in call_strategy for cohort", cohort, ":", e$message))
      return(NULL)
    })
    
    # Skip if strategies are not available
    if (is.null(strategies)) {
      next
    }
    
    # Define output PDF file for the cohort
    pdf_file <- file.path(output_folder, paste0("cohort", cohort, "_all_paths", ".pdf"))
    
    # Create PDF and plot all paths for the cohort
    tryCatch({
      pdf(file = pdf_file)
      message("Saving plots for cohort ", cohort, " to ", pdf_file)
      
      for (i in seq_along(experiment$metrics)) {
        tryCatch({
          plot_path(experiment$metrics[[i]], title = paste0(experiment$metrics[[i]]$id, " - ", strategies$calls[i, "name"]))
        }, error = function(e) {
          warning(paste("Error plotting path for cohort", cohort, "metric index", i, ":", e$message))
        })
      }
      
      dev.off()
      message("Finished saving plots for cohort ", cohort)
    }, error = function(e) {
      warning(paste("Error generating PDF for cohort", cohort, ": ", e$message))
    })
  }
  message("All plots saved successfully.")
  
}



# Bulk Processing a Whole Experiment --------------------------------------

bulk_strategy_calling <- function() {
  experiment_file <- all_rats_desc_fp
  trials_dir <- "/Users/miasponseller/Desktop/Lab/Rtrack/All Trials"
  
  experiment <<- Rtrack::read_experiment(experiment_file, data.dir = trials_dir)
  strategies <<- Rtrack::call_strategy(experiment)
  list(experiment = experiment, strategies = strategies)
}


# Strategy plots ----------------------------------------------------------

strategy_plots <- function() {
  # if bulk_strategy_calling() has not alrady been run, run it
  if (!exists("experiment", envir = .GlobalEnv) || !exists("strategies", envir = .GlobalEnv)) {
    bulk_strategy_calling()
  }
  
  # Strategy plot, across all rats
  Rtrack::plot_strategies(strategies, experiment = experiment)
  
  # Strategy plots, by Cohort
  #Rtrack::plot_strategies(strategies, experiment = experiment, factor = "Cohort")
  
  # Strategy plots, by Age
  Rtrack::plot_strategies(strategies, experiment = experiment, factor = "Age")
}


# Analysis of selected metrics --------------------------------------------
selected_metrics_plot <- function() {
  # if bulk_strategy_calling() has not alrady been run, run it
  if (!exists("experiment", envir = .GlobalEnv) || !exists("strategies", envir = .GlobalEnv)) {
    bulk_strategy_calling()
  }
  
  # Available metrics
  print(experiment$summary.variables)
  
  # plot selected metrics by factor (first item listed)
  Rtrack::plot_variable("velocity", experiment = experiment, factor = "Age",
                        factor.colors = c(young = "green", old = "purple"),
                        lwd = 1.5)
  
  message("Selected metrics have been plotted successfully.")
}


# Bulk density maps -------------------------------------------------------
bulk_density_map <- function() {
  # if bulk_strategy_calling() has not alrady been run, run it
  if (!exists("experiment", envir = .GlobalEnv) || !exists("strategies", envir = .GlobalEnv)) {
    bulk_strategy_calling()
  }
  
  old_metrics = experiment$metrics[experiment$factors$Age == "old"]
  young_metrics = experiment$metrics[experiment$factors$Age == "young"]
  
  par(mfrow = c(1,2))
  Rtrack::plot_density(old_metrics, title = "Old")
  Rtrack::plot_density(young_metrics, title = "Young")
  
  message("Density maps have been plotted successfully.")
}

# Run Functions ------------------------------------------------------------

#all_rats_experiment_desc_file(cohort_list, cohort_file_paths, rat_list, rtrack_folder)
process_single_path(13, 5)
#bulk_process_experiment()
#plot_all_paths(cohort_list)
#bulk_strategy_calling()
#strategy_plots()
#selected_metrics_plot()
#bulk_density_map()
#strat_export()



