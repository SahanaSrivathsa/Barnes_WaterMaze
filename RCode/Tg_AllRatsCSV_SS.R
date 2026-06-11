library(readxl)
library(dplyr)
library(stringr)
library(purrr)
library(writexl)
library(tidyr)
### SS MODIFICATION ####
# I changed the files to get barnes id, sex ,treatment type/genotype from the sheet
# Tg_AllRats_Spatial.xlsx in the base folder. ie master copy

###### MIA ORIGINAL #####
# Folder with all ANY-Maze spreadsheets, which should include Key and Spatial sheets.
# Key columns: BarnesID or CowenID, Age
# not necessary, but can include for Key: Sex, treatment type
# Spatial columns: Test, Animal, Trial
# not necessary, but can include for Spatial: Duration, Distance, Mean speed, Path efficiency, CIPL
############################################################################################
all_rat_file <- "D:/NARP_Data/RTrack_NARPMale/TgF344_Aging_AllRats.xlsx"
spatial_sheets_folder <- "D:/NARP_Data/RTrack_NARPMale/Spatial_Sheets"
output_file <- "D:/NARP_Data/RTrack_NARPMale/Tg_AllRats_Spatial.csv"

# Read all rat file
all_rat_df <- read_excel(all_rat_file, sheet = "ALL", skip = 1)


# Clean column names inline
names(all_rat_df) <- names(all_rat_df) %>%
  str_replace_all("\\s+", "") %>%
  str_replace_all("_", "") %>%
  str_replace_all("/", "")

# Rename cleaned columns: Because R cant read thsese
# Barnes ID -> BarnesID
# Cowen ID-> CowenID
# JHU ID -> JHUID
# APP/WT -> APPWT
# Age at Behavior Start -> AgeatBehaviorStart

id_cols <- intersect(c("BarnesID", "CowenID", "JHUID"), colnames(all_rat_df))

all_rat_df <- all_rat_df %>%
  select(any_of(c(
    id_cols,
    "Sex",
    "APPWT",
    "AgeatBehaviorStart"
  ))) %>%
  rename(
    Age = AgeatBehaviorStart,
    Genotype = APPWT
  )

all_rat_long <- all_rat_df %>%
  pivot_longer(
    cols = all_of(id_cols),
    names_to = "IDtype",
    values_to = "Animal"
  ) %>%
  mutate(
    Animal = as.character(Animal),
    Animal = str_trim(Animal)
  ) %>%
  filter(!is.na(Animal), Animal != "") %>%
  distinct(Animal, .keep_all = TRUE)

# List all .xlsx files in the folder
spatial_files <- list.files(spatial_sheets_folder, pattern = "\\.xlsx$", full.names = TRUE)
print(spatial_files)

# Empty list to store dataframes
all_data <- list()

for (file_path in spatial_files) {
  file_name <- basename(file_path)
  
  # Get cohort number from file name
  cohort_num <- str_match(file_name, "Coh([^-]+)-")[, 2]
  if (is.na(cohort_num)) {
    message(paste("No cohort found for", file_name))
    cohort_num <- "Unknown"
  }
  
  tryCatch({
    
    # Read only Spatial sheet from cohort file
    spatial_df <- read_excel(file_path, sheet = "Spatial")
    
    # Clean Spatial column names inline
    names(spatial_df) <- names(spatial_df) %>%
      str_replace_all("\\s+", "") %>%
      str_replace_all("_", "") %>%
      str_replace_all("/", "")
    # Find CIPL column even if named "Platform : CIPL", "Platform : CIPL (m·s)", etc.
    cipl_col <- names(spatial_df)[str_detect(names(spatial_df), regex("CIPL", ignore_case = TRUE))]
    
    if (length(cipl_col) == 0) {
      stop(paste("No CIPL column found in", file_name))
    }
    
    # Use the first CIPL column if more than one exists
    spatial_df <- spatial_df %>%
      rename(CIPL = all_of(cipl_col[1]))
    spatial_df <- spatial_df %>%
      select(any_of(c(
        "Test",
        "Animal",
        "Trial",
        "Pathefficiency",
        "CIPL"
      ))) %>%
      mutate(
        Animal = as.character(Animal),
        Animal = str_trim(Animal)
      )
    
    # Merge Spatial sheet with master rat metadata
    merged_df <- spatial_df %>%
      left_join(all_rat_long, by = "Animal") %>%
      mutate(Cohort = cohort_num) %>%
      select(
        Test,
        Animal,
        Trial,
        Pathefficiency,
        Cohort,
        Sex,
        Genotype,
        Age,
        CIPL
      )
    
    all_data[[length(all_data) + 1]] <- merged_df
    
  }, error = function(e) {
    message(paste("Error reading", file_name, ":", e$message))
  })
}

# Combine all data and write to file
if (length(all_data) == 0) {
  stop("No files were successfully read. Check Spatial sheet names and Animal IDs.")
}

final_df <- bind_rows(all_data)

write.csv(final_df, output_file, row.names = FALSE)

cat("All data combined and saved to", output_file, "\n")

