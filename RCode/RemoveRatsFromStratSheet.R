library(readxl)
library(writexl)
library(dplyr)

# Strategy sheet from Rtrack
strat_sheet <- '/Users/miasponseller/Desktop/Lab/Rtrack/Tg/Tg_MWM_results_06-11-2025.xlsx'

# Output file path
output_file <- '/Users/miasponseller/Desktop/Lab/Rtrack/Tg/rats_removed_Tg_MWM_results_06-11-2025.xlsx'

# Animal numbers to remove
rat_nums <- c(11055, 11066, 11067, 11069, 11070, 11072, 11074, 11077, 11078, 11082, 
              11084, 10927, 11019, 11029, 11031, 11035, 11022, 10997, 10998, 11001,
              11003, 11004, 11005, 11006, 11007, 11015, 11017, 075, 020, 079, 084, 
              021, 071, 078, 050, 102)

df <- read_excel(strat_sheet)

cleaned_df <- df %>% 
  filter(!(`_TargetID` %in% rat_nums))

write_xlsx(cleaned_df, output_file)

message('Updated file saved to: ', output_file)