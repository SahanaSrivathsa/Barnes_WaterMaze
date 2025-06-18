library(readr)
library(dplyr)

# AllRatsCSV
all_rats_csv <- read_csv('/Users/miasponseller/Desktop/Lab/Rtrack/Tg/Tg_AllRats_Spatial.csv')

# Output file path
output_file <- '/Users/miasponseller/Desktop/Lab/Rtrack/Tg/Tg_AllRats_Spatial_cleaned.csv'

# Animal numbers to remove
rat_nums <- c(
  "11055", '11056', "11066","11067","11069","11070","11072","11074","11077","11078","11082",
  "11084","10927","11019","11029","11031","11035","11022","10997","10998","11001",
  "11003","11004","11005","11006","11007","11015","11017",
  "075","020","079","084","021","071","078","050","102", '11076', '11077'
)

cleaned_df <- all_rats_csv %>% 
  mutate(Animal = as.character(Animal)) %>% 
  filter(!(Animal %in% rat_nums))

write_csv(cleaned_df, output_file)

message('Updated file saved to: ', output_file)

