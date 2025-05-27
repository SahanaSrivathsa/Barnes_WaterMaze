import pandas as pd

import os
import re

# folder with CAS files and output file path
cas_folder = '/Users/miasponseller/Desktop/Lab/Rtrack/CAS/Split Spatial Sheets'
output_file = '/Users/miasponseller/Desktop/CAS_AllRats_Spatial.csv'

# performance sheets
young = pd.read_excel('/Users/miasponseller/Desktop/2023-01-18_Rats_New_age_6.xlsx')
middle = pd.read_excel('/Users/miasponseller/Desktop/2023-01-18_Rats_New_age_15.xlsx')
old = pd.read_excel('/Users/miasponseller/Desktop/2023-01-18_Rats_New_age_23.xlsx')

# ID performance dictionaries
young_perf_dict = young.set_index('ID')['performance'].to_dict()
middle_perf_dict = middle.set_index('ID')['performance'].to_dict()
old_perf_dict = old.set_index('ID')['performance'].to_dict()
all_perf_dict = {**young_perf_dict, **middle_perf_dict, **old_perf_dict}

# list all CAS files in the folder
cas_files = [f for f in os.listdir(cas_folder) if f.endswith('.xlsx')]

# empty list to store DataFrames
all_data = []

for file in cas_files:
    file_path = os.path.join(cas_folder, file)

    # cohort number from filename
    cohort_match = re.search(r'Coh([^-]+)-', file)
    cohort_number = cohort_match.group(1) if cohort_match else 'Unknown'

    try:
        key_df = pd.read_excel(file_path, sheet_name='Key')
        spatial_df = pd.read_excel(file_path, sheet_name='Spatial')
    except Exception as e:
        print(f"Error reading {file}: {e}")
        continue

    key_df = key_df[['BarnesID', 'AnyMaze#', 'Age']]

    # merge on Animal and BarnesID
    merged_df = pd.merge(spatial_df, key_df, left_on='Animal', right_on='BarnesID', how='left')
    merged_df['Cohort'] = cohort_number


    def get_performance(row):
        age = str(row['Age'])
        animal_id = row['Animal']

        if '6' in age:
            return young_perf_dict.get(animal_id, all_perf_dict.get(animal_id))
        elif '15' in age:
            return middle_perf_dict.get(animal_id, all_perf_dict.get(animal_id))
        elif '23' in age:
            return old_perf_dict.get(animal_id, all_perf_dict.get(animal_id))
        else:
            return all_perf_dict.get(animal_id)

    merged_df['Performance'] = merged_df.apply(get_performance, axis=1)
    merged_df.drop(columns=['BarnesID'], inplace=True)

    all_data.append(merged_df)

final_df = pd.concat(all_data, ignore_index=True)

final_df.to_csv(output_file, index=False)

print(f"All data combined and saved to {output_file}")