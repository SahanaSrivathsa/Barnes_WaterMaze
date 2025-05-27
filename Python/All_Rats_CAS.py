import pandas as pd

# Load the Excel file
input_file = '/Users/miasponseller/Desktop/Lab/Rtrack/CAS/CAS_AllRats_Spatial.csv'
df = pd.read_excel(input_file)

# Create new columns
df['_TrackID'] = 'Coh' + df['Cohort'].astype(str) + '_test' + df['Test'].astype(str)
df['_TargetID'] = df['Animal']
df['_Trial'] = df['Trial']

# Map trial number to day (1–6 -> 1, 7–12 -> 2, etc.)
df['_Day'] = ((df['Trial'] - 1) // 6 + 1).clip(upper=4)

df['_TrackFileFormat'] = 'anymaze.csv'
df['Age'] = df['Age']
df['Performance'] = df['Performance']
df['Cohort'] = df['Cohort']
df['_Arena'] = 'Cohort' + df['Cohort'].astype(str) + 'Arena.txt'
df['_TrackFile'] = 'Coh' + df['Cohort'].astype(str) + '_Trial' + df['Test'].astype(str)

# Select and order columns as specified
output_columns = [
    '_TrackID', '_TargetID', '_Trial', '_Day',
    '_TrackFileFormat', 'Age', 'Cohort', '_Arena', '_TrackFile'
]
output_df = df[output_columns]

# Save to new Excel file
output_file = '/Users/miasponseller/Desktop/Lab/Rtrack/CAS/CAS_exp_desc.xlsx'
output_df.to_excel(output_file, index=False)

print(f"Saved transformed data to {output_file}")
