import pandas as pd

# Load the Excel file
input_file = '/Users/miasponseller/Desktop/Lab/Rtrack/CAS/CAS_AllRats_Spatial.csv'
df = pd.read_csv(input_file)

# Filter animals with exactly 24 trials
trial_counts = df.groupby('Animal')['Trial'].count()
animals_with_24_trials = trial_counts[trial_counts == 24].index

df = df[df['Animal'].isin(animals_with_24_trials)]
removed_rats = trial_counts[trial_counts != 24].index.tolist()

# Create new columns
df['_TrackID'] = 'Coh' + df['Cohort'].astype(str) + '_test' + df['Test'].astype(str)
df['_TargetID'] = df['Animal']
df['_Trial'] = df['Trial']

# Map trial number to day (1–6 -> 1, 7–12 -> 2, etc.)
df['_Day'] = ((df['Trial'] - 1) // 6 + 1).clip(upper=4)

df['_TrackFileFormat'] = 'anymaze.csv'
df['_Arena'] = 'Arena Files/Cohort' + df['Cohort'].astype(str) + 'Arena.txt'
df['_TrackFile'] = 'Coh' + df['Cohort'].astype(str) + '_Trial' + df['Test'].astype(str) + '.csv'

# Select and order columns as specified
output_columns = [
    '_TrackID', '_TargetID', '_Trial', '_Day',
    '_TrackFileFormat', 'Age', 'Cohort', '_Arena', '_TrackFile', 'Performance'
]
output_df = df[output_columns]

# Save to new Excel file
output_file = '/Users/miasponseller/Desktop/Lab/Rtrack/CAS/CAS_exp_desc.xlsx'
output_df.to_excel(output_file, index=False)

print(f"Saved transformed data to {output_file}")

# Print removed rats and number of trials each had
print("Removed rats (not having exactly 24 trials):")
for rat in removed_rats:
    count = trial_counts[rat]
    print(f"Rat {rat} had {count} trials.")
