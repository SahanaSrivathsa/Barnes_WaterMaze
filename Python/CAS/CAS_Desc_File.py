import pandas as pd

# Load the Excel file
input_file = '/Users/miasponseller/Desktop/Lab/Rtrack/CAS/CAS_AllRats_Spatial.csv'
df = pd.read_csv(input_file)

# Create new columns
df['_TrackID'] = 'Coh' + df['Cohort'].astype(str) + '_test' + df['Test'].astype(str)
df['_TargetID'] = df['Animal']
df['_Trial'] = df['Trial']
df['_Day'] = ((df['Trial'] - 1) // 6 + 1).clip(upper=4)
df['_TrackFileFormat'] = 'anymaze.csv'
df['_Arena'] = 'Arena Files/Cohort' + df['Cohort'].astype(str) + 'Arena.txt'
df['_TrackFile'] = 'Coh' + df['Cohort'].astype(str) + '_Trial' + df['Test'].astype(str) + '.csv'

# Define expected trial labels
expected_trials = set(range(1, 25))

# Identify valid animals (have exactly trials 1-24)
trial_sets = df.groupby('_TargetID')['_Trial'].apply(lambda x: set(x.dropna()))
valid_animals = trial_sets[trial_sets == expected_trials].index
invalid_animals = trial_sets[trial_sets != expected_trials]

# Filter dataframe to keep only valid animals
filtered_df = df[df['_TargetID'].isin(valid_animals)]

# Save filtered data
output_columns = [
    '_TrackID', '_TargetID', '_Trial', '_Day',
    '_TrackFileFormat', 'Age', 'Cohort', '_Arena', '_TrackFile', 'Performance'
]
output_df = filtered_df[output_columns]

output_file = '/Users/miasponseller/Desktop/Lab/Rtrack/CAS/CAS_exp_desc.xlsx'
output_df.to_excel(output_file, index=False)

print(f"Experiment description file saved to {output_file}")

# Print animals with missing or misnumbered trials
print("\nAnimals excluded (missing trials 1–24 or trial labels incorrect):")
for rat, trials in invalid_animals.items():
    print(f"Rat {rat} Trials: {sorted(trials)}")