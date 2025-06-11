import pandas as pd
from pandas import read_excel

# Read the Excel file
fp = read_excel('/Users/miasponseller/Desktop/Lab/Rtrack/CAS/CAS_MWM_results_05-28-2025.xlsx')

# Group by _TargetID and collect the unique trial numbers
trial_sets = fp.groupby('_TargetID')['_Trial'].apply(lambda x: set(x.dropna()))

# Define the expected full set of trials
expected_trials = set(range(1, 25))

# Check which animals have complete and incomplete trial sets
complete = trial_sets[trial_sets == expected_trials]
incomplete = trial_sets[trial_sets != expected_trials]

print("\nAnimals missing one or more trials:")
for animal, trials in incomplete.items():
    present = sorted(trials)
    missing = sorted(expected_trials - trials)
    print(f"{animal}:")
    print(f"  Trials present: {present}")
    print(f"  Trials missing: {missing}")
