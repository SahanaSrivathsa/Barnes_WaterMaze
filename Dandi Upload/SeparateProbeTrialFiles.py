import pandas as pd
import numpy as np

def split_on_platform(csv_file, x_plat, y_plat, radius,
                      before_file='/Users/miasponseller/Desktop/Lab/Rtrack/WMaze/Cohorts/Cohort1/before_platform.csv',
                      after_file='/Users/miasponseller/Desktop/Lab/Rtrack/WMaze/Cohorts/Cohort1/after_platform.csv'):

    # load trial path data
    trial_path = pd.read_csv(csv_file)

    # Rename columns (needed for raw.csv format in Rtrack)
    trial_path = trial_path.rename(columns={'Centre position X': 'X', 'Centre position Y': 'Y'})

    # platform area distance
    trial_path['distance'] = np.sqrt((trial_path['X'] - x_plat)**2 + (trial_path['Y'] - y_plat)**2)

    # find first point where the distance is <= radius
    inside = trial_path['distance'] <= radius
    if not inside.any():
        print("No points inside the platform area.")
        trial_path.drop(columns=['distance']).to_csv(before_file, index=False)
        return

    entry_idx = inside[inside].index[0]                # first index inside platform
    entry_pos = trial_path.index.get_loc(entry_idx)    # numeric location

    # Split trial path
    before_df = trial_path.iloc[:entry_pos].drop(columns=['distance'])
    after_df = trial_path.iloc[entry_pos:].drop(columns=['distance'])

    # Save to CSVs
    before_df.to_csv(before_file, index=False)
    after_df.to_csv(after_file, index=False)
    print(f"Data split and saved to {before_file} and {after_file}.")

split_on_platform('/Users/miasponseller/Desktop/Lab/Rtrack/WMaze/Cohorts/Cohort1/Probe Trials/Coh1_Trial290.csv', x_plat = 415, y_plat = 169.5, radius = 23.3)