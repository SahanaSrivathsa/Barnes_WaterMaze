import os
import shutil

# Set the path to the main folder that contains subfolders
source_folder = '/Users/miasponseller/Desktop/Lab/Rtrack/CAS/CAS Tracks'

# Set the path to the destination folder where all files will be collected
destination_folder = '/Users/miasponseller/Desktop/Lab/Rtrack/CAS/All CAS Tracks'

# Create the destination folder if it doesn't exist
os.makedirs(destination_folder, exist_ok=True)

# Walk through all subdirectories and files
for root, dirs, files in os.walk(source_folder):
    for file in files:
        source_path = os.path.join(root, file)
        destination_path = os.path.join(destination_folder, file)

        # Avoid overwriting files with the same name by renaming if needed
        if os.path.exists(destination_path):
            base, ext = os.path.splitext(file)
            i = 1
            while os.path.exists(destination_path):
                destination_path = os.path.join(destination_folder, f"{base}_{i}{ext}")
                i += 1

        # Copy the file
        shutil.copy2(source_path, destination_path)

print("All files have been copied.")
