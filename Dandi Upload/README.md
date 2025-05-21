This converts an Excel file into a NWB file for upload to DANDI.  For other format conversions, see https://neuroconv.readthedocs.io/en/stable/conversion_examples_gallery/index.html

1. Set the file path for the excel file that you are wanting to convert to a .nwb as the excel_path.  The output file path should be structured as sub-#_ses-#.nwb
2. You need to set a session_start_time.  I just set this as the date that the file was generated.
3. For the Subject metadata, you must include subject_id, species (has to be latin binomial such as "Rattus norvegicus"), sex, and age (date_of_birth can be used instead)
4. Once code is ran, you should have a NWB file created that can be found at the nwb_path
5. You can validate that the NWB file will be accepted by DANDI by typing "dandi validate ." into the terminal.

To upload to DANDI using NWB Guide:
1. Install NWB Guide (https://nwb-guide.readthedocs.io/en/stable/installation.html)
2. Navigate to Upload tab
3. Select the NWB file(s) you want to upload
4. If you have an existing Dandiset, you can look it up by ID number.  Otherwise, create a new Dandiset.
5. Click Upload Files.  Once complete, you should be able to the file(s) in your Dandiset under Files.
