from datetime import datetime
from zoneinfo import ZoneInfo
from pathlib import Path
from neuroconv.datainterfaces import ExcelTimeIntervalsInterface
from pynwb import NWBHDF5IO
from pynwb.file import Subject

# Input and output files
excel_path = Path("/Users/miasponseller/Desktop/MWM_results_04-06-2025_for_dandi.xlsx")
nwb_path = Path("/Users/miasponseller/Desktop/dandi_upload/sub-01/sub-01_ses-01_behavior.nwb")

# Convert Excel to NWB with metadata
interface = ExcelTimeIntervalsInterface(file_path=excel_path, verbose=False)
metadata = interface.get_metadata()
metadata["NWBFile"] = dict(
    session_start_time=datetime(2025, 4, 6, 0, 0, 0, tzinfo=ZoneInfo("America/Phoenix"))
)

interface.run_conversion(nwbfile_path=nwb_path, metadata=metadata)
print(f"NWB file created at: {nwb_path}")

# Add Subject metadata
with NWBHDF5IO(nwb_path, mode="r+") as io:
    nwbfile = io.read(write=True)

    if nwbfile.subject is None:
        nwbfile.subject = Subject(
            subject_id="sub-01",
            species="Rattus norvegicus", # has to be latin binomial
            sex="M",
            age="P6M/P22M",  # between 6 months and 22 months
        )
        io.write(nwbfile)
        print("Subject metadata added.")
    else:
        print("Subject already exists.")
