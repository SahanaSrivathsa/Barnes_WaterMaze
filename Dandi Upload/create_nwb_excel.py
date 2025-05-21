from datetime import datetime
from zoneinfo import ZoneInfo
from pathlib import Path
from neuroconv.datainterfaces import ExcelTimeIntervalsInterface

file_path = "/Users/miasponseller/Desktop/MWM_results_04-06-2025_for_dandi.xlsx"
interface = ExcelTimeIntervalsInterface(file_path=file_path, verbose=False)

metadata = interface.get_metadata()
session_start_time = datetime(2025, 4, 6, 0, 0, 0, tzinfo=ZoneInfo("America/Phoenix"))
metadata["NWBFile"] = dict(session_start_time=session_start_time)

nwbfile_path = "/Users/miasponseller/Desktop/MWM_results_04-06-2025_for_dandi.nwb"
nwbfile = interface.run_conversion(nwbfile_path=nwbfile_path, metadata=metadata)