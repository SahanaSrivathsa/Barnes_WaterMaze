from pynwb import NWBHDF5IO
from pynwb.file import Subject

nwb_path = "/Users/miasponseller/Desktop/dandi_upload/sub-01/sub-01_ses-01_behavior.nwb"

with NWBHDF5IO(nwb_path, mode="r+") as io:
    nwbfile = io.read()

    if nwbfile.subject is None:
        nwbfile.subject = Subject(
            subject_id="sub-01",
            species="Rattus norvegicus",
            sex="M",
            age="P6M/P22M", # between 6 months and 22 months
        )
        io.write(nwbfile)
        print("Subject added.")
    else:
        print("Subject already exists.")
