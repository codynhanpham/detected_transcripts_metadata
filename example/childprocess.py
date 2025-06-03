import json, subprocess

json_string = subprocess.run(
    [
        "../target/release/detected_transcripts_metadata.exe", # or path to the executable
        "-i",
        "/path/to/detected_transcripts.csv",
        "-o",
        "-",
        "-c",
        "512", # in KiB
        "-q",
    ],
    shell=False, capture_output=True
)

data = json.loads(json_string.stdout)
print(data)