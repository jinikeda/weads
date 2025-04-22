import os

def fileexists(file):
    if not isinstance(file, str) or not file.strip():
        raise ValueError(f"Invalid file path: {file!r}")

    if not os.path.exists(file):
        raise FileNotFoundError(f"File not found: {file}")
# ----------------------------------------------------------
