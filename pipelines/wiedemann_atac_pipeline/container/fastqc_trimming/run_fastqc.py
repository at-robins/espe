#!/usr/bin/python
"""This module runs the after trimming QC process."""

import json
import os
import sys
from contextlib import suppress

BASE_COMMAND = "perl -- /FastQC/fastqc"
MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = next(iter(MOUNT_PATHS["dependencies"].values())) + "/"

# If a specific environment variable is set, appends the respective option.
options = ""

with suppress(Exception):
    adapters = MOUNT_PATHS["globals"]["ADAPTERS"]
    options += f" --adapters {adapters}/qc_adapters.txt"

kmers = os.environ.get("KMERS")
if kmers is not None:
    options += f" --kmers {kmers}"

svg = os.environ.get("SVG")
if svg is not None and svg == "true":
    options += " --svg"

if not options:
    print("Running with default options.")
else:
    print("Specified options:" + options)

# Iterates over all sample directories and processes them conserving the directory structure.
for root, dirs, files in os.walk(INPUT_FOLDER):
    for file in files:
        if file.casefold().endswith(".fq.gz") or file.casefold().endswith(".fastq.gz"):
            file_input_path = os.path.join(root, file)
            print(f"Processing file {file_input_path}...")
            folder_output_path = os.path.join(
                MOUNT_PATHS["output"],
                root.removeprefix(INPUT_FOLDER)
            )
            full_command = (f"{BASE_COMMAND}{options} "
            f"--outdir={folder_output_path} "
            f"{file_input_path}")
            os.makedirs(folder_output_path, exist_ok = True)
            exit_code = os.waitstatus_to_exitcode(os.system(full_command))
            if exit_code != 0:
                sys.exit(exit_code)
