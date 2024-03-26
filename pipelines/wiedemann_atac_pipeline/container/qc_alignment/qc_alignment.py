#!/usr/bin/python
"""This module runs the post alignment QC."""

import json
import math
import multiprocessing
import os
import sys

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = next(iter(MOUNT_PATHS["dependencies"].values())) + "/"

# If a specific environment variable is set, appends the respective option.
options = " -outformat html -c"

threads = math.floor(multiprocessing.cpu_count() * 0.8)
if threads > 0:
    options += f" -nt {threads}"

phred = os.environ.get("GLOBAL_ORGANISM")
if phred == "mouse":
    options += " -gd MOUSE"
elif phred == "human":
    options += " -gd HUMAN"
else:
    print(f"Unknown organsim option: {phred}", file=sys.stderr)
    sys.exit(1)

print("Specified options:" + options)

if not options:
    print("Running with default step options.")
else:
    print("Specified step options:" + options)

# Iterates over all sample directories and processes them conserving the directory structure.
BAM_SUFFIX = ".bam"
for root, dirs, files in os.walk(INPUT_FOLDER):
    for file in files:
        if file.endswith(BAM_SUFFIX):
            file_path = os.path.join(root, file)
            output_path = os.path.join(
                MOUNT_PATHS["output"],
                root.removeprefix(INPUT_FOLDER),
            )
            full_command = (
                f"qualimap bamqc -bam {file_path} "
                f"-outdir {output_path} "
                f"{options} --java-mem-size=4G"
            )
            print(f"Running command: {full_command}", flush=True)
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            exit_code = os.waitstatus_to_exitcode(os.system(full_command))
            if exit_code != 0:
                sys.exit(exit_code)
