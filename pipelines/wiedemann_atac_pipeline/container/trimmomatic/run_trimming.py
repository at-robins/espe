#!/usr/bin/python
"""This module runs the trimming process."""

import json
import math
import multiprocessing
import os
import sys

BASE_COMMAND = "java -jar /Trimmomatic-0.39/trimmomatic-0.39.jar PE"
STEP_OPTIONS = ("ILLUMINACLIP:/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10:2:True "
"LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36")
MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = MOUNT_PATHS["input"] + "/"

# If a specific environment variable is set, appends the respective option.
options = ""

threads = math.floor(multiprocessing.cpu_count() * 0.8)
if threads > 0:
    options += f" -threads {threads}"

phred = os.environ.get("PHRED")
if phred is not None:
    if phred == "PHRED33":
        options += " -phred33"
    elif phred == "PHRED64":
        options += " -phred64"
    else:
        print(f"Unknown PHRED score option: {phred}", file=sys.stderr)

if not options:
    print("Running with default options.")
else:
    print("Specified options:" + options)

# Iterates over all sample directories and processes them conserving the directory structure.
for root, dirs, files in os.walk(INPUT_FOLDER):
    if len(files) > 0:
        for file in files:
            input_files = ""
            file_base_name = ""
            file_base_input_path = ""
            if file.casefold().endswith("_1.fq.gz"):
                file_base_name = file.removesuffix("_1.fq.gz")
                file_base_input_path = os.path.join(root, file_base_name)
                input_files = f"{file_base_input_path}_1.fq.gz {file_base_input_path}_2.fq.gz"
            elif file.casefold().endswith("_1.fastq.gz"):
                file_base_name = file.removesuffix("_1.fastq.gz")
                file_base_input_path = os.path.join(root, file_base_name)
                input_files = f"{file_base_input_path}_1.fastq.gz {file_base_input_path}_2.fastq.gz"

            if input_files:
                file_base_output_path = os.path.join(
                    MOUNT_PATHS["output"],
                    file_base_input_path.removeprefix(INPUT_FOLDER)
                )
                full_command = (f"{BASE_COMMAND}{options} {input_files} "
                f"{file_base_output_path}_1_paired.fq.gz "
                f"{file_base_output_path}_1_unpaired.fq.gz "
                f"{file_base_output_path}_2_paired.fq.gz "
                f"{file_base_output_path}_2_unpaired.fq.gz "
                f"{STEP_OPTIONS}")
                os.makedirs(os.path.dirname(file_base_output_path), exist_ok = True)
                exit_code = os.waitstatus_to_exitcode(os.system(full_command))
                if exit_code != 0:
                    sys.exit(exit_code)
