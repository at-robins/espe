#!/usr/bin/python
"""This module removes blacklisted regions."""

import json
import os
import subprocess

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = MOUNT_PATHS["dependencies"]["seqoutbias"]
BLACKLIST_PATH = os.path.join(MOUNT_PATHS["globals"]["GENOME"], "blacklist.bed")
BG_SUFFIX = ".bedGraph"

# Iterates over all sample directories and processes them conserving the directory structure.
for root, dirs, files in os.walk(INPUT_FOLDER):
    for file in files:
        if file.endswith(BG_SUFFIX):
            file_base_name = file.removesuffix(BG_SUFFIX)
            file_base_path = os.path.join(root, file_base_name)
            file_base_output_path = os.path.join(
                MOUNT_PATHS["output"],
                file_base_path.removeprefix(INPUT_FOLDER + "/"),
            )
            full_command = (
                f"bedtools intersect -a {file_base_path}{BG_SUFFIX} -b {BLACKLIST_PATH} -wa -v "
                f"> {file_base_output_path}{BG_SUFFIX}"
            )
            print(f"Running command: {full_command}")
            os.makedirs(os.path.dirname(file_base_output_path), exist_ok=True)
            subprocess.run(
                full_command,
                shell=True,
                check=True,
            )
