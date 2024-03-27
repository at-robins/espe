#!/usr/bin/python
"""This module merges peak regions."""

import json
import os
import shutil
import subprocess

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = next(iter(MOUNT_PATHS["dependencies"].values()))
NARROW_SUFFIX = ".narrowPeak"
BROAD_SUFFIX = ".broadPeak"

narrow_combined_path = os.path.join(MOUNT_PATHS["output"], f"combined{NARROW_SUFFIX}")
broad_combined_path = os.path.join(MOUNT_PATHS["output"], f"combined{BROAD_SUFFIX}")
narrow_sorted_path = os.path.join(
    MOUNT_PATHS["output"], f"combined_sorted{NARROW_SUFFIX}"
)
broad_sorted_path = os.path.join(
    MOUNT_PATHS["output"], f"combined_sorted{BROAD_SUFFIX}"
)
narrow_merged_path = os.path.join(MOUNT_PATHS["output"], f"merged{NARROW_SUFFIX}")
broad_merged_path = os.path.join(MOUNT_PATHS["output"], f"merged{BROAD_SUFFIX}")


def copy_into(in_file, out_file):
    """
    Copies the contents of the input file into the output file.
    """
    with open(out_file, mode="at", encoding="utf-8") as o_f:
        with open(in_file, mode="rt", encoding="utf-8") as i_f:
            shutil.copyfileobj(i_f, o_f)


# Combines the peak files.
for root, dirs, files in os.walk(INPUT_FOLDER):
    for file in files:
        input_file = os.path.join(root, file)
        if file.endswith(NARROW_SUFFIX):
            print(f"Processing file {os.path.join(root, file)}...")
            copy_into(in_file=input_file, out_file=narrow_combined_path)
        elif file.endswith(BROAD_SUFFIX):
            print(f"Processing file {os.path.join(root, file)}...")
            copy_into(in_file=input_file, out_file=broad_combined_path)

print("Sorting files...", flush=True)
subprocess.run(
    (
        f"sort -k1,1V -k2,2n {narrow_combined_path} "
        '| awk \'{print $0"\\tPeak"sprintf("%09d", NR);}\' '
        f"> {narrow_sorted_path}"
    ),
    cwd=MOUNT_PATHS["output"],
    shell=True,
    check=True,
)
subprocess.run(
    (
        f"sort -k1,1V -k2,2n {broad_combined_path} "
        '| awk \'{print $0"\\tPeak"sprintf("%09d", NR);}\' '
        f"> {broad_sorted_path}"
    ),
    cwd=MOUNT_PATHS["output"],
    shell=True,
    check=True,
)
print("Merging files...", flush=True)
subprocess.run(
    f"bedtools merge -i {narrow_sorted_path} > {narrow_merged_path}",
    cwd=MOUNT_PATHS["output"],
    shell=True,
    check=True,
)
subprocess.run(
    f"bedtools merge -i {broad_sorted_path} > {broad_merged_path}",
    cwd=MOUNT_PATHS["output"],
    shell=True,
    check=True,
)
