#!/usr/bin/python
"""This module merges peak regions."""

import csv
import json
import os
import shutil
import subprocess

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = next(iter(MOUNT_PATHS["dependencies"].values()))
NARROW_SUFFIX = ".narrowPeak"
BROAD_SUFFIX = ".broadPeak"
COMBINED_SUFFIX = ".mergedPeak"

narrow_combined_path = os.path.join(MOUNT_PATHS["output"], f"combined{NARROW_SUFFIX}")
broad_combined_path = os.path.join(MOUNT_PATHS["output"], f"combined{BROAD_SUFFIX}")
combined_path = os.path.join(MOUNT_PATHS["output"], f"combined{COMBINED_SUFFIX}")
narrow_sorted_path = os.path.join(
    MOUNT_PATHS["output"], f"combined_sorted{NARROW_SUFFIX}"
)
broad_sorted_path = os.path.join(
    MOUNT_PATHS["output"], f"combined_sorted{BROAD_SUFFIX}"
)
sorted_path = os.path.join(MOUNT_PATHS["output"], f"combined_sorted{COMBINED_SUFFIX}")
narrow_merged_path = os.path.join(MOUNT_PATHS["output"], f"merged{NARROW_SUFFIX}")
broad_merged_path = os.path.join(MOUNT_PATHS["output"], f"merged{BROAD_SUFFIX}")
merged_path = os.path.join(MOUNT_PATHS["output"], f"merged{COMBINED_SUFFIX}")


def copy_into(in_file, out_file):
    """
    Copies the contents of the input file into the output file.
    """
    with open(out_file, mode="at", encoding="utf-8") as o_f:
        with open(in_file, mode="rt", encoding="utf-8") as i_f:
            shutil.copyfileobj(i_f, o_f)


def fix_merged_file(out_file):
    """
    Adds the relevant BED columns to the merged file.
    """
    in_file = f"{out_file}_tmp"
    os.rename(out_file, in_file)
    with open(in_file, newline="", mode="rt", encoding="utf-8") as merged_in:
        with open(out_file, newline="", mode="wt", encoding="utf-8") as merged_out:
            merged_reader = csv.reader(
                merged_in, dialect="unix", delimiter="\t", quotechar='"'
            )
            merged_writer = csv.writer(
                merged_out,
                dialect="unix",
                delimiter="\t",
                quotechar='"',
                quoting=csv.QUOTE_MINIMAL,
            )
            for i, row in enumerate(merged_reader):
                merged_writer.writerow(
                    [
                        row[0],
                        row[1],
                        row[2],
                        f"peak_{i}",
                        row[3],
                        ".",  # The strand does not matter so we do not set one.
                    ]
                )
    os.remove(in_file)

def remove_narrow_peak_column(out_file):
    """
    Removes the narrow peak specific last column from the merged file.
    """
    in_file = f"{out_file}_tmp"
    os.rename(out_file, in_file)
    with open(in_file, newline="", mode="rt", encoding="utf-8") as merged_in:
        with open(out_file, newline="", mode="wt", encoding="utf-8") as merged_out:
            merged_reader = csv.reader(
                merged_in, dialect="unix", delimiter="\t", quotechar='"'
            )
            merged_writer = csv.writer(
                merged_out,
                dialect="unix",
                delimiter="\t",
                quotechar='"',
                quoting=csv.QUOTE_MINIMAL,
            )
            for i, row in enumerate(merged_reader):
                merged_writer.writerow(
                    row[0:9]
                )
    os.remove(in_file)


# Combines the peak files.
for root, dirs, files in os.walk(INPUT_FOLDER):
    for file in files:
        input_file = os.path.join(root, file)
        if file.endswith(NARROW_SUFFIX):
            print(f"Processing file {os.path.join(root, file)}...")
            copy_into(in_file=input_file, out_file=narrow_combined_path)
            copy_into(in_file=input_file, out_file=combined_path)
        elif file.endswith(BROAD_SUFFIX):
            print(f"Processing file {os.path.join(root, file)}...")
            copy_into(in_file=input_file, out_file=broad_combined_path)
            copy_into(in_file=input_file, out_file=combined_path)

print("Removing narrow peak specific columns...", flush=True)
remove_narrow_peak_column(combined_path)

print("Sorting files...", flush=True)
subprocess.run(
    (f"sort -k1,1V -k2,2n {narrow_combined_path} " f"> {narrow_sorted_path}"),
    cwd=MOUNT_PATHS["output"],
    shell=True,
    check=True,
)
subprocess.run(
    (f"sort -k1,1V -k2,2n {broad_combined_path} " f"> {broad_sorted_path}"),
    cwd=MOUNT_PATHS["output"],
    shell=True,
    check=True,
)
subprocess.run(
    (f"sort -k1,1V -k2,2n {combined_path} " f"> {sorted_path}"),
    cwd=MOUNT_PATHS["output"],
    shell=True,
    check=True,
)

print("Merging files...", flush=True)
subprocess.run(
    f"bedtools merge -c 5 -o mean -i {narrow_sorted_path} > {narrow_merged_path}",
    cwd=MOUNT_PATHS["output"],
    shell=True,
    check=True,
)
subprocess.run(
    f"bedtools merge -c 5 -o mean -i {broad_sorted_path} > {broad_merged_path}",
    cwd=MOUNT_PATHS["output"],
    shell=True,
    check=True,
)
subprocess.run(
    f"bedtools merge -c 5 -o mean -i {sorted_path} > {merged_path}",
    cwd=MOUNT_PATHS["output"],
    shell=True,
    check=True,
)

print("Adding additional columns...", flush=True)
fix_merged_file(narrow_merged_path)
fix_merged_file(broad_merged_path)
fix_merged_file(merged_path)
