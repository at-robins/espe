#!/usr/bin/python
"""This module runs the de novo motif discovery."""

import csv
import json
import math
import multiprocessing
import os
import subprocess

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = next(iter(MOUNT_PATHS["dependencies"].values()))

threads = math.floor(multiprocessing.cpu_count() * 0.8)
if threads < 1:
    threads = 1

print("Loading genome database...")
organism_env = os.environ.get("GLOBAL_ORGANISM")
if organism_env is not None and organism_env == "human":
    print("\tUsing human data...", flush=True)
    genome_id = "hg38"
else:
    print("\tUsing mouse data...", flush=True)
    genome_id = "mm10"

print("Searching differentially accessibility analysis output...", flush=True)
for root, dirs, files in os.walk(INPUT_FOLDER):
    for file in files:
        if file.startswith("annotated"):
            file_path_annotated = os.path.join(root, file)
            print(f"Detected {file_path_annotated}")
            print("\tPre-processing...", flush=True)
            output_directory_name = file.removeprefix(
                "annotated_differential_accessibility_"
            ).removesuffix(".csv")
            directory_path_output = os.path.join(
                MOUNT_PATHS["output"], output_directory_name
            )
            os.makedirs(directory_path_output, exist_ok=True)
            file_path_peaks = os.path.join(
                directory_path_output, "differential_accessible_peaks.tsv"
            )
            # Creates the input file with all differentially accessible peaks.
            with open(
                file_path_annotated, newline="", mode="rt", encoding="utf-8"
            ) as annotated_file:
                with open(
                    file_path_peaks, newline="", mode="wt", encoding="utf-8"
                ) as peak_file:
                    annotated_reader = csv.DictReader(
                        annotated_file, dialect="unix", delimiter=",", quotechar='"'
                    )
                    peak_writer = csv.writer(
                        peak_file,
                        dialect="unix",
                        delimiter="\t",
                        quotechar='"',
                        quoting=csv.QUOTE_MINIMAL,
                    )
                    # Writes the header.
                    peak_writer.writerow(
                        [
                            "PeakID",
                            "Chr",
                            "Start",
                            "End",
                            "Strand",
                        ]
                    )
                    for row in annotated_reader:
                        # Copies relevant values of each row.
                        if row["padj"] != "" and float(row["padj"]) <= 0.05:
                            peak_writer.writerow(
                                [
                                    row["PeakID"],
                                    row["Chr"],
                                    row["Start"],
                                    row["End"],
                                    row["Strand"],
                                ]
                            )

            print("\tRunning motif discovery...", flush=True)
            subprocess.run(
                (
                    f"findMotifsGenome.pl {file_path_peaks} {genome_id} "
                    f"{directory_path_output} -size given -mask -p {threads}"
                ),
                cwd=directory_path_output,
                shell=True,
                check=True,
            )
