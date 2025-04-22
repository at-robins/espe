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
MERGED_MOTIF_FILE_NAME = "merged.motif"
MOTIF_ANNOTATED_PEAK_FILE_NAME = "motif_annotated_peaks.tsv"
PEAK_REGION_FILE_NAME = "differentially_accessible_peaks.tsv"

threads = math.floor(multiprocessing.cpu_count() * 0.8)
if threads < 1:
    threads = 1


def merge_motifs(motif_input_folder):
    """
    Merges motif files.
    """
    print("\tMerging motif files...", flush=True)
    for root, dirs, files in os.walk(motif_input_folder):
        for file in files:
            if file.endswith(".motif"):
                subprocess.run(
                    (
                        f'cat "{os.path.join(root, file)}" >> '
                        f'"{os.path.join(root, MERGED_MOTIF_FILE_NAME)}"'
                    ),
                    cwd=motif_input_folder,
                    shell=True,
                    check=True,
                )


def filter_p(pvalue):
    """
    Applies the p-value filter.
    """
    return pvalue != "" and float(pvalue) <= 0.05


def fitler_lfc(lfc, filter):
    """
    Applies the log-2-fold-change filter if set.
    """
    return (
        filter is None
        or (lfc != "" and float(lfc) >= 0.0 and filter == "+")
        or (lfc != "" and float(lfc) <= 0.0 and filter == "-")
    )


def create_region_file(input_file_path, output_directory_path, filter=None):
    if filter == "+":
        output_directory_path_final = os.path.join(
            output_directory_path, "more accessible regions"
        )
    elif filter == "-":
        output_directory_path_final = os.path.join(
            output_directory_path, "less accessible regions"
        )
    else:
        output_directory_path_final = os.path.join(output_directory_path, "all regions")

    os.makedirs(output_directory_path_final, exist_ok=True)
    file_path_peaks = os.path.join(output_directory_path_final, PEAK_REGION_FILE_NAME)

    # Creates the input file with all differentially accessible peaks.
    with open(
        input_file_path, newline="", mode="rt", encoding="utf-8"
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
                if filter_p(pvalue=row["padj"]) and fitler_lfc(
                    lfc=row["log2FoldChange"], filter=filter
                ):
                    peak_writer.writerow(
                        [
                            row["PeakID"],
                            row["Chr"],
                            row["Start"],
                            row["End"],
                            row["Strand"],
                        ]
                    )


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
            create_region_file(
                input_file_path=file_path_annotated,
                output_directory_path=directory_path_output,
                filter=None,
            )
            create_region_file(
                input_file_path=file_path_annotated,
                output_directory_path=directory_path_output,
                filter="+",
            )
            create_region_file(
                input_file_path=file_path_annotated,
                output_directory_path=directory_path_output,
                filter="-",
            )

for directory_path_output, dirs, files in os.walk(MOUNT_PATHS["output"]):
    for file in files:
        if file == PEAK_REGION_FILE_NAME:
            file_path_peaks = os.path.join(directory_path_output, file)
            print(f"Running motif discovery for {file_path_peaks}...", flush=True)
            subprocess.run(
                (
                    f'findMotifsGenome.pl "{file_path_peaks}" {genome_id} '
                    f'"{directory_path_output}" -size given -mask -p {threads}'
                ),
                cwd=directory_path_output,
                shell=True,
                check=True,
            )

            print("\tMerging known motifs...", flush=True)
            directory_known_motifs = os.path.join(directory_path_output, "knownResults")
            file_known_motifs = os.path.join(
                directory_known_motifs, MERGED_MOTIF_FILE_NAME
            )
            file_known_motifs_annotated = os.path.join(
                directory_known_motifs, MOTIF_ANNOTATED_PEAK_FILE_NAME
            )
            merge_motifs(directory_known_motifs)

            print("\tMerging de novo motifs...", flush=True)
            directory_denovo_motifs = os.path.join(
                directory_path_output, "homerResults"
            )
            file_denovo_motifs = os.path.join(
                directory_denovo_motifs, MERGED_MOTIF_FILE_NAME
            )
            file_denovo_motifs_annotated = os.path.join(
                directory_denovo_motifs, MOTIF_ANNOTATED_PEAK_FILE_NAME
            )
            merge_motifs(directory_denovo_motifs)

            if os.path.isfile(file_known_motifs):
                print("\tAnnotating peaks with discovered known motifs...", flush=True)
                subprocess.run(
                    (
                        f'findMotifsGenome.pl "{file_path_peaks}" {genome_id} '
                        f'"{directory_path_output}" -size given -mask -p {threads} '
                        f'-find "{file_known_motifs}" > "{file_known_motifs_annotated}"'
                    ),
                    cwd=directory_path_output,
                    shell=True,
                    check=True,
                )
            else:
                print(
                    "\tNo known motifs discovered. Skipping peak annotation...",
                    flush=True,
                )

            if os.path.isfile(file_denovo_motifs):
                print(
                    "\tAnnotating peaks with discovered de novo motifs...", flush=True
                )
                subprocess.run(
                    (
                        f'findMotifsGenome.pl "{file_path_peaks}" {genome_id} '
                        f'"{directory_path_output}" -size given -mask -p {threads} '
                        f'-find "{file_denovo_motifs}" > "{file_denovo_motifs_annotated}"'
                    ),
                    cwd=directory_path_output,
                    shell=True,
                    check=True,
                )
            else:
                print(
                    "\tNo de novo motifs discovered. Skipping peak annotation...",
                    flush=True,
                )
