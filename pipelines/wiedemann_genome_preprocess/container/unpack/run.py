#!/usr/bin/python
"""This module downloads the genome files."""

import csv
import json
import logging
import os
import shutil
import sys

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
TEMPORARY_OUTPUT_PATH = os.path.join(MOUNT_PATHS["output"], "tmp")
DOWNLOAD_PATH_GENOME = os.path.join(TEMPORARY_OUTPUT_PATH, "genome.zip")
GENOME_INFO_DIRECTORY = "/genome_info"

# Setup.
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
os.makedirs(TEMPORARY_OUTPUT_PATH, exist_ok=True)


def rename_chromosomes(input_path, output_path, rename_info):
    """
    Renames all instances of chromoses based on the specified information and saves the
    file to the specified output path
    """
    logging.info(f"Renaming {input_path} to {output_path}...")
    with open(input_path, mode="rt", encoding="utf-8") as file_in:
        with open(output_path, mode="wt", encoding="utf-8") as file_out:
            for line in file_in:
                replaced_line = line
                for original, replacement in rename_info.items():
                    replaced_line = replaced_line.replace(original, replacement)
                file_out.write(replaced_line)


logging.info("Loading chromosome information...")
genome_info_path = os.path.join(GENOME_INFO_DIRECTORY, "info_mouse.tsv")
genome_info = {}
with open(genome_info_path, newline="", encoding="utf-8") as tsvfile:
    info_reader = csv.DictReader(tsvfile, dialect="unix", delimiter="\t", quotechar='"')
    for row in info_reader:
        old_chrom_name = row["RefSeq seq accession"]
        new_chrom_name = row["UCSC style name"]
        if old_chrom_name and new_chrom_name:
            logging.info(f"\tRenaming {old_chrom_name} to {new_chrom_name}...")
            genome_info[old_chrom_name] = new_chrom_name

logging.info("Downloading genome and annotations...")
command_download_genome = (
    "wget "
    "--no-verbose "
    f"-O {DOWNLOAD_PATH_GENOME} "
    "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000001635.26/download"
    "?include_annotation_type=GENOME_FASTA,GENOME_GFF,GENOME_GTF"
)
exit_code_download_genome = os.waitstatus_to_exitcode(
    os.system(command_download_genome)
)
if exit_code_download_genome != 0:
    sys.exit(exit_code_download_genome)

logging.info("Extracting genome and annotations...")
command_unzip_genome = f"unzip -d {TEMPORARY_OUTPUT_PATH} {DOWNLOAD_PATH_GENOME}"
exit_code_unzip_genome = os.waitstatus_to_exitcode(os.system(command_unzip_genome))
if exit_code_unzip_genome != 0:
    sys.exit(exit_code_unzip_genome)

logging.info("Renaming files and chromosomes...")
genome_fasta_path_in = os.path.join(
    TEMPORARY_OUTPUT_PATH,
    "ncbi_dataset",
    "data",
    "GCF_000001635.26",
    "GCF_000001635.26_GRCm38.p6_genomic.fna",
)
genome_annotations_gff_path_in = os.path.join(
    TEMPORARY_OUTPUT_PATH,
    "ncbi_dataset",
    "data",
    "GCF_000001635.26",
    "genomic.gff",
)
genome_annotations_gtf_path_in = os.path.join(
    TEMPORARY_OUTPUT_PATH,
    "ncbi_dataset",
    "data",
    "GCF_000001635.26",
    "genomic.gtf",
)
genome_fasta_path_out = os.path.join(MOUNT_PATHS["output"], "genome.fa")
genome_annotations_gff_path_out = os.path.join(MOUNT_PATHS["output"], "annotations.gff")
genome_annotations_gtf_path_out = os.path.join(MOUNT_PATHS["output"], "annotations.gtf")
rename_chromosomes(
    input_path=genome_fasta_path_in,
    output_path=genome_fasta_path_out,
    rename_info=genome_info,
)
rename_chromosomes(
    input_path=genome_annotations_gff_path_in,
    output_path=genome_annotations_gff_path_out,
    rename_info=genome_info,
)
rename_chromosomes(
    input_path=genome_annotations_gtf_path_in,
    output_path=genome_annotations_gtf_path_out,
    rename_info=genome_info,
)

logging.info("Cleaning up temporary files...")
shutil.rmtree(TEMPORARY_OUTPUT_PATH)
logging.info("Done.")
