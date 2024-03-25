#!/usr/bin/python
"""This module downloads the genome files."""

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

logging.info("Downloading genome and annotations...")
command_download_genome = (
    "wget "
    "--no-verbose "
    f"-O {DOWNLOAD_PATH_GENOME} "
    "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000001635.26/download"
    "?include_annotation_type=GENOME_FASTA,GENOME_GFF"
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
genome_annotations_path_oin = os.path.join(
    TEMPORARY_OUTPUT_PATH,
    "ncbi_dataset",
    "data",
    "GCF_000001635.26",
    "genomic.gff",
)
genome_fasta_path_out = os.path.join(MOUNT_PATHS["output"], "genome.fa")
genome_annotations_path_out = os.path.join(MOUNT_PATHS["output"], "annotations.gff")

with open("log.txt") as infile:
    for line in infile:
        print(line)

logging.info("Renaming chromosomes...")
genome_info_path = os - path.join(GENOME_INFO_DIRECTORY, "chromosome_names_mouse.tsv")
with open(genome_info_path, newline="", encoding="utf-8") as tsvfile:
    info_reader = csv.DictReader(tsvfile, dialect="unix", delimiter="\t", quotechar='"')
    for row in info_reader:
        old_chrom_name = row["RefSeq seq accession"]
        new_chrom_name = row["UCSC style name"]
        if old_chrom_name and new_chrom_name:
            logging.info("\tRenaming {old_chrom_name} to {new_chrom_name}...")

        sample_type = sample_info_map.get(row["sample type"])
        sample_directory = row["sample directory"]
        if sample_type is None:
            sample_info_map[row["sample type"]] = [sample_directory]
        else:
            sample_type.append(sample_directory)

logging.info("Compressing files...")
command_zip = f"gzip -k {os.path.join(MOUNT_PATHS['output'], 'genome.fa')}"
exit_code_zip = os.waitstatus_to_exitcode(os.system(command_zip))
if exit_code_zip != 0:
    sys.exit(exit_code_zip)

logging.info("Cleaning up temporary files...")
shutil.rmtree(TEMPORARY_OUTPUT_PATH)
logging.info("Done.")
