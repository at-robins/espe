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
DOWNLOAD_PATH_ANNOTATION_GTF = os.path.join(TEMPORARY_OUTPUT_PATH, "annotation_gtf.gz")
DOWNLOAD_PATH_ANNOTATION_GFF = os.path.join(TEMPORARY_OUTPUT_PATH, "annotation_gff.gz")
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


def annotation_url(organism, release, gtf):
    """
    Returns the URL to download the specific annotation file.
    """
    if gtf:
        extension_name = "gtf"
    else:
        extension_name = "gff3"

    return (
        f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_{organism}"
        f"/release_{release}/gencode.v{release}.annotation.{extension_name}.gz"
    )


organism_env = os.environ.get("GLOBAL_ORGANISM")
if organism_env is not None and organism_env == "human":
    print("\tUsing human data...", flush=True)
    genome_id = "GCF_000001405.40"
    genome_name = "GRCh38.p14"
    annotation_release_id = "48"
    annotation_organism_id = "human"
    info_tsv = "info_human.tsv"
else:
    print("\tUsing mouse data...", flush=True)
    genome_id = "GCF_000001635.26"
    genome_name = "GRCm38.p6"
    annotation_release_id = "M25"
    annotation_organism_id = "mouse"
    info_tsv = "info_mouse.tsv"

logging.info("Loading chromosome information...")
genome_info_path = os.path.join(GENOME_INFO_DIRECTORY, info_tsv)
genome_info = {}
with open(genome_info_path, newline="", encoding="utf-8") as tsvfile:
    info_reader = csv.DictReader(tsvfile, dialect="unix", delimiter="\t", quotechar='"')
    for row in info_reader:
        old_chrom_name = row["RefSeq seq accession"]
        new_chrom_name = row["UCSC style name"]
        if old_chrom_name and new_chrom_name:
            logging.info(f"\tRenaming {old_chrom_name} to {new_chrom_name}...")
            genome_info[old_chrom_name] = new_chrom_name

logging.info("Downloading genome...")
command_download_genome = (
    "wget "
    "--no-verbose "
    f"-O {DOWNLOAD_PATH_GENOME} "
    f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{genome_id}/download"
    "?include_annotation_type=GENOME_FASTA"
)
exit_code_download_genome = os.waitstatus_to_exitcode(
    os.system(command_download_genome)
)
if exit_code_download_genome != 0:
    sys.exit(exit_code_download_genome)

logging.info("Extracting genome...")
command_unzip_genome = f"unzip -d {TEMPORARY_OUTPUT_PATH} {DOWNLOAD_PATH_GENOME}"
exit_code_unzip_genome = os.waitstatus_to_exitcode(os.system(command_unzip_genome))
if exit_code_unzip_genome != 0:
    sys.exit(exit_code_unzip_genome)


logging.info("Downloading annotations...")
command_download_annotations = (
    "wget "
    "--no-verbose "
    f"-O {DOWNLOAD_PATH_ANNOTATION_GTF} "
    f"{annotation_url(organism=annotation_organism_id, release=annotation_release_id, gtf=True)} && "
    "wget "
    "--no-verbose "
    f"-O {DOWNLOAD_PATH_ANNOTATION_GFF} "
    f"{annotation_url(organism=annotation_organism_id, release=annotation_release_id, gtf=False)}"
)
exit_code_download_annotations = os.waitstatus_to_exitcode(
    os.system(command_download_annotations)
)
if exit_code_download_annotations != 0:
    sys.exit(exit_code_download_annotations)

logging.info("Extracting annotations...")
command_unzip_annotations = (
    f"gunzip {DOWNLOAD_PATH_ANNOTATION_GTF} && "
    f"gunzip {DOWNLOAD_PATH_ANNOTATION_GFF}"
)
exit_code_unzip_annotations = os.waitstatus_to_exitcode(
    os.system(command_unzip_annotations)
)
if exit_code_unzip_annotations != 0:
    sys.exit(exit_code_unzip_annotations)


logging.info("Renaming files and chromosomes...")
genome_fasta_path_in = os.path.join(
    TEMPORARY_OUTPUT_PATH,
    "ncbi_dataset",
    "data",
    genome_id,
    f"{genome_id}_{genome_name}_genomic.fna",
)
genome_annotations_gff_path_in = DOWNLOAD_PATH_ANNOTATION_GFF.removesuffix(".gz")
genome_annotations_gtf_path_in = DOWNLOAD_PATH_ANNOTATION_GTF.removesuffix(".gz")
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
