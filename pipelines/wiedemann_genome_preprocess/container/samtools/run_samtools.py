#!/usr/bin/python
"""This module runs the samtools indexing."""

import json
import logging
import os
import sys

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))

# Setup.
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")

logging.info("Creating indices...")
samtools_command = ("samtools faidx "
                    f"{MOUNT_PATHS['dependencies']['get_genome']}/genome.fa "
                    f"-o {MOUNT_PATHS['output']}/genome.fa.fai")

exit_code = os.waitstatus_to_exitcode(os.system(samtools_command))
if exit_code != 0:
    sys.exit(exit_code)

logging.info("Calculating genome sizes...")
cut_command = ("cut -f 1,2 "
                f"{MOUNT_PATHS['output']}/genome.fa.fai "
                f"> {MOUNT_PATHS['output']}/genome.chrom.sizes")

exit_code = os.waitstatus_to_exitcode(os.system(cut_command))
if exit_code != 0:
    sys.exit(exit_code)

logging.info("Done.")
