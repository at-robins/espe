#!/bin/sh

FASTQC_OPTIONS=""
# If a specific environment variable is set, appends the respective option.
if [[ ! -z "${SVG}" ]]
then
    FASTQC_OPTIONS="${FASTQC_OPTIONS} --svg";
fi

# Prints the specified options.
if [[ ! -z "${FASTQC_OPTIONS}" ]]
then
    echo "Specified options: ${FASTQC_OPTIONS}";
else
    echo Running with default options.
fi

# Iterates over all sample directories and processes them conserving the directory structure.
for directory in /input/samples/*/
do
    perl -- /FastQC/fastqc --outdir=/output/$(basename $directory) $FASTQC_OPTIONS ${directory}*.fastq.gz
done