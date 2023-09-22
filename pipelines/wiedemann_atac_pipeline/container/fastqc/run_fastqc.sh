#!/bin/sh

FASTQC_OPTIONS=""
# If a specific environment variable is set, appends the respective option.
if [[ ! -z "${ADAPTERS}" ]]
then
    FASTQC_OPTIONS="${FASTQC_OPTIONS} --adapters /input/globals/${ADAPTERS}/adapters.txt";
fi

if [[ ! -z "${KMERS}" ]]
then
    FASTQC_OPTIONS="${FASTQC_OPTIONS} --kmers ${KMERS}";
fi

if [[ ! -z "${SVG}" ]] && [[ "${SVG}" == 'true' ]]
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
    mkdir -p /output/$(basename $directory)
    perl -- /FastQC/fastqc --outdir=/output/$(basename $directory) $FASTQC_OPTIONS ${directory}*.fastq.gz
    perl -- /FastQC/fastqc --outdir=/output/$(basename $directory) $FASTQC_OPTIONS ${directory}*.fq.gz
done