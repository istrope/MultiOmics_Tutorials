#!/bin/bash

# Set the CellRanger path (update this with your actual CellRanger installation path)
CELLRANGER_PATH="cellranger-8.0.0/"

# Set the base directory for your Fastq files
FASTQ_DIR="fastq"

# Set the directory for your CellRanger output
OUTPUT_BASE_DIR="cellranger_out"

# Set the transcriptome reference
TRANSCRIPTOME="refdata-gex-GRCh38-2024-A"

# Loop through each sample directory in the FASTQ_DIR
for SAMPLE_DIR in $FASTQ_DIR/*; do
    # Extract the sample ID from the directory name
    SAMPLE_ID=$(basename $SAMPLE_DIR)

    # Set the output directory for the current sample
    OUTPUT_DIR="$OUTPUT_BASE_DIR/${SAMPLE_ID}"

    # Run CellRanger count
    $CELLRANGER_PATH/cellranger count --id=${SAMPLE_ID} \
                                      --transcriptome=$TRANSCRIPTOME \
                                      --fastqs=$SAMPLE_DIR \
                                      --localmem=120 \
                                      --localcores=32 \
                                      --create-bam=false \
                                      --output-dir=$OUTPUT_DIR
done

