#!/bin/bash
# Demultiplex the basecalled fastq file 

# System variables
READS_FOLDER="03_basecalled"

# User-defined variables
THREADS="7"

# Demultiplex 



read_fast5_basecaller.py \
    --input $INPUT_FOLDER \
    --worker_threads $THREADS \
    --save_path $OUTPUT_FOLDER \
    --flowcell $FLOWCELL \
    --kit $KIT \
    -r \
    --barcoding \
    -o fastq

