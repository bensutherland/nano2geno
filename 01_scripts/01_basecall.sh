#!/bin/bash
# Basecall folders within the 02_raw_data folder

# System variables
INPUT_FOLDER="02_raw_data"
OUTPUT_FOLDER="03_basecalled"

# User-defined variables
THREADS="7"
FLOWCELL="FLO-MIN106"
KIT="SQK-LSK108"


# Basecall 
read_fast5_basecaller.py \
    --input $INPUT_FOLDER \
    --worker_threads $THREADS \
    --save_path $OUTPUT_FOLDER \
    --flowcell $FLOWCELL \
    --kit $KIT \
    -r \
    --barcoding \
    -o fastq

