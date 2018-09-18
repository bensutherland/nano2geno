#!/bin/bash
# Join the two types of barcode-called data, forward and reverse 

# System variables
READS_FOLDER="03_basecalled"
DEMULTIPLEXED_FOLDER="03b_demultiplexed"
FORWARD_FOLDER="forward_demultiplexed"
REVERSE_FOLDER="reverse_demultiplexed"

SAMPLES_FOLDER="04_samples"

BARCODE_FOLDER="00_archive"
FORWARD_BARCODES="IonCode768.fa"
REVERSE_BARCODES="IonCode768_revcomp.fa"

# User-defined variables
THREADS="7"
FASTQ="all_reads.fastq"

# Merge samples
for i in $(ls $DEMULTIPLEXED_FOLDER/$FORWARD_FOLDER/) ; do cat $DEMULTIPLEXED_FOLDER/$FORWARD_FOLDER/$i $DEMULTIPLEXED_FOLDER/$REVERSE_FOLDER/$i > $SAMPLES_FOLDER/$i ; done

# Keep only the final untrimmed, not the first round untrimmed
rm $SAMPLES_FOLDER/untrimmed.fastq
cp $DEMULTIPLEXED_FOLDER/$REVERSE_FOLDER/reverse_untrimmed.fastq $SAMPLES_FOLDER
