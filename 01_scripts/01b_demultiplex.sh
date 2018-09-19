#!/bin/bash
# Demultiplex the basecalled fastq file 

# System variables
READS_FOLDER="03_basecalled"
DEMULTIPLEXED_FOLDER="03b_demultiplexed"
FORWARD_FOLDER="forward_demultiplexed"
REVERSE_FOLDER="reverse_demultiplexed"

BARCODE_FOLDER="00_archive"
FORWARD_BARCODES="IonCode768.fa"
REVERSE_BARCODES="IonCode768_revcomp.fa"



# User-defined variables
THREADS="7"
FASTQ="all_reads.fastq"


# Demultiplex with forward barcodes 
cutadapt -a file:$BARCODE_FOLDER/$FORWARD_BARCODES \
    --untrimmed-output $DEMULTIPLEXED_FOLDER/$FORWARD_FOLDER/forward_untrimmed.fastq \
    -o $DEMULTIPLEXED_FOLDER/$FORWARD_FOLDER/trimmed-{name}.fastq \
    $READS_FOLDER/$FASTQ

# Secondly, demultiplex the untrimmed with reverse complement barcodes
cutadapt -a file:$BARCODE_FOLDER/$REVERSE_BARCODES \
    --untrimmed-output $DEMULTIPLEXED_FOLDER/$REVERSE_FOLDER/reverse_untrimmed.fastq \
    -o $DEMULTIPLEXED_FOLDER/$REVERSE_FOLDER/trimmed-{name}.fastq \
    $DEMULTIPLEXED_FOLDER/$FORWARD_FOLDER/forward_untrimmed.fastq # this is the file containing untrimmed reads from the first pass.

