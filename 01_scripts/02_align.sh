#!/bin/bash
# Align the fastq against the reference genome 

# System variables
FASTQ_FOLDER="03_basecalled"

# User-defined variables
THREADS="7"
GENOME="/home/ben/Documents/genomes/GCA_002872995.1_Otsh_v1.0_genomic.fna.gz"
FASTQ="all_reads.fastq"

# Basecall 
minimap2 -ax map-ont $GENOME $FASTQ_FOLDER/$FASTQ > aln_new.sam 


