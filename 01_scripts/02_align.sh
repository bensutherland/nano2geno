#!/bin/bash
# Align the fastq against the reference genome 

# System variables
FASTQ_FOLDER="03_basecalled"
MAPPED_FOLDER="04_mapped"

# User-defined variables
THREADS="7"
GENOME="/home/ben/Documents/genomes/GCF_002872995.1_Otsh_v1.0_genomic.fna"
FASTQ="all_reads.fastq"
BAM="all_reads.sorted.bam"

# Align with minimap2 (requires that ref is indexed) 
minimap2 -ax map-ont $GENOME $FASTQ_FOLDER/$FASTQ | 
    samtools view -bS - | samtools sort -o $MAPPED_FOLDER/$BAM - 

# Index the output with samtools
samtools index $MAPPED_FOLDER/$BAM
