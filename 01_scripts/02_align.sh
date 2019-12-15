#!/bin/bash
# Align the fastq against the reference genome 

# System variables
SAMPLE_FOLDER="04_samples"
MAPPED_FOLDER="04_mapped"

# User-defined variables
THREADS="8"
#GENOME="/home/ben/Documents/genomes/GCF_002872995.1_Otsh_v1.0_genomic.fna"
#GENOME="/home/ben/Documents/01_nanopore/nano2geno/02b_genome/GCF_002021735.1_Okis_V1_genomic.fna"
GENOME="02b_genome/Oki_amps.fa"

# Align with minimap2 (requires that ref is indexed) 
for file in $(ls -1 $SAMPLE_FOLDER/*.fastq )
do
    # 
    echo "Aligning file $file"
    name=$(basename $file)

    # Align fastq
    bwa mem $GENOME $SAMPLE_FOLDER/$name | 
        samtools view -bS -q 2 - | 
        samtools sort -o $MAPPED_FOLDER/"${name%.fastq}".bam -

    # Index bam
    samtools index $MAPPED_FOLDER/"${name%.fastq}".bam

    # Collect stats on the alignment
    samtools stats $MAPPED_FOLDER/"${name%.fastq}".bam > $MAPPED_FOLDER/"${name%.fastq}".ali.stats.txt
    grep "^COV" $MAPPED_FOLDER/"${name%.fastq}".ali.stats.txt > $MAPPED_FOLDER/"${name%.fastq}".cov.stats.txt

    # 
done
