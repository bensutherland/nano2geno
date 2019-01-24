#!/bin/bash
# Align the fastq against the reference genome 

# System variables
SAMPLE_FOLDER="04_samples"
MAPPED_FOLDER="04_mapped"

# User-defined variables
THREADS="7"
GENOME="02b_genome/Ots_amplicons_rubias_no-N.fa"

# Align with marginAlign 
for file in $(ls -1 $SAMPLE_FOLDER/*.fastq )
do
    # Initializing 
    echo "Aligning file $file"
    name=$(basename $file)

    # Align 
    marginAlign $SAMPLE_FOLDER/$name $GENOME $MAPPED_FOLDER/"${name%.fastq}".sam --jobTree $SAMPLE_FOLDER/jobTree_"${name%.fastq}" --maxThreads=$THREADS

    # Collect stats on the alignment
    samtools stats $MAPPED_FOLDER/"${name%.fastq}".sam > $MAPPED_FOLDER/"${name%.fastq}"_ali_stats.txt
    grep "^COV" $MAPPED_FOLDER/"${name%.fastq}"_ali_stats.txt > $MAPPED_FOLDER/"${name%.fastq}"_cov_stats.txt

done

#todo remove jobtree files

