#!/bin/bash
# Align the fastq against the reference genome 

# System variables
MAPPED_FOLDER="04_mapped"

# Align with minimap2 (requires that ref is indexed) 
for file in $(ls -1 $MAPPED_FOLDER/*.bam )
do
    # 
    echo "Stats for $file"
    name=$(basename $file)

    # Count alignments per chromosome (chromosome-only, not unanchored)
    samtools view $file | awk '{ print $3 }' - | sort -n | uniq -c | grep -vE 'NW' -  > 05_results/${name%.bam}"_align_per_chr.txt"

done

