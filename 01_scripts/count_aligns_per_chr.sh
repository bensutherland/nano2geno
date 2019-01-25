#!/bin/bash
# Count number alignments per chromosome

# System variables
MAPPED_FOLDER="04_mapped"

# User-defined variables
CHR_ID="linkage_group_LG"

# Count alignments per chromosome with samtools
for file in $(ls -1 $MAPPED_FOLDER/*.bam )
do
    # Reporting per sample 
    echo "Stats for $file"
    name=$(basename $file)

    # Count alignments per chromosome (chromosome-only, not unanchored)
    samtools view $file | awk '{ print $3 }' - | sort -n | uniq -c | grep "$CHR_ID" -  > 05_results/${name%.bam}"_align_per_chr.txt"

done

