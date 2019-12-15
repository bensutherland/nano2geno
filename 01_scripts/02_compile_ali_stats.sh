#!/bin/bash
# Align the fastq against the reference genome 

# System variables
MAPPED_FOLDER="04_mapped/bwa"


# Align with minimap2 (requires that ref is indexed) 
for file in $(ls -1 $MAPPED_FOLDER/*.ali.stats.txt)
do
    # 
   
    echo $(basename $file) >> $MAPPED_FOLDER/alistats.txt
    grep "average quality:" $file >> $MAPPED_FOLDER/alistats.txt

    # 
done
