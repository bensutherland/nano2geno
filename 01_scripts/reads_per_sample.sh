#!/bin/bash
# Join the two types of barcode-called data, forward and reverse 

# System variables
SAMPLES_FOLDER="03b_demultiplexed/all_decat_cat"

# User-defined variables
THREADS="7"

for i in $(ls $SAMPLES_FOLDER/*.fastq) ; do echo $i ; awk '{s++}END{print s/4}' $i ; done > reads_per_sample.txt

# Put into good format
awk '{printf "%s%s", $0, (NR%2?FS:RS)}' reads_per_sample.txt > reads_per_sample2.txt

# Reporting
echo "### Total reads in samples ###"
grep -v 'untrimmed' reads_per_sample2.txt | awk '{ SUM += $2 } END { print SUM }'

echo "### Total reads unassigned ###"
grep 'reverse_untrimmed.fastq'  reads_per_sample2.txt | awk '{ print $2 }' - 
