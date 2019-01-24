#!/bin/bash
# Using the alignments, call SNPs using marginCaller

# System variables
MAPPED_FOLDER="04_mapped"
VCF_FOLDER="06_vcf"

# User-defined variables
THREADS="7"
GENOME="02b_genome/Ots_amplicons_rubias_no-N.fa"

# Use the alignment and produce a vcf
for file in $(ls -1 $MAPPED_FOLDER/*.sam )
do
    # Initializing 
    echo "Calling SNPS in file $file"
    name=$(basename $file)

    # Align 
    marginCaller $MAPPED_FOLDER/$name $GENOME $VCF_FOLDER/"${name%.sam}".vcf --jobTree $VCF_FOLDER/jobTree_"${name%.sam}" --maxThreads=$THREADS

done
