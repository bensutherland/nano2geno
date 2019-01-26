#!/bin/bash
# Call SNPs using pysamstats per locus for each sample

# System variables
INPUT_FOLDER="04_mapped"
OUTPUT_FOLDER="06_vcf"

# User-defined variables
THREADS="8"
GENOME="02b_genome/Ots_amplicons_rubias_no-N.fa"
LOCUS="00_archive/Ots_amplicon_names_rubias-SNP_location.txt"


# Align with minimap2 (requires that ref is indexed) 
for file in $(ls -1 $INPUT_FOLDER/*.bam )
do
    # Reporting 
    echo "Calling variants in $file"
    name=$(basename $file)

    # pysamstats to identify all loci
    pysamstats -t variation -f $GENOME $INPUT_FOLDER/$name > $OUTPUT_FOLDER/"${name%.bam}"_var.txt

    # Keep only selected loci
    grep -f $LOCUS $OUTPUT_FOLDER/"${name%.bam}"_var.txt > $OUTPUT_FOLDER/"${name%.bam}"_sel_var.txt

done
