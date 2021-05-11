#!/bin/bash
# Call SNPs using pysamstats per locus for each sample

# System variables
INPUT_FOLDER="04_mapped"
OUTPUT_FOLDER="06_vcf"

# User-defined variables
THREADS="8"
GENOME="02b_genome/Oki_amps.fa"

# Calling all nucleotide position variants using pysamstats 
for file in $(ls -1 $INPUT_FOLDER/*.bam )
do
    # Reporting 
    echo "Calling variants in $file"
    name=$(basename $file)

    # pysamstats to identify all loci
    pysamstats -t variation -f $GENOME $INPUT_FOLDER/$name > $OUTPUT_FOLDER/"${name%.bam}"_var.txt

done
