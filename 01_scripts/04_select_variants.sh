#!/bin/bash
# Add index columns, create a hotspot name file for selection purposes
# Selects the specific loci out of a variant file (pysamstats output)

# System variables
INPUT_FOLDER="06_vcf"

# User-defined variables
THREADS="8"
GENOME="02b_genome/Ots_amplicons_rubias_no-N.fa"
CHR_POS_FILE="00_archive/locus_selection_w_index.txt"
CHR_POS_NAME_ONLY="00_archive/locus_selection_w_index_name_only.txt"

# Add an index column to the variant file 
for file in $( ls -1 $INPUT_FOLDER/*_var.txt )
do
    # Reporting 
    echo "Making unique ID column $file"
    name=$(basename $file)

    # Make a column in every _var.txt file that can be used for data selection
    awk '{ print $1"_"$2 "\t" $0 }' $INPUT_FOLDER/$name > $INPUT_FOLDER/"${name%.txt}"_w_index.txt

done

# Do selection on datafiles
for file in $( ls -1 $INPUT_FOLDER/*_w_index.txt )
do
    # Reporting
    echo "Selecting lines out of $file"
    name=$(basename $file)

    # Use the locus selection file to only retain selected loci
    cat $CHR_POS_NAME_ONLY | 
        awk '{ print $1"\t" }' - |
        xargs -I{} grep -P {} $INPUT_FOLDER/$name > $INPUT_FOLDER/"${name%.txt}"_selected.txt 

done

