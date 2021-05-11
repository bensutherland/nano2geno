#!/bin/bash
# Prepares the locus selection file. Must have the following format:   
# chrom \t pos \t ref \t var 

# User-defined variables
LOCUS_FILE="Oki_amps_name_pos_ref_var.txt"

# System variables (don't change)
INPUT_FOLDER="06_vcf"
CHR_POS_FILE="00_archive/locus_selection_w_index.txt"
CHR_POS_NAME_ONLY="00_archive/locus_selection_w_index_name_only.txt"

# Create new column in locus file for selecting data
awk '{ print $1"_"$2 "\t" $0 }' 00_archive/$LOCUS_FILE > $CHR_POS_FILE

# Create a locus selecting file that has only the index column
# This file is what is used for grep 
awk '{ print $1 }' $CHR_POS_FILE > $CHR_POS_NAME_ONLY 

