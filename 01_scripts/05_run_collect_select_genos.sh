#!/bin/bash
# Runs through samples into Rscript

# System variables
INPUT_FOLDER="06_vcf"

# If running individual samples, for debugging
#sample="none" # testing mode
#Rscript --vanilla 01_scripts/04_collect_select_genos.R $sampl


for file in $(ls -1 $INPUT_FOLDER/*"_var_w_index_selected.txt")
do
    # Reporting
    echo "Working on file $file"
    sample=$( basename $file | awk -F_ '{ print $1 }' - )
    echo "Your sample is $sample"

    # Run Rscript
    Rscript --vanilla 01_scripts/04_collect_select_genos.R $sample 
    
done

