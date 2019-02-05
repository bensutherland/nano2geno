#!/bin/bash
# Runs through samples into Rscript

# User-defined variables
#sample="none"
sample="BC10"

Rscript --vanilla 01_scripts/04_collect_select_genos.R $sample 

