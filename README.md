# nano2geno
SNP genotyping and stock ID from nanopore data after concatenation of amplicon libraries     
C: Ben J. G. Sutherland and Christoph Deeg     

Currently in **Development mode only**     
Disclaimer: This is pipeline comes as is, with no guarantees. At the moment it is purely for the authors' use to explore the suitability of in-field SNP genotyping     

#### Requirements and stock ID using the nanopore platform.       
Minknow: Data collection and basecalling (alternatively use Albacore to create fastq files from fast5)     
Porechop https://github.com/rrwick/Porechop      
minimap2    
samtools    
pysamstats    
R: Packages: rubias, tidyverse, ggplot2    

#### Inputs
Basecalled nanopore reads in fastq format. In our case we are working with concatenated barcoded reads:     
Amplicons are barcoded with ONT 96 barcode kit (individual ID barcodes), these amplicons are linked to each other using a custom concatenation adapter.     
In the future, these reads can be barcoded again with a secondary barcode to identify different libraries (library ID barcodes) to run more than 96 individuals on the same flow cell.     
Reference data: Reference genome or reference amplicon sequences. With defined SNP loci with known reference and alternate nucleotides in a "chromosome position file".    
Barcodes and concatenation adapter sequences: ONT barcodes included in "adapters.py" provided by porechop. Add custom adapters and barcodes if needed.    
Rubias baseline data of SNP distribution in reference populations    


### 1. Prepare sequencing data
Copy files from ONT sequencing run folder fastq_pass into 03_basecalled.       
`cp xxx/fastq_pass/*.fastq 03_basecalled`    

Combine all data from a single run:      
`cat 03_basecalled/*.fastq > 03_basecalled/all_reads.fastq`      

How many reads are present?    
`grep -cE '^\+$' 03_basecalled/all_reads.fastq`     

If using limited resources (e.g. laptop in field), split reads into subfiles in subfolder (sub_x) to preserve RAM during porechop deconcatenation and demultiplexing.
In our case, Ubuntu 14.06, 31.2 GiB RAM 7700K CPU @ 4.20GHz × 8, files with read numbers exeeding 550k need to be split into sub folders of aproximately 100 files (4000 reads per file). 
`ls -Q fastq_pass/ | head -100 | xargs -i mv fastq_pass/{} fastq_pass/sub_1/`     

Then concatenate as above:     
`cat 03_basecalled/sub_x/*.fastq > 03_basecalled/sub_x/all_sub_x_reads.fastq`     

If subsetting, adjust accordingly below (e.g. run commands in each sub folder).    


### 2. Data quality check
Run fastqc on this data:    
```
fastqc 03_basecalled/all_reads.fastq -o 03_basecalled/fastqc -t 5   
multiqc -o 03_basecalled/fastqc/ 03_basecalled_fastqc    
```

### 3. De-multiplex and deconcatenation with Porechop
Disclaimer: We are currently using different adapter files (adapters-1.py, adapters-2.py, etc; within porechop) and rename (adapters.py) before each step (might get automated in future).     

#### 3a. De-multiplex using external library barcodes (optional, if sequencing >96 samples)
First we want to de-multiplex concatenated reads by library, according to the external library ID barcode.      
To do this, edit the adapter file to remove all of the adapters that are not used as external library ID adapters.      
To exclued all individual ID barcodes is paramount to avoid porechop from discarding concatenated reads due to middle adapters.      

To prepare for de-multiplex step one, rename adapters-1.py (containing only the library ID barcodes) to adapters.py     
`mv ~/Programs/Porechop/porechop/adapters-1.py ~/Programs/Porechop/porechop/adapters.py`     

Demultiplex using the external library ID barcodes:      
`~/Programs/Porechop/porechop-runner.py -i 03_basecalled/all_reads.fastq -b 03a_demultiplex_library -t 8 --adapter_threshold 90 --end_threshold 75 --extra_trim_end 0`       
The results will be in the folder specified by -b. (here `03a_demultiplex_library`)      
This will bin the concatenated reads by library ID barcode after trimming the barcodes from the end.       

Then rename the used adapter file back to it's original name:    
`mv ~/Programs/Porechop/porechop/adapters.py ~/Programs/Porechop/porechop/adapters-1.py`     


#### 3b. De-concatenate library using concatenation adapter 
To prepare for individual ID de-multiplex (step 3c), rename adapters-2.py to adapters.py. Adapters-2.py contains only the custom concatenation adapter.
`mv ~/Programs/Porechop/porechop/adapters-2.py ~/Programs/Porechop/porechop/adapters.py`     

Then, run porechop without the binning option to deconcatenate and inflate the read number (RAM intensive!)
`~/Programs/Porechop/porechop-runner.py -i 03a_demultiplex_library/BC96.fastq -o 03a_demultiplex_library/BC96_deconcat.fastq -t 16 --middle_threshold 75 --min_split_read_size 100 --extra_middle_trim_bad_side 0 --extra_middle_trim_good_side 0`     
Then rename the used adapters.py file back to it's original name:    
`mv ~/Programs/Porechop/porechop/adapters.py ~/Programs/Porechop/porechop/adapters-2.py`     

#### 3c. De-multiplex individuals from the de-concatenated library
To prepare for the final de-multiplexing and binning step (de-multiplex by individual ID barcodes), rename adapters-3.py to adapters.py. Adapters-3.py is the standard adapters.py file provided by porechop and contains all barcodes used by ONT.
`mv ~/Programs/Porechop/porechop/adapters-3.py ~/Programs/Porechop/porechop/adapters.py`     

Now run porechop to de-multiplex and bin samples based on all 96 barcodes. Search 100000 reads to identify all used barcodes.      
`~/Programs/Porechop/porechop-runner.py -i 03a_demultiplex_library/BC96_deconcat.fastq -b 03b_demultiplexed/ -t 16 --adapter_threshold 90 --end_threshold 75 --check_reads 100000`      

Then rename the used barcode file back to its original name:     
`mv ~/Programs/Porechop/porechop/adapters.py ~/Programs/Porechop/porechop/adapters-3.py`     


#### 3d. Conflate subfolders (optional)
If using subfolders to preserve RAM, run "01c_cat_subsets.sh" to compile all binned reads with the same barcode into combined files (adjust 01c_cat_subsets.sh to fit number of subdirectories).     
`mkdir 03b_demultiplexed/all_decat_cat`
`./01_scripts/01c_cat_subsets.sh`

#### 3e. Prepare files for next steps
Move the demultiplexed and binned reads into the appropriate folder.    
`cp -l 03b_demultiplexed/*.fastq 04_samples`      

### 4. Align against reference regions
The demultiplexed fastq files are in `04_samples`.     

Index your reference genome or amplicon sequence file using minimap2    
`minimap2 -d Ots_subset.mmi 02b_genome/genome_or amplicon_ref.fasta`    

Set the genome variable and use minimap2 and samtools to align fastq samples in `04_samples` against the reference genome or amplicon sequence files:    
`./01_scripts/02_align.sh`      
This will produce stat files as well as indexed and sorted bam files.    

### 5. Call variants 
#### 5.a. Identify nucleotides at all loci with pysamstats
This will use the reference and all bam files in `04_mapped`.     
Set the reference variable to the appropriate file and run pysamstats and output the nucleotides found at every position in all aligned sample files.    
`01_scripts/03_call_SNP_pysamstats.sh`      

#### 5.b. Prepare chromosome selection file 
This step will use a chromosome position file in `00_archive` that specifies the chromosome and the position of the SNP, as well as the according nucleotides needed for the analysis.   
The chromosome position file will have the following format:     
`chrom \t pos \t ref \t var`       

Run the following script to prepare this chromosome position file into the formats needed for the analysis:     
`./01_scripts/04a_prepare_locus_selection_file.sh`    

This will produce a file called `00_archive/locus_selection_w_index.txt` that will be used later. 

#### 5.c. Select loci of interest
Use the following script to add an index column to the data file and then select out the specific lines of interest containing the SNPs from the data files.     
`./01_scripts/04_select_variants.sh`    

This will output files entitled `<sample>_var_w_index_selected.txt` to be used downstream.    

#### 6. Collect specific SNPs 
Use the Rscript interactively to adjust thresholds for nucleotide calling and hete-cuttoff as needed (e.g. 10x and 25%):    
`04_collect_select_genos.R`     

If you want to automatically run all samples, use the shell script, which will operate on all the file in `06_vcf/*"_var_w_index_selected.txt"`:      
`./01_scripts/05_run_collect_select_genos.sh`    
Note that this will also require that you have the guide file in `00_archive/locus_selection_w_index.txt`       

The results will end up in `07_rubias` in wide format, one file for each sample, with titles `<sample>_rubias_input.txt`      

Then run the following script to collect all input data and put into a single dataframe      
`01_scripts/06_merge_all_rubias_input_files.R`

#### 7. Rubias

Now run the rubias R script for stock ID `01_scripts/07_Rubias_GSI-script-NA.R`      
