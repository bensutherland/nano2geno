# nano2geno
Genotyping from nanopore data

Currently in **Development mode only**
Disclaimer: this is a simple pipeline that comes with no guarantees. At the moment it is purely for the author's use to better understand nanopore data.   

#### Requirements
Albacore basecaller     
Porechop https://github.com/rrwick/Porechop      
minimap2    
samtools    

#### Inputs
fast5     
genome.fa or reference sequences of amplicons   
inner barcodes (.csv or .fa)    

### 1. Basecalling
If you have fast5 files to basecall, move the folder of fast5 files to `02_raw_data`. Note that fast5 files will be needed to index the fastq for SNP calling later.     
Note: for just getting stats on a fastq file without doing genotyping, collect fastq files and put in `03_basecalled/all_reads.fastq`, and skip to the next step. (e.g. `cp -r ../raw_data/PBT_trial*/*.fastq ./03_basecalled` then cat these files together into a fastq)         

Run `read_fast5_basecaller.py` recursively, with demultiplexing of nanopore barcodes on fast5 files within `02_raw_data` use the following:    
`./01_scripts/01_basecall.sh`     

Output will be in `03_basecalled/workspace`, with reads inside folders `pass/fail/calibration_strands`, and reads in the `Pass` folder in different barcode subfolders containing fastq files for each barcode.    

If you have an external library barcode, or don't care what barcode was on your library, combine the fastq files into `03_basecalled/all_reads.fastq`     
`find ./03_basecalled/workspace/pass -name "*.fastq" | while read file ; do cat $file >> 03_basecalled/all_reads.fastq ; done`
(make sure to remove any previous version of all_reads.fastq or else this will combine all to this)   
(#todo: describe how to split this into sections to run in parallel)    

(# new instructions #) 
Combine all data from a single run 
`cat 03_basecalled/*.fastq > 03_basecalled/all_reads.fastq`      

If read number exeeds 550k, split into subfiles in subfolder (sub_x) of aproximately 100 files (4000 reads per file) to preserve RAM during porechop demultiplexing. 
`ls -Q fastq_pass/ | head -100 | xargs -i mv fastq_pass/{} fastq_pass/sub_1/`
Run porechop accordingly on subfiles

How many reads are present?    
`grep -cE '^\+$' 03_basecalled/all_reads.fastq`     

### 2. Data quality check
Run fastqc on this data:    
`fastqc 03_basecalled/all_reads.fastq -o 03_basecalled/fastqc -t 5`   
`multiqc -o 03_basecalled/fastqc/ 03_basecalled_fastqc`    

### 3A. Demultiplexing the inner library with cutadapt
#### a. Create fasta from .csv, and also make a reverse-complement
My data is kept in `IonCode768.csv`     

Custom: make csv file a fasta with named adapters:    
`grep -vE '^index' ./IonCode768.csv | awk -F"," '{ print ">"$2 "\n" $3 $5 "\n" }' -  | grep -vE '^$' - > IonCode768.fa`

Make reverse complement of barcode file (from Pierre Lindenbaum, Biostars : https://www.biostars.org/p/189325/)   
`cat IonCode768.fa | while read L; do echo $L; read L; echo "$L" | rev | tr "ATGC" "TACG" ; done | sed -e "s/^M//" > IonCode768_revcomp.fa`      

#### b. Demultiplex
Demultiplex using the forward adapter, then on the unidentified file from the forward adapter, demultiplex again using the reverse complement adapter file:     
`./01_scripts/01b_demultiplex.sh`     
(#todo: this doesn't yet use parallel effectively, could implement w/ stacks_workflow cutadapt parallel approach, but this requires different input folders, so would have to break up into sections, or do each library separately (probably the best bet)).    

Per sample, combine the files demultiplexed by the forward adapter with those demultiplexed by the reverse adapter:     
`01_scripts/collect_samples.sh`
(note: currently expect warnings for combinations that weren't identified in the reverse complement adapter round, as there were much fewer reads remaining).     
(#todo: could this be done with a single adapter file with both forward and reverse adapters? Probably!)

#### c. Evaluate demultiplex results
Calculate the number of reads per sample, to produce `reads_per_sample2.txt`:     
`01_scripts/reads_per_sample.sh`
(note: contains code from: 'moving every second row to a new column with awk')

Then use the script `01_scripts/plot_demultiplex_result.R` that works on the read per sample table produced above. This will generate a horizontal barplot per sample 05_results/reads_per_sample.pdf    

(#todo: still may need to remove the reverse complement adapter, perhaps with a full cutadapt run).   

### 3B. De-multiplex with Porechop
currently using Adapters 1, 2, 3 and renaming each step (#todo, make a script to automate renaming)

#### i. De-multiplex using external barcodes
First we want to de-multiplex each concatenated read with Porechop, using the library ID barcodes.      
To do this, edit the adapter file to remove all of the adapters that you are not using as external library ID adapters, otherwise the program looks for those and starts to cut the read at these locations. This step is done later in the pipeline. (Here remove the first 94 PCR barcodes, leaving only 95 and 96 as library barcodes).      

If an internal barcode is seen at this step, the read is thrown out (#todo: fix this)     

To prepare for de-multiplex step one, rename adapters-1 to adapters.py
`mv ~/Programs/Porechop/porechop/adapters-1.py ~/Programs/Porechop/porechop/adapters.py`     

Demultiplex using external adapters:      
`~/Programs/Porechop/porechop-runner.py -i 03_basecalled/all_reads.fastq -b 03a_demultiplex_library -t 8 --adapter_threshold 90 --end_threshold 75`       
The results will be in the folder specified by -b. (here `03a_demultiplex_library`)      
This will produce a file per supplied external barcode, and a file containing reads without any adapter (none.fastq).       

Then rename the used adapter file back to it's original name:    
`mv ~/Programs/Porechop/porechop/adapters.py ~/Programs/Porechop/porechop/adapters-1.py`     

(#todo add 'extra-trim-end 0' to avoid cutting out cat adapter)

#### ii. De-concatenate using concatenation adapter within each library
To prepare for de-multiplex step two (de-concatenation only with concatenation adapter), rename adapters-2 to adapters.py
`mv ~/Programs/Porechop/porechop/adapters-2.py ~/Programs/Porechop/porechop/adapters.py`     

(#todo: automate for all barcode files)

`~/Programs/Porechop/porechop-runner.py -i 03a_demultiplex_library/BC96.fastq -o 03a_demultiplex_library/BC96_deconcat.fastq -t 16 --middle_threshold 75 --min_split_read_size 100 --extra_middle_trim_bad_side 0 --extra_middle_trim_good_side 0`     

Then rename the used adapter file back to it's original name:    
`mv ~/Programs/Porechop/porechop/adapters.py ~/Programs/Porechop/porechop/adapters-2.py`     

#### iii. De-multiplex the de-concatenated library
To prepare for de-multiplex step three (de-multiplex with sample ID barcodes), rename adapters-3 to adapters.py
`mv ~/Programs/Porechop/porechop/adapters-3.py ~/Programs/Porechop/porechop/adapters.py`     

Now run porechop to de-multiplex samples based on all 96 barcodes, search extra reads to identify all the barcodes.      
`~/Programs/Porechop/porechop-runner.py -i 03a_demultiplex_library/BC96_deconcat.fastq -b 03b_demultiplexed/ -t 16 --adapter_threshold 90 --end_threshold 75 --check_reads 100000`      


Then rename the used barcode file back to it's original name:     
`mv ~/Programs/Porechop/porechop/adapters.py ~/Programs/Porechop/porechop/adapters-3.py`     

If using subfolders, use "01c_cat_subsets" to compile all binned reads into combined folder (adjust to fit number of subdirectories).
`mkdir all_decat_cat`
`01_scripts/01c_cat_subsets.sh`

Move the demultiplexed samples into the appropriate folder.    
`cp -l 03b_demultiplexed/*.fastq 04_samples`      


### 4. Limiting genome to amplicons only (optional)
This is only necessary if you do not have a hotspot file already. It will pull the regions needed from a reference genome.    
Get the range in a bed file, then run the following to get an amplicon file of just the expected amplicons to align against
`GENOME="ch_WG00004_7.20170208.fasta"; bedtools getfasta -fi $GENOME -bed ch_WG00004_9.20170224.designed.bed -fo ch_WG00004_7.20170208_extracted.fa`


### 5. Align against reference regions
The demultiplexed fastq files are in `04_samples`.     

Index your reference genome using minimap2    
`minimap2 -d Otsh_subset.mmi path/to/the/genome_assembly.fasta`    

Set the genome variable and use minimap2 and samtools to align fastq samples in `04_samples` against the reference genome:    
`01_scripts/02_align.sh`      
Will produce stat files as well as indexed and sorted bam files.    

Run the following to generate alignments per chromosome in `05_results/sampleID_align_per_chr.txt`:    
`01_scripts/count_aligns_per_chr.sh`      

Then use the Rscript to generate figures:     
`01_scripts/plot_alignment_coverage.R`     
This will produce files `05_results/per_nucleotide_coverage.pdf` and `05_results/sampleID_align_per_chr.txt`, which requires the reads per sample table, the alignments per chromosome, as well as coverage statistics from the alignment.    

### 6. Call variants 

#### 6.a. Identify nucleotides at all loci
This will use the genome and all bam files in `04_mapped`.     
Set the genome variable to run pysamstats and output the nucleotides found at every position in all aligned sample files.    
`01_scripts/03_call_SNP_pysamstats.sh`       

#### 6.b. Prepare chromosome selection file 
This step will use a chromosome position file that specifies the chromosome and the position needed for the analysis.   
The chromosome position file will have the following format:     
`chrom \t pos \t ref \t var`       

Run the following script to prepare this chromosome position file into the formats needed for the analysis:     
`./01_scripts/prepare_locus_selection_file.sh`    

This will produce a file called `00_archive/locus_selection_w_index.txt` that will be used later. 

#### 6.c. Select loci of interest
Use the following script to add an index column to the data file and then select out the specific lines of interest from the data files.     
`./01_scripts/04_select_variants.sh`    

This will output files entitled `<sample>_var_w_index_selected.txt` to be used downstream.    

#### 7. Collect specific SNPs 
Use the Rscript interactively:    
`04_collect_select_genos.R`     

If you want to automatically run all samples, use the shell script, which will operate on all the file in `06_vcf/*"_var_w_index_selected.txt"`:      
`./01_scripts/05_run_collect_select_genos.sh`    
Note that this will also require that you have the guide file in `00_archive/locus_selection_w_index.txt`       

The results will end up in `07_rubias` in wide format, one file for each sample, with titles `<sample>_rubias_input.txt`      

Then run the following script to collect all input data and put into a single dataframe      
`01_scripts/06_merge_all_rubias_input_files.R`

#### 8. Rubias
