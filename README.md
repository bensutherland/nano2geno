# nano2geno
Genotyping from nanopore data

Currently in **Development mode only**
Disclaimer: this is a simple pipeline that comes with no guarantees. At the moment it is purely for the author's use to better understand nanopore data.   

#### Requirements
Albacore basecaller     
Nanopolish
minimap2    
samtools    

#### Inputs
fast5     
genome.fa    
inner barcodes (.csv or .fa)    

### 1. Basecalling
If you have fast5 files to basecall, move the folder of fast5 files to `02_raw_data`. Note that fast5 files will be needed to index the fastq for SNP calling later.     
Note: for just getting stats on a fastq file without doing genotyping, collect fastq files and put in `03_basecalled/all_reads.fastq`, and skip to the next step. (e.g. `cp -r ../raw_data/PBT_trial*/*.fastq ./03_basecalled` then cat these files together into a fastq)         

Run `read_fast5_basecaller.py` recursively, with demultiplexing of nanopore barcodes on fast5 files within `02_raw_data` use the following:    
`./01_scripts/01_basecall.sh`     

Output will be in `03_basecalled/workspace`, with reads inside folders `pass/fail/calibration_strands`, and reads in the `Pass` folder in different barcode subfolders containing fastq files for each barcode.    

If you don't care what nanopore barcode was on your library, combine the fastq files into `03_basecalled/all_reads.fastq`     
`find ./03_basecalled/workspace/pass -name "*.fastq" | while read file ; do cat $file >> 03_basecalled/all_reads.fastq ; done`
(make sure to remove any previous version of all_reads.fastq or else this will combine all to this)   
(#todo: describe how to split this into sections to run in parallel)    

### 2. Data quality check
Run fastqc on this data:    
`fastqc 03_basecalled/all_reads.fastq -o 03_basecalled/fastqc -t 5`   
`multiqc -o 03_basecalled/fastqc/ 03_basecalled_fastqc`    

### 3. Demultiplexing the inner library
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

Then use the script `01_scripts/plot_demultiplex_result.R` that works on the read per sample table produced above. This will generate a horizontal barplot per sample `05_results/reads_per_sample.pdf`    

(#todo: still may need to remove the reverse complement adapter, perhaps with a full cutadapt run).   

### 4. Limiting genome to amplicons only
(#todo: not yet applied)
Get the range in a bed file, then run the following to get an amplicon file of just the expected amplicons to align against
`GENOME="ch_WG00004_7.20170208.fasta"; bedtools getfasta -fi $GENOME -bed ch_WG00004_9.20170224.designed.bed -fo ch_WG00004_7.20170208_extracted.fa`

### 5. Align against reference genome
The demultiplexed fastq files are in `04_samples`.     
Note: if want to only analyze a couple of files, move any unwanted into `04_samples/temp_storage`.   

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

### 6. Call SNPs
*this section still under development*
Use nanopolish index on the sample.fastq file using the sequencing_summary.txt file from basecalling and the folder with fast5 files. (#todo: is seq summary txt file necessary?)   
`nanopolish index -d /path/to/fast5/or/at/least/symlinks 04_samples/your_sample.fastq`    
All reads should be accounted for if this worked correctly. 

Then run nanopolish variants to make a vcf with SNPs:    
`nanopolish variants --progress -t 2 --reads 04_samples/your_sample.fastq --genome /path/to/genome.fa --ploidy 2 --bam 04_mapped/sample.bam -w chr:1-200 > 06_vcf/your_sample.vcf`    

Note: you will need to run nanopolish on each amplicon individually, so it will be run a total of the number of regions of interest. See details on applying xargs for this purpose here: https://github.com/jts/nanopolish/issues/224  

(#todo: make a script to automate the basecalling per region for all samples, then combine together)    
(#todo: make a file to be used for the -w flag in the above)    
