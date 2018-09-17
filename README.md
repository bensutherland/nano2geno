# nano2geno
Genotyping from nanopore data

Disclaimer: this is a simple pipeline that comes with no guarantees. At the moment it is purely for the author's use to better understand nanopore data.   

#### Requirements
Albacore basecaller
Nanopolish

### Basecalling
If doing local basecalling, move the folder containing the fast5 files to be basecalled to `02_raw_data`, and run the script `01_scripts/01_basecall.sh`. This script will call the python script `read_fast5_basecaller.py`.     
Your output will be put into the folder `03_basecalled`, specifically into the subfolder `workspace`, which will contain reads labeled as `pass/fail/calibration_strands`. Within the Pass folder will be different barcode subfolders, and within each of these will be the called fastq files.      

If you are going to do separate analyses of different libraries, these will matter, but if not, then you can just take all of the files and combine them into one big fastq file.

This will export all of the files as `03_basecalled/all_reads.fastq`     
`find ./03_basecalled/workspace/pass -name "*.fastq" | while read file ; do cat $file >> 03_basecalled/all_reads.fastq ; done`

### De-multiplexing the inner library
To add

### Calling SNPs
#### 1. Data preprocessing
As described by @jts, one needs to index the output of the albacore basecaller:   
`nanopolish index -d ./02_raw_data/skip -s 03_basecalled/sequencing_summary.txt 03_basecalled/all_reads.fastq`   
Where the -d flag directs towards the fast5 files, and the -s flag points towards the output of albacore sequence summary, and the final output of the albacore basecaller. 
All reads should be accounted for if this worked correctly. 

Importantly! the fasta file you are looking for is entitled: `03_basecalled/all_reads.fastq.index.gzi`    





This requires that you align the fastq against your genome first to produce a bam file. 
minimap2 seems to be the go-to for nanopore data currently. 

Index:    
`minimap2 -d Otsh_subset.mmi ch_WG00004_7.20170208.fasta`    

Align:
`minimap2 -ax map-ont ch_WG00004_7.20170208.fasta nano2geno/03_basecalled/all_reads.fastq > aln.sam`

Remember, you can use 
samtools flagstat aln.sam 

and to go from sam to bam:
`samtools view -Sb aln.sam > aln.bam`


Further following jts:    
`bwa index ch_WG00004_7.20170208.fasta`
Align basecalled reads to the draft genome:    
``


