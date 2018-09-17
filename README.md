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


