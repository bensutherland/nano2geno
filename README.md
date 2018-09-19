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
fast5 or fastq    
genome    
inner barcodes    

### 1. Basecalling
If you have fast5 files to basecall, move the folder of fast5 files to `02_raw_data`.    
If you have fastq files already, collect these and put in `03_basecalled/all_reads.fastq`, and skip to the next step. (e.g. `cp -r ../raw_data/PBT_trial*/*.fastq ./03_basecalled` then cat)         

Run `read_fast5_basecaller.py` using the following:    
`./01_scripts/01_basecall.sh`     

Output will be in `03_basecalled/workspace`   
Reads will be in subfolders depending on `pass/fail/calibration_strands` 
In the Pass folder will be different barcode subfolders containing called fastq files for that barcode.   

If you don't care what nanopore barcode was on your library, take the files and combine into one large fastq as `03_basecalled/all_reads.fastq`     
`find ./03_basecalled/workspace/pass -name "*.fastq" | while read file ; do cat $file >> 03_basecalled/all_reads.fastq ; done`
(make sure to remove any previous version of all_reads.fastq)   

### 2. Data quality check
Do some fastqc on this data:    
`fastqc`
`fastqc 03_basecalled/all_reads.fastq -o 03_basecalled/fastqc -t 5`   
`multiqc -o 03_basecalled/fastqc/ 03_basecalled_fastqc`    

### Demultiplexing the inner library
Create a fasta file with the adapters named
`grep -vE '^index' ./IonCode768.csv | awk -F"," '{ print ">"$2 "\n" $3 $5 "\n" }' -  | grep -vE '^$' - > IonCode768.fa`

Run cutadapt using this file:   
`cutadapt -a file:./../IonCode768.fa --untrimmed-output untrimmed.fastq -o 03b_demultiplexed/trimmed-{name}.fastq 03_basecalled/all_reads.fastq`

Then re-run it on the untrimmed.fastq with the reverse complement:    
`~/Documents/01_nanopore/nano2geno$ cutadapt -a file:./../IonCode768_revcomp.fa --untrimmed-output untrimmed_both.fastq -o 03b_demultiplexed/reverse_demultiplex/trimmed-{name}.fastq 03b_demultip
lexed/untrimmed.fastq`    

And finally join these two files:    
`01_scripts/collect_samples.sh`
(note: currently expect warnings for those that weren't found in the second round)

Audit the number of reads per sample, will produce `reads_per_sample2.txt`:
`01_scripts/reads_per_sample.sh`
(note: contains code from: 'moving every second row to a new column with awk')

(note this does not include a reverse complement adapters, but it should catch the barcode in the other side. May need to make a reverse complement adapter set to trim the adapters off)
To reverse complement, from Pierre Lindenbaum Biostars : https://www.biostars.org/p/189325/ 
` cat input.fa | while read L; do  echo $L; read L; echo "$L" | rev | tr "ATGC" "TACG" ; done`

`cat IonCode768.fa | while read L; do echo $L; read L; echo "$L" | rev | tr "ATGC" "TACG" ; done | sed -e "s/^M//" > IonCode768_revcomp.fa`


http://porecamp.github.io/2016/tutorials/mappingtute.html







## NEEDS CORRECTION


### Calling SNPs
#### 1. Data preprocessing
As described by @jts, one needs to index the output of the albacore basecaller:   
`nanopolish index -d ./02_raw_data/skip -s 03_basecalled/sequencing_summary.txt 03_basecalled/all_reads.fastq`   
Where the -d flag directs towards the fast5 files, and the -s flag points towards the output of albacore sequence summary, and the final output of the albacore basecaller. 
All reads should be accounted for if this worked correctly. 

Probably not the file you are looking for is entitled: `03_basecalled/all_reads.fastq.index.gzi`    





This requires that you align the fastq against your genome first to produce a bam file. 
minimap2 seems to be the go-to for nanopore data currently. 

Index:    
`minimap2 -d Otsh_subset.mmi the_genome_assembly.fasta`    

Align:
`minimap2 -ax map-ont ch_WG00004_7.20170208.fasta nano2geno/03_basecalled/all_reads.fastq > aln.sam`

Use this script to automate:    
`.01_scripts/02_align.sh`

Observe the stats on the alignment:   
`samtools stats 04_mapped/all_reads.sorted.bam > 04_mapped/all_reads.ali.stats.txt`

Collect information on coverage:   
`grep "^COV" 04_mapped/all_reads.ali.stats.txt > 04_mapped/all_reads.coverage.tx`


Remember, you can use 
`samtools flagstat aln.sam`

and to go from sam to bam:
`samtools view -Sb aln.sam > aln.bam`

Sort
`samtools sort aln_new.bam -o aln_new.sorted.bam`


Make sure to Nanopolish index your fasta file as well.  

Use Nanopolish to compute the consensus sequence
`python /home/ben/Programs/nanopolish/scripts/nanopolish_makerange.py`

Computes the consensus sequence of the genome assembly based on the nanopore reads: 
`python /home/ben/Programs/nanopolish/scripts/nanopolish_makerange.py ~/Documents/genomes/GCF_002021735.1_Okis_V1_genomic.fna | parallel --results nanopolish.results -P 8 nanopolish variants --consensus -o polished.{1}.vcf -w {1} -r 03_basecalled/all_reads.fa -b aln_new.sorted.bam -g ~/Documents/genomes/GCF_002021735.1_Okis_V1_genomic.fna -t 4 --min-candidate-frequency 0.1`



Following the porecamp.github.io tutorial:    
Variant calling with nanopolish:   
Three steps:    
1. align reads w/ aligner (done, this is `04_mapped/all_reads.sorted.bam`)
2. align events w/ nanopolish eventalign
3. call vcf with nanopolished variants


Nanopolish needs a fasta file, so use Fastq to fasta:
`fastq_to_fasta.py 03_basecalled/all_reads.fastq 03_basecalled/all_reads.fa`

And this also needs to be indexed by nanopolish:    
`nanopolish index -d ./02_raw_data/skip -s 03_basecalled/sequencing_summary.txt 03_basecalled/all_reads.fa`

Then align events w/ nanopolish eventalign:   
`nanopolish eventalign --reads 03_basecalled/all_reads.fa -b 04_mapped/all_reads.sorted.bam -g ~/Documents/genomes/GCF_002872995.1_Otsh_v1.0_genomic.fna --sam | samtools view -bS - | samtools sort -o 04_mapped/all_reads.eventalign.bam -`

Index the new file that nanopolish eventalign produced:    
`samtools index 04_mapped/all_reads.eventalign.bam`   

Then run the following, adding the window of interest:   
`nanopolish variants --progress -t 2 --reads 03_basecalled/all_reads.fa -o all_reads.eventalign.vcf -b 04_mapped/all_reads.sorted.bam -g ~/Documents/genomes/GCF_002872995.1_Otsh_v1.0_genomic.fna --snp --ploidy 4 -w NC_037115.1:39000000-3902000`

However, this produces a vcf without anything in it. I need to be able to find a section of the bam with SNPs actually present. Would I do this per sample? 
