# This script allows one to plot some information about per nucleotide or per chromosome coverage from alignments
# 

# Clear space
# rm(list=ls())
par(mfrow=c(1,3))
setwd("~/Documents/01_nanopore/nano2geno")

# Install Packages
# install.packages("ggplot2")
library(ggplot2)

#### Input reads per sample table ####
reads.per.sample <- read.delim2(file = "reads_per_sample2.txt", sep = " ", header = F
                                , col.names = c("sample", "num.reads"))
head(reads.per.sample)
dim(reads.per.sample)

# Remove untrimmed samples
reads.per.sample <- reads.per.sample[grep(pattern = "untrimmed", x = reads.per.sample$sample, ignore.case = T, perl = F, invert = T),]

# Clean up names
reads.per.sample$sample <- gsub(pattern = "04_samples/", replacement = "", x = reads.per.sample$sample)
# reads.per.sample$sample <- gsub(pattern = "04_samples/", replacement = "", x = reads.per.sample$sample)
reads.per.sample$sample <- gsub(pattern = ".fastq", replacement = "", x = reads.per.sample$sample)

head(reads.per.sample)

## Set up a loop
targets <- c("trimmed-IonCode_0135", "trimmed-IonCode_0151", "trimmed-IonCode_0366")
# 
# target.sample <- "trimmed-IonCode_0135"
# # target.sample <- "trimmed-IonCode_0151"
# # target.sample <- "trimmed-IonCode_0366"

for(i in 1:length(targets)){
  target.sample <- targets[i]

# Import coverage statistics per nucleotide
filename <- paste("04_mapped/", target.sample, ".cov.stats.txt", sep = "")
cov <- read.table(filename, sep = "\t")
dim(cov)

## ggplot option:
# ggplot(cov, aes(x=V3, y=V4)) + geom_bar(stat='identity') + xlab('coverage') + ylab('count')

# How many reads were present in this sample?
reads.in.target.sample <- reads.per.sample[which(reads.per.sample$sample==target.sample) ,2]

# I think what this is showing is the number of bases with the different depths of coverage... 
# It must be bases b/c it is far larger than the number of reads
plot(cov$V4 ~ cov$V3
     , xlim = c(0,500)
     , ylim = c(0,10000)
     , las = 1
     , xlab = "coverage"
     , ylab = "count")
text(x = 350, y = 8000, labels = paste("# reads = ", reads.in.target.sample))
# save out as 7 x 4, entitled 05_results/per_nucleotide_coverage.pdf
}

##### Plot alignment per chromosome ####
pdf(file = "05_results/alignments_per_chromosome.pdf", width = 11.5, height = 5)
par(mfrow = c(1,3), mar=c(10,4,3,3))
target.samples <- c("trimmed-IonCode_0135", "trimmed-IonCode_0151", "trimmed-IonCode_0366")

for(i in 1:length(target.samples)){
  target.sample <- target.samples[i]
  print(target.sample)
  # target.sample <- "trimmed-IonCode_0135"
# # target.sample <- "trimmed-IonCode_0151"
# # target.sample <- "trimmed-IonCode_0366"

# reporting
reads.in.target.sample <- reads.per.sample[which(reads.per.sample$sample==target.sample) ,2]

filename <- paste0("05_results/", target.sample, "_align_per_chr.txt")
align.per.chr <- read.table(file = filename, header = F, col.names = c("alignments","chr"))
head(align.per.chr)

# NOT YET FUNCTIONING
# remove the one chr with very low coverage (NC_002980.1): 
# align.per.chr <- align.per.chr[-c(which(align.per.chr$chr=="NC_002980.1")), ]
#grep(align.per.chr, pattern = "NC_002980.1", invert = T)

# sort (probably not necessary)
# align.per.chr[order(align.per.chr$chr),]

align.sum <- sum(align.per.chr$alignments)

barplot(height = align.per.chr$alignments, names.arg = align.per.chr$chr, las = 3
        , ylim = c(0,5000)
        , ylab = "# alignments")
text(x = 30, y = 4000, labels = paste("# reads = ", reads.in.target.sample))
text(x = 30, y = 3800, labels = paste("# alignments = ", align.sum))

}

dev.off()
