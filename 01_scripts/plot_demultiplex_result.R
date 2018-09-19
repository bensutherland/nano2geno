# Plotting the number of reads per sample
# requires input of 'reads_per_sample2.txt'

# Clear space
# rm(list=ls())

setwd("~/Documents/01_nanopore/nano2geno")

reads.per.sample <- read.delim2(file = "reads_per_sample2.txt", sep = " ", header = F
                                , col.names = c("sample", "num.reads"))
head(reads.per.sample)
dim(reads.per.sample)

# Remove untrimmed samples
reads.per.sample <- reads.per.sample[grep(pattern = "untrimmed", x = reads.per.sample$sample, ignore.case = T, perl = F, invert = T),]

# Clean up names
reads.per.sample$sample <- gsub(pattern = "04_samples/trimmed-IonCode_", replacement = "", x = reads.per.sample$sample)
reads.per.sample$sample <- gsub(pattern = "04_samples/", replacement = "", x = reads.per.sample$sample)
reads.per.sample$sample <- gsub(pattern = ".fastq", replacement = "", x = reads.per.sample$sample)


head(reads.per.sample)


per.sample.min.depth <- 447 * 40


# Plotting
pdf(file = "05_results/reads_per_sample.pdf", width = 6, height = 4 )
options(scipen = 99999999)
par(mar=c(4,4,3,3))
barplot(height = reads.per.sample$num.reads
        , names.arg = reads.per.sample$sample
        # , xaxt = "n"
        #, las = 3
        , cex.names = 0.5
        , horiz = T
        , xlab = "reads/sample"
        , xlim = c(0,600000)
        , yaxt = "n"
        , ylab = "samples"
        )
#axis(side = 1, at = seq(1:length(reads.per.sample$num.reads)), labels = reads.per.sample$sample, las = 3)

# text(x = 8000, y = 25, labels = "Untrimmed reads")
text(x = 450000, y = 700, labels = paste("# samples =", nrow(reads.per.sample)))
text(x = 450000, y = 650, labels = paste("med. = ", 
      round( median(x = reads.per.sample$num.reads))) )
text(x = 450000, y = 600, labels = paste("sd. = ", 
      round( sd(x = reads.per.sample$num.reads))) )
text(x = 450000, y = 550, labels = paste("min. = ", 
      round( min(x = reads.per.sample$num.reads))) )
text(x = 450000, y = 500, labels = paste("min. = ", 
     round( max(x = reads.per.sample$num.reads))) )


abline(v = per.sample.min.depth, lty = 2, col = "red")
dev.off()
