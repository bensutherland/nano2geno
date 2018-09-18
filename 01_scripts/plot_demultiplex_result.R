setwd("~/Documents/01_nanopore/nano2geno")

reads.per.sample <- read.delim2(file = "reads_per_sample2.txt", sep = " ", header = F
                                , col.names = c("sample", "num.reads"))
head(reads.per.sample)
dim(reads.per.sample)

# remove 04_samples/
reads.per.sample$sample <- gsub(pattern = "04_samples/trimmed-IonCode_", replacement = "", x = reads.per.sample$sample)
reads.per.sample$sample <- gsub(pattern = "04_samples/", replacement = "", x = reads.per.sample$sample)


head(reads.per.sample)

# Plotting 
barplot(height = reads.per.sample$num.reads
        , names.arg = reads.per.sample$sample
        # , xaxt = "n"
        , las = 1
        , cex.names = 0.5
        , horiz = T
        , xlab = "num.reads / sample"
        , xlim = c(0,15000)
        , yaxt = "n"
        , main = paste(nrow(reads.per.sample), "samples")
        )
#axis(side = 1, at = seq(1:length(reads.per.sample$num.reads)), labels = reads.per.sample$sample, las = 3)

text(x = 8000, y = 25, labels = "Untrimmed reads")
