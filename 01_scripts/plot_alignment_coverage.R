setwd("~/Documents/01_nanopore/nano2geno")
library(ggplot2)

# 
cov <- read.table("04_mapped/all_reads.coverage.txt", sep = "\t")
head(cov)
cov[1,]
ggplot(cov, aes(x=V3, y=V4)) + geom_bar(stat='identity') + xlab('coverage') + ylab('count')



# I think what this is showing is the number of bases with the different depths of coverage... 
# It must be bases b/c it is far larger than the number of reads
plot(cov$V4 ~ cov$V3, xlim = c(0,100), ylim = c(0,10000)
     , las = 1
     , xlab = "coverage"
     , ylab = "count")
text(x = 80, y = 8000, labels = "cov, all indiv")
