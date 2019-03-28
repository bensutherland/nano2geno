### Nano_GSI -> SNP analysis R script
## 1: Sequence Library on minION
## 2: Generate fastq data
## 3: Quality control with fastqc or multiqc
## 4: Demultiplex
## 5: SNP calling with marginAlign
## 6: transform marginAlign DB to match rubias input format
## 7: GSI using rubias

# create new directory for the sequencing run
setwd("/home/cdeeg/Documents/01_nanopore/nano2geno/")

#confirm dir
getwd()

# load the packages
library(rubias)
library(tidyverse)
library(ggplot2)

#### make the csv into a character vector!!!!!!!

## load the database -> adjust for Oki or Ots ##

reference_DB <- read.table("00_archive/Oki_rubias_baseline.txt",  sep="\t", header=TRUE, stringsAsFactors=FALSE)
reference_DB <-reference_DB %>% mutate_all(as.character)

## If going from scractch: ##
#mixed_stock <- read.csv("xxx.csv", header=TRUE, stringsAsFactors=FALSE)
#mixed_stock <- mixed_stock %>% mutate_all(as.character)
#loci <- colnames(reference_DB)
#mixed_stock <- mixed_stock[,loci]
#mixed_stock <- mixed_stock %>% mutate_all(as.character)

## If continuing from nano2geno after '01_scripts/06_merge_all_rubias_input_files.R`:
mixed_stock <- rubias.all.df
mixed_stock <- mixed_stock %>% mutate_all(as.character)
loci <- colnames(reference_DB)
mixed_stock <- mixed_stock[,loci]
mixed_stock <- mixed_stock %>% mutate_all(as.character)




## analyse the genetic mixture of the stock ##

mix_est <- infer_mixture(reference = reference_DB, mixture = mixed_stock, gen_start_col = 5)

## calculate mixing proportions; adjust output for Oki or Ots ##
rep_mix_ests <- mix_est$mixing_proportions %>%
  group_by(mixture_collection, repunit) %>%
  summarise(repprop = sum(pi))  
# write.csv(rep_mix_ests,"Oki_rep_mixture_est.csv")

## calculate individuals posteriors; adjust output for Oki or Ots ##
rep_indiv_ests <- mix_est$indiv_posteriors %>%
  group_by(mixture_collection, indiv, repunit) %>%
  summarise(rep_pofz = sum(PofZ))
# write.csv(rep_indiv_ests,"Oki_rep_individual_est.csv")

## find the top 5 most abundant: ##
top5 <- rep_mix_ests %>%
  filter(mixture_collection == "lab") %>% 
  arrange(desc(repprop)) %>%
  slice(1:5)

## check how many MCMC sweeps were done: ##
nsweeps <- max(mix_est$mix_prop_traces$sweep)

## keep only rec1, then discard the first 200 sweeps as burn-in, ##
## and then aggregate over reporting units ##
## and then keep only the top5 from above ##
trace_subset <- mix_est$mix_prop_traces %>%
filter(mixture_collection == "lab", sweep > 200) %>%
group_by(sweep, repunit) %>%
summarise(repprop = sum(pi)) %>% 
filter(repunit %in% top5$repunit)

## make pie chart of comunity structure composition ##
#ggplot(trace_subset, aes(x = "", y = repprop, colour = repunit)) + geom_col()
#barp <- ggplot(trace_subset, aes(x = "", y = repprop, colour = repunit)) + geom_col()
#pie <- barp + coord_polar("y", start=0)
#plot(pie)

## Check rep geom dens ##
ggplot(trace_subset, aes(x = repprop, colour = repunit)) + geom_density()

