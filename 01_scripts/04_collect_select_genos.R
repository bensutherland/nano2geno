#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE) # enables use of arguments from command line

# Call SNPS
# Simply call SNPs based on allelic ratio

# Clear space
# rm(list=ls())

# Set working directory
setwd("~/Documents/01_nanopore/nano2geno")

# Set user variables
threshold <- 10 # the total number of reads required to call a genotype

#### Run interactive or automated ####
# sample <- "none" # interactive mode
# sample <- "BC10" # debugging, interactive mode
sample <- args[1] # automated mode

#### 01. Input Data ####
# Guide file
selecting.dat <- read.delim2(file = "00_archive/locus_selection_w_index.txt"
                             , header = T)
colnames(selecting.dat)
head(selecting.dat)
dim(selecting.dat) #296 records

# Empirical data file

### LOOP SECTION 
# ### TODO ### AUTOMATE MULTIPLE FILES
# rubias.all.df <- NULL; rubias.all <- NULL
# 
# #### START LOOP #####
# data.files <- list.files(path = "06_vcf", pattern = "var_w_index_selected.txt")
# data.files

# sample <- NULL
# for(d in 1:length(data.files)){
#   sample <- gsub(data.files[d], pattern = "\\_.*", replacement = "")
#   print(sample)

### END LOOP SECTION #####


### Import sample file
input.filename <- paste0("06_vcf/", sample, "_var_w_index_selected.txt")
print(input.filename)

vcf.dat <- read.delim2(file = input.filename, header = T)
colnames(vcf.dat)
vcf.dat <- vcf.dat[, c("chrom_pos", "chrom", "pos", "ref", "reads_all", "A", "C", "T", "G")]
colnames(vcf.dat)
dim(vcf.dat) #311 records

# Join the two files
data <- merge(x = vcf.dat, y = selecting.dat, by = "chrom_pos", all.y = F) # all.y will keep all records in guide file, we don't want that here
dim(data) 

colnames(data)
head(data)


#### 02. Identifying the major and minor allele ####
# Solely based on the two nucleotides with the most reads at a site (max alleles)
# Collect nucleotides per locus and sort in order of nucleotide count 

# Identify two max alleles
genos <- NULL; arr.genos <- list() ; chr <- NULL

for(i in 1:nrow(data)){
  print(i)
  
  # Per row, pull out the nucleotide coverage
  genos <- data[i, c("A","C","T","G")]
  print(genos)
  
  # Transpose to sort by count
  genos <- t(genos) 
  print(genos)
  
  # Collect the chromosome name of interest
  chr <- data$chrom_pos[i]
  chr <- as.character(chr)
  print(chr)
  
  # Save ordered nucleotides
  arr.genos[[ chr ]] <- genos[order(genos),] 
}

arr.genos


# Collect data from list into dataframe
top.genos.df <- NULL; chr <- NULL; temp <- NULL; top.genos <- NULL
top.counts <- NULL; top.counts.df <- NULL

for(i in 1:length(arr.genos)){
  
  # define chr
  chr <- names(arr.genos)[i] 
  print(chr)
  
  # obtain the genos for each chr
  per.geno  <- arr.genos[[chr]]
  print(per.geno)
  
  # Identify top two nucleotides by count
  top.genos <- names(per.geno)[3:4]
  print(top.genos)
  
  # Collect data and add to dataframe
  top.genos.df <- rbind(top.genos.df, top.genos)
  
  # rename row with chromosome of interest
  rownames(top.genos.df)[i] <- chr
  
  # obtain the counts for each chr (read depth)
  top.counts <- arr.genos[[chr]][3:4]
  print(top.counts)
  top.counts.df <- rbind(top.counts.df, top.counts)
  colnames(top.counts.df) <- c("1","2")
  rownames(top.counts.df)[i] <- chr
  
}


head(top.genos.df)
head(top.counts.df)

# Combine these two dataframes by rowname
all.data <- merge(x = top.genos.df, y = top.counts.df, by ="row.names")
rownames(all.data) <- all.data[,"Row.names"]
all.data <- all.data[ , -(which(colnames(all.data)=="Row.names"))]

# Rename columns
colnames(all.data) <- c("min", "maj", "r.min", "r.maj")
head(all.data)

# Convert to dataframe
all.data.df <- as.data.frame(as.matrix(all.data), stringsAsFactors = F)
head(all.data.df)
str(all.data.df)

# Counts as numeric
all.data.df$r.min <- as.numeric(all.data.df$r.min)
all.data.df$r.maj <- as.numeric(all.data.df$r.maj)

## Re-order columns
all.data.df <- all.data.df[ , c("maj", "r.maj", "min", "r.min") ]
head(all.data.df)


#### 03. Calculate allelic ratio (min / maj) ####
#### 
all.data.df$al.ratio <- round(x = 
                                all.data.df[, "r.min"] / (all.data.df[, "r.maj"] + all.data.df[, "r.min"])
                              , digits = 2)
head(all.data.df)


#### 04. Bin genotypes into homozygote or heterozygote ####
# Based on allelic ratio, identify heterozygous or homozygous calls
snp.call <- NULL; 

for(i in 1:nrow(all.data.df)){
  
  # Add option for if NA
  if(is.na(all.data.df$al.ratio[i])==T){
    snp.call <- NA
  } else if(all.data.df$al.ratio[i] > 0.75){
    snp.call <- "homo.min"
  } else if(all.data.df$al.ratio[i] < 0.25){
    snp.call <- "homo.maj"
  } else {
    snp.call <- "het"
  }
  print(snp.call)
  
  # add to df
  all.data.df$snp.call[i] <- snp.call
  
}

head(all.data.df, n = 10)
# Note: this will never produce a homozygote minor, due to the ordering of alleles by read count above


#### 05. Filter based on low read count ####
# Set minimum total read filter
print(paste0("your chosen threshold is = ", threshold))

# Set true/false whether locus is kept or removed based on threshold
for(i in 1:nrow(all.data.df)){
  
  # If reads are greater than threshold
  if(all.data.df$r.maj[i] + all.data.df$r.min[i] > threshold){
    all.data.df$keep[i] <- "TRUE"
    
  } else all.data.df$keep[i] <- "FALSE"
}

all.data.df

# add rownames as column
all.data.df$chr <- rownames(all.data.df)

# Reorder to put chromosome name first
all.data.df <- all.data.df[, c("chr", "maj", "r.maj", "min", "r.min", "al.ratio", "snp.call", "keep")]
head(all.data.df)

# How many loci are removed (i.e. FALSE) # note: issue w/ NA
table(all.data.df$keep)


#### 06. Make genotypes column ####
# Make two vectors, genos1 and genos2 that contain the allele based on the allelic ratio call above
genos1 <- NULL; genos2 <- NULL

for(i in 1:nrow(all.data.df)){
  # Add option for missing data
  if(is.na(all.data.df$snp.call[i])==TRUE){
    genos1[i] <- NA # give NAs to genos
    genos2[i] <- NA
    
  } else if(all.data.df$snp.call[i]=="homo.maj"){
    genos1[i] <- all.data.df$maj[i] # put major allele as genos1
    genos2[i] <- all.data.df$maj[i] # put major allele as genos2
    
  } else if(all.data.df$snp.call[i]=="het"){
    genos1[i] <- all.data.df$maj[i] # put major allele as genos1
    genos2[i] <- all.data.df$min[i] # put minor allele as genos2
  }
}


genos1
genos2

# Add these calls to df
all.data.df$genos1 <- genos1
all.data.df$genos2 <- genos2

head(all.data.df)

##### 07. Identify actual genotypes ####
# Connect to the guide file
all.data.w.guide.df <- merge(x = all.data.df, y = selecting.dat, by.x = "chr", by.y  = "chrom_pos")
dim(all.data.w.guide.df) # 296 records
head(all.data.w.guide.df)


##### Give true geno score w/ hotspot file
# Use the exact nucleotide to identify the ref or var SNP for homozygous calls
hotspot.call <- NULL

for(i in 1:nrow(all.data.w.guide.df)){
  
  # Make filtered out markers NA
  if(all.data.w.guide.df$keep[i] == FALSE ){
    hotspot.call[i] <- NA
  
    # If the two genos aren't the same, it's a het   
  } else if(all.data.w.guide.df$genos1[i]!=all.data.w.guide.df$genos2[i]){
    hotspot.call[i] <- "het"
    
    # If two genos are the same, check if homozygous ref or var depending on guide file
  } else if( (all.data.w.guide.df$genos1[i]==all.data.w.guide.df$genos2[i]) # is a homozygote 
             && all.data.w.guide.df$genos1[i]==all.data.w.guide.df$ref[i]){ # is a homo.ref
    hotspot.call[i] <- "homo.ref"
    
  } else if( (all.data.w.guide.df$genos1[i]==all.data.w.guide.df$genos2[i]) # is a homozygote 
             && all.data.w.guide.df$genos1[i]==all.data.w.guide.df$var[i]){ # is a homo.var
    hotspot.call[i] <- "homo.var"
  }
}

# TODO: confirm this het is the same as the expected het, otherwise throw NA #

hotspot.call

table(hotspot.call) # shows the number of markers with each call

table(is.na(hotspot.call)) # shows the number of dropped due to low coverage

all.data.w.guide.w.hotspot.call.df <- cbind(all.data.w.guide.df, hotspot.call)
head(all.data.w.guide.w.hotspot.call.df, n = 10)


#### 08. Add dummy rows for markers in guide file not in vcf ####
dim(all.data.w.guide.w.hotspot.call.df)
dim(selecting.dat)
data.w.dummy.rows <- merge(x = all.data.w.guide.w.hotspot.call.df, y = selecting.dat
                           , by.x = "chr", by.y = "chrom_pos"
                           , all.y = T) # all.y will keep all records in guide file
dim(data.w.dummy.rows) 
tail(data.w.dummy.rows)

# Rename for downstream
all.data.w.guide.w.hotspot.call.df <- data.w.dummy.rows

##### 09. Convert to Rubias format ####
# Provide scoring necessary for rubias input
rubias1 <- NULL; rubias2 <- NULL
#scoring: homo.ref = 1,1 ; het = 1,2 ; homo.var = 2,2 

# Convert to rubias format based on hotspot.call
for(i in 1:nrow(all.data.w.guide.w.hotspot.call.df)){
  
  # Any locus that was removed by filter or not found at all
  if(is.na(all.data.w.guide.w.hotspot.call.df$hotspot.call[i]) == TRUE){
    rubias1[i] <- NA
    rubias2[i] <- NA
    
    # otherwise use the hotspot call to code in rubias format
  } else if(all.data.w.guide.w.hotspot.call.df$hotspot.call[i]=="homo.ref"){
    rubias1[i] <- 1
    rubias2[i] <- 1  
  } else if(all.data.w.guide.w.hotspot.call.df$hotspot.call[i]=="het"){
    rubias1[i] <- 1
    rubias2[i] <- 2
  } else if(all.data.w.guide.w.hotspot.call.df$hotspot.call[i]=="homo.var"){
    rubias1[i] <- 2
    rubias2[i] <- 2
  }
}

complete.df <- cbind(all.data.w.guide.w.hotspot.call.df, rubias1, rubias2)

head(complete.df)


#### 10. Convert long form to wide form ####
### Make long form into horizontal form ###
## TODO ## Fix this section up
rubias.only.1 <- complete.df[, c("chrom.y", "rubias1")]
rubias.only.2 <- complete.df[, c("chrom.y", "rubias2")]

rownames(rubias.only.1) <- rubias.only.1$chrom
rownames(rubias.only.2) <- paste0(rubias.only.1$chrom, "_1")

head(rubias.only.1)
head(rubias.only.2)

# Transpose
t.x.1 <- t(rubias.only.1)
t.x.2 <- t(rubias.only.2)

# View data
t.x.1[1:2, 1:5]
t.x.2[1:2, 1:5]

test <- cbind(t.x.1, t.x.2)
dim(test)

test[1:2,1:6]
test[1:2,570:582]

test.df <- as.data.frame(test)
test.df[1:2,1:6]
test.df[1:2,570:582]

output.df <- test.df["rubias1" ,]
output.df
output.df.t <- t(output.df)

output.df.t

output.df.t <- output.df.t[ order(row.names(output.df.t)), ]
head(output.df.t)

str(output.df.t)

names(output.df.t)
final.out.df <- rbind(names(output.df.t), output.df.t)
dim(final.out.df)
final.out.df[1:2,1:5]


# # Put into alphabetic order
# test <- t(final.out.df)
# colnames(test) <- c("chr", "geno")
# test[order(test$chr)
# final.out.df



#### 11. Add metadata to rubias output ####
## Add header info: sample_type, collection, repunit, indiv
final.header <- NULL; sample.info <- NULL
final.header <- c("sample_type", "collection", "repunit", "indiv")
sample.info <- c("mixed", "lab", NA, sample) ## Change sample.x to be the rep of the loop

final.leader <- data.frame(final.header, sample.info, stringsAsFactors=FALSE)
t.final.leader <- t(final.leader)

# Add metadata section to the front of the output dataframe
rubias.out.df <- cbind(t.final.leader, final.out.df)

rubias.out.df[1:2,1:10]

# move the first row into the header
colnames(rubias.out.df) <- rubias.out.df[1,]
rubias.out.df[1:2,1:10]

out.filename <- paste0("07_rubias/", sample, "_rubias_out.txt")
write.table(x = rubias.out.df, file = out.filename, sep = ","
            , col.names = F)


### BOTTOM OF LOOP SECTION
# }
### END BOTTOM OF LOOP SECTION

### TODO ### Next need to bring back all of these rubias files and join them together
# test <- read.delim2(file = out.filename, header = T, sep = ",", stringsAsFactors = F)
# str(test)

### FRAGMENTS ###
# only keep the sample row and the header
# rubias.out.df <- rubias.out.df["sample.info",]
# str(as.data.frame(rubias.out.df))
# names(rubias.out.df[1:10])
# rubias.combine.fraction.1 <- t(rubias.out.df)
# head(rubias.combine.fraction.1)
# merge(x = rubias.combine.fraction.1, rubias.all.df)


##### LOOP OPTION
##### THEN BIND THIS ON TO THE PRECEDING RUN #####

# write.csv2(x = rubias.out.df, file = "")
# 
# write.table(x = rubias.out.df, file = "rubias_inputBC73.csv", sep=",", row.names = F, col.names = F)
# 
# # Export results
# write.csv(x = all.data.df, file = "06_vcf/genotypes.csv", row.names = F)
# write.csv(x = all.data.w.guide.df, file = "temp.csv")
