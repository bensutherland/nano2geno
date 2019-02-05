# Call SNPS
# Simply call SNPs based on allelic ratio

# Clear space
# rm(list=ls())

# Set working directory
setwd("~/Documents/01_nanopore/nano2geno")

# Set user variables
threshold <- 10 # the total number of reads required to call a genotype

#### 01. Input Data ####
# guide file
selecting.dat <- read.delim2(file = "00_archive/locus_selection_w_index.txt"
                             , header = T)
colnames(selecting.dat)
head(selecting.dat)
dim(selecting.dat) #296 records

# empirical data file
### TODO ### AUTOMATE MULTIPLE FILES
vcf.dat <- read.delim2(file = "06_vcf/none_var_w_index_selected.txt"
                             , header = T)
colnames(vcf.dat)
vcf.dat <- vcf.dat[, c("chrom_pos", "chrom", "pos", "ref", "reads_all", "A", "C", "T", "G")]
colnames(vcf.dat)
dim(vcf.dat) #311 records

# Join the two files
data <- merge(x = vcf.dat, y = selecting.dat, by = "chrom_pos", all.y = TRUE) # all.y will keep all records in guide file
dim(data) 

colnames(data)
head(data)


#### 02. Formatting data and sorting into major and minor allele ####
# Collect nucleotides per locus and sort in order of nucleotide count 
# set nulls
genos <- NULL; arr.genos <- list() ; chr <- NULL

for(i in 1:nrow(data)){
  print(i)
  
  # select the two max alleles
  genos <- data[i, c("A","C","T","G")]
  print(genos)
  genos <- t(genos) # transpose to sort by count
  print(genos)
  
  # save the chromosome name of interest
  chr <- data$chrom_pos[i]
  chr <- as.character(chr)
  print(chr)
  
  arr.genos[[ chr ]] <- genos[order(genos),] 
}

arr.genos


# Collect data out of list
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


#top.genos.df
head(top.genos.df)

#top.counts.df
head(top.counts.df)

# Combine these two dataframes by rowname
all.data <- merge(x = top.genos.df, y = top.counts.df, by ="row.names")
rownames(all.data) <- all.data[,"Row.names"]
all.data <- all.data[ , -(which(colnames(all.data)=="Row.names"))]

# Rename columns
colnames(all.data) <- c("min", "maj", "r.min", "r.maj")
head(all.data)

# Make data.frame
all.data.df <- as.data.frame(as.matrix(all.data))
head(all.data.df)
str(all.data.df)

## All values are factors, so...
# Coerce counts to numeric
all.data.df$r.min <- as.numeric(as.character(all.data.df$r.min))
all.data.df$r.maj <- as.numeric(as.character(all.data.df$r.maj))

# Coerce factor to character
all.data.df$min <- as.character(all.data.df$min)
all.data.df$maj <- as.character(all.data.df$maj)

str(all.data.df)

head(all.data.df)


## Order dataframe to show major then minor allele as largest then second largest
all.data.df <- all.data.df[ , c("maj", "r.maj", "min", "r.min") ]
head(all.data.df)


#### 03. Calculate allelic ratio (min / maj) ####
#### 
all.data.df$al.ratio <- round(x = 
                                all.data.df[, "r.min"] / (all.data.df[, "r.maj"] + all.data.df[, "r.min"])
                              , digits = 2)
head(all.data.df)


#### 04. Bin genotypes into homozygote or heterozygote ####
snp.call <- NULL; 

for(i in 1:nrow(all.data.df)){
  # 
  if(all.data.df$al.ratio[i] > 0.75){
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
  if(all.data.df$r.maj[i] + all.data.df$r.min[i] > threshold){
    all.data.df$keep[i] <- "TRUE"
  } else all.data.df$keep[i] <- "FALSE"
}

### TODO ### This throws an error due to NA values, need to fix that. 

all.data.df

# add rownames as column
all.data.df$chr <- rownames(all.data.df)

# Reorder to put chromosome name first
all.data.df <- all.data.df[, c("chr", "maj", "r.maj", "min", "r.min", "al.ratio", "snp.call", "keep")]
head(all.data.df)

# How many loci are removed (i.e. FALSE) # note: issue w/ NA
table(all.data.df$keep)


#### 06. Make genotypes column ####
genos1 <- NULL; genos2 <- NULL

for(i in 1:nrow(all.data.df)){
  if(all.data.df$snp.call[i]=="homo.maj"){
    #print("homo.maj")
    genos1[i] <- all.data.df$maj[i] # put major allele as genos1
    genos2[i] <- all.data.df$maj[i] # put major allele as genos2
    
  } else if(all.data.df$snp.call[i]=="het"){
    #print("het")
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


# Connect back to the guide file
all.data.w.guide.df <- merge(x = all.data.df, y = selecting.dat, by.x = "chr", by.y  = "chrom_pos")
dim(all.data.w.guide.df) # 296 records
head(all.data.w.guide.df)


##### Give true geno score w/ hotspot file
hotspot.call <- NULL

for(i in 1:nrow(all.data.w.guide.df)){
  if(all.data.w.guide.df$keep[i] == FALSE ){
    hotspot.call[i] <- NA
     
  } else if(all.data.w.guide.df$genos1[i]!=all.data.w.guide.df$genos2[i]){
    hotspot.call[i] <- "het"
    
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
all.data.w.guide.w.hotspot.call.df <- cbind(all.data.w.guide.df, hotspot.call)
head(all.data.w.guide.w.hotspot.call.df, n = 10)


##### 07. Convert to Rubias format ####
rubias1 <- NULL; rubias2 <- NULL
#scoring: homo.ref = 1,1 ; het = 1,2 ; homo.var = 2,2 

# Loop 
for(i in 1:nrow(all.data.w.guide.w.hotspot.call.df)){
  if(all.data.w.guide.w.hotspot.call.df$keep[i] == FALSE){
    rubias1[i] <- NA
    rubias2[i] <- NA ### MAY REMOVE THIS, NOT WORKING PROPERLY
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

### Hacky below ###
rubias.only.1 <- complete.df[, c("chrom", "rubias1")]
rubias.only.2 <- complete.df[, c("chrom", "rubias2")]

rownames(rubias.only.1) <- rubias.only.1$chrom
rownames(rubias.only.2) <- paste0(rubias.only.1$chrom, "_1")

head(rubias.only.1)
head(rubias.only.2)

t.x.1 <- t(rubias.only.1)
t.x.2 <- t(rubias.only.2)

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

### Missing: sample_type, collection, repunit, indiv
final_header <- c('sample_type', 'collection', 'repunit', 'indiv')
BC73_info <- c('mixed','lab',NA, 'BC73')
final_leader <- data.frame(final_header, BC73_info, stringsAsFactors=FALSE)
t_final_leader <- t(final_leader)
final_rubias_out.df <- cbind(t_final_leader, final.out.df)


### Data available in output.df.t
write.table(x = final_rubias_out.df, file = "rubias_inputBC73.csv", sep=",", row.names = F, col.names = F)




# Export results
write.csv(x = all.data.df, file = "06_vcf/genotypes.csv", row.names = F)
write.csv(x = all.data.w.guide.df, file = "temp.csv")









