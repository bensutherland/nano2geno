# Call SNPS
# Simply call SNPs based on allelic ratio

# Clear space
# rm(list=ls())

# Set working directory
setwd("~/Documents/01_nanopore/nano2geno")

#### 01. Input Data ####
# guide file
selecting.dat <- read.delim2(file = "00_archive/Ots_amplicon_names_rubias_SNP_location_ref_var_match.txt"
                             , header = T)

colnames(selecting.dat)
dim(selecting.dat) #296 records

# empirical data file
vcf.dat <- read.delim2(file = "06_vcf/none_sel_var_match.txt"
                             , header = T)
colnames(vcf.dat)
vcf.dat <- vcf.dat[, c("chrom_pos", "chrom", "pos", "ref", "reads_all", "A", "C", "T", "G")]
colnames(vcf.dat)
dim(vcf.dat) #311 records

# Join the two files
data <- merge(x = vcf.dat, y = selecting.dat, by = "chrom_pos")
dim(data) # 291 records

colnames(data)
head(data)


#### 02. Collect and Format Data ####
# set nulls
genos <- NULL; arr.genos <- list() ; chr <- NULL

for(i in 2:nrow(data)){
  print(i)
  
  # find the two max alleles
  genos <- data[i, c(6:9)]
  print(genos)
  genos <- t(genos)
  print(genos)
  
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
  
  top.genos <- names(per.geno)[3:4]
  print(top.genos)
  
  top.genos.df <- rbind(top.genos.df, top.genos)
  
  # rename row
  rownames(top.genos.df)[i] <- chr
  
  # obtain the counts for each chr
  top.counts <- arr.genos[[chr]][3:4]
  print(top.counts)
  top.counts.df <- rbind(top.counts.df, top.counts)
  rownames(top.counts.df)[i] <- chr
  
  
}


top.genos.df
top.counts.df

#### Combine data (needs merge #todo)
all.data <- cbind(top.genos.df, top.counts.df)
all.data

colnames(all.data) <- c("alt", "ref", "r.alt", "r.ref")
all.data

all.data.df <- as.data.frame(as.matrix(all.data))
all.data.df

str(all.data.df)

# Coerce to numeric
all.data.df$r.alt <- as.numeric(as.character(all.data.df$r.alt))
all.data.df$r.ref <- as.numeric(as.character(all.data.df$r.ref))

all.data.df$alt <- as.character(all.data.df$alt)
all.data.df$ref <- as.character(all.data.df$ref)

str(all.data.df)

all.data.df


## Order as largest then second largest
all.data.df <- all.data.df[ , c(2, 4, 1, 3) ]
all.data.df


#### 03. Additional Stats ####
#### Calculate alleleic ratio (alt / ref)
all.data.df$al.ratio <- round(x = 
                                all.data.df[, "r.alt"] / (all.data.df[, "r.ref"] + all.data.df[, "r.alt"])
                              , digits = 2)

## Backup
# all.data.df.bck <- all.data.df


#### 04. Bin Genos ####
snp.call <- NULL; 

for(i in 1:nrow(all.data.df)){
  # 
  if(all.data.df$al.ratio[i] > 0.75){
    snp.call <- "homo.alt"
  } else if(all.data.df$al.ratio[i] < 0.25){
    snp.call <- "homo.ref"
  } else {
    snp.call <- "het"
  }
  print(snp.call)
  
  # add to df
  all.data.df$snp.call[i] <- snp.call
  
}
all.data.df

#### 05. Filter based on low read count ####
# Set minimum total read filter
threshold <- 20

for(i in 1:nrow(all.data.df)){
  if(all.data.df$r.ref[i] + all.data.df$r.alt[i] > threshold){
    all.data.df$keep[i] <- "TRUE"
  } else all.data.df$keep[i] <- "FALSE"
}

all.data.df

# add rownames as column
all.data.df$chr <- rownames(all.data.df)
all.data.df <- all.data.df[, c(8, 1:7)]
all.data.df

# Export results
write.csv(x = all.data.df, file = "06_vcf/genotypes.csv", row.names = F)

