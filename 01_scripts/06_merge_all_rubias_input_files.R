# Take the rubias input files from 07_rubias and make them into a single file

# Clear space
# rm(list=ls())

# Set working directory
setwd("~/Documents/01_nanopore/nano2geno")

#### 01. Input Data ####
data.files <- list.files(path = "07_rubias/", pattern = "_rubias_input.txt")
data.files

# Set nulls
rubias.all.df <- NULL; rubias.sample <- NULL
sample <- NULL

# Import each file and combine into a single dataframe
for(i in 1:length(data.files)){
  sample <- gsub(data.files[d], pattern = "\\_.*", replacement = "")
  print(sample)
  
  # Set filename
  filename <- paste0("07_rubias/", data.files[i])
  
  # import the data
  rubias.sample <- read.table(file = filename, sep = ",", header = TRUE)
  
  # Combine into the rest of the data
  rubias.all.df <- rbind(rubias.all.df, rubias.sample)
}

# str(rubias.all.df)

# Set format requirements
write.table(x = rubias.all.df, file = "07_rubias/rubias_all_input.txt", sep = ",", quote = F)
