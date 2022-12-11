library(tidyverse)
library(reshape2)

work.dir <- getwd()

#get HISAT2 log files
fileList <- Sys.glob(file.path(work.dir, "slurm", "slurm_hisat2*.log"))

#create df for storing alignment info
df <- data.frame(matrix(ncol = 4, nrow = length(fileList)))
names(df) <- c("sample","zero","one",">one")

#read data and add to df
for (i in 1:length(fileList)){
  file <- fileList[i]
  
  #get sample name
  sample <- system(paste("sed -n '1p'",file), intern = TRUE)
  
  #get alignment rates
  zero <- system(paste("grep 'aligned concordantly 0 times'",file), intern = TRUE)
  zero <- as.numeric(substr(strsplit(zero, " ")[[1]][6], start=2, stop=nchar(strsplit(zero, " ")[[1]][6])-2))
  
  one <- system(paste("grep 'aligned concordantly exactly 1 time'",file), intern = TRUE)
  one <- as.numeric(substr(strsplit(one, " ")[[1]][6], start=2, stop=nchar(strsplit(one, " ")[[1]][6])-2))
  
  more.one <- system(paste("grep 'aligned concordantly >1 times'",file), intern = TRUE)
  more.one <- as.numeric(substr(strsplit(more.one, " ")[[1]][6], start=2, stop=nchar(strsplit(more.one, " ")[[1]][6])-2))
  
  #add data to df
  df[i,1] <- sample
  df[i,2] <- zero
  df[i,3] <- one
  df[i,4] <- more.one
}
