library(tidyverse)
library(reshape2)

#get parsed arguments
args <- commandArgs(trailingOnly = TRUE)
work.dir <- args[1]

#get HISAT2 log files
fileList <- Sys.glob(file.path(work.dir, "slurm", "slurm_hisat2*[0-9].log"))

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

#order df
df <- df[order(df$sample), ] 

#create df for plotting
df.melt <- melt(df, value.name = "rate")

p <- ggplot(df.melt, aes(sample, rate)) +   
  geom_bar(aes(fill = variable), 
           position = "dodge", 
           stat = "identity",
           color = "black") +
  theme_bw() +
  theme(text=element_text(size=16),
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust=1)) +
  ylab("Alignment rate (%)") +
  xlab(NULL) +
  guides(fill=guide_legend("Aligned:")) +
  scale_fill_discrete(labels=c("0 times",
                                "1 time",
                                ">1 times"))

ggsave(file.path(work.dir, "bam", "alignment-rates.pdf"), p,
       width = 7,
       height = 3.5)





