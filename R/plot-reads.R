library(tidyverse)
library(reshape2)

#get parsed arguments
args <- commandArgs(trailingOnly = TRUE)
work.dir <- args[1]

#save read counts to df
df <- read.csv(file.path(work.dir,"bam","read-counts.csv"))
df <- df[order(df$sample), ] 

#create df for plotting
df.melt <- melt(df, value.name = sample)

#plot data
p <- ggplot(df.melt, aes(sample, Parent_Input)) +   
  geom_bar(aes(fill = variable), 
           position = "dodge", 
           stat = "identity",
           color = "black") +
  theme_bw() +
  theme(text=element_text(size=16),
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust=1)) +
  ylab("Read count") +
  xlab(NULL) +
  guides(fill=guide_legend(NULL)) +
  scale_fill_discrete(labels=c("pre-deduplication",
                               "post-deduplication"))

#save plot
ggsave(file.path(work.dir, "bam", "read-count_deduplication.pdf"), p,
       width = 7,
       height = 3.5)





