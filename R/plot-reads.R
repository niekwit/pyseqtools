library(tidyverse)
library(reshape2)

#get parsed arguments
args <- commandArgs(trailingOnly = TRUE)
work.dir <- args[1]

df <- read.csv(file.path(work.dir,"bam","read-counts.csv"))