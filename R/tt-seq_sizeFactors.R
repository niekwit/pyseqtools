library(DESeq2)

#get parsed arguments
args <- commandArgs(trailingOnly=TRUE)
work.dir <- args[1]
script.dir <- args[2]
gtf <- args[3]
genome <- args[4]