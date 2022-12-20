library(tidyverse)
library(DiffBind)
library(rtracklayer)
library(ChIPseeker)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggupset)
library(ReactomePA)


#get parsed arguments
args <- commandArgs(trailingOnly = TRUE)
work.dir <- args[1]
genome <- args[2]

#load appropriate genome db
if (genome == "hg38"){library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  } else if (genome == "hg19"){library(TxDb.Hsapiens.UCSC.hg19.knownGene)}

#load sample info
sampleInfo <- read.csv(file.path(work.dir,"samples.csv"))

#get bam files
bamFiles <- Sys.glob(file.path(work.dir,"bam",genome,"*bam"))

#check if bam files are deduplicated
test <- str_detect(bamFiles,"dedupl")
dedup <- TRUE %in% test

if (dedup == TRUE){
  bamFiles <- Sys.glob(file.path(work.dir,"bam",genome,"*-dedupl.bam"))
}


###create diffbind sample sheet
ipSamples <- sort(unlist(sampleInfo[sampleInfo$type == "ip",]["sample"]))
inputSamples <- sort(unlist(sampleInfo[sampleInfo$type == "input",]["sample"]))

diffbind.sheet <- data.frame(matrix(ncol = 11, nrow = length(ipSamples)))
names(diffbind.sheet) <- c("SampleID","Tissue","Factor","Condition","Treatment","Replicate","bamReads",
                           "ControlID","bamControl","Peaks","PeakCaller")

#add ip sample names to diffbind.sheet
diffbind.sheet$SampleID <- ipSamples

treatment <- unlist(sampleInfo[sampleInfo$type == "ip",]["treatment"])
diffbind.sheet$Treatment <- treatment

factors <- unlist(sampleInfo[sampleInfo$type == "ip",]["factor"])
diffbind.sheet$Factor <- factors

#add conditions to diffbind.sheet
conditions <- unique(sampleInfo$genotype)
for (i in conditions){
  for (j in 1: nrow(diffbind.sheet)){
    if (str_detect(diffbind.sheet[j,1], i)){
      diffbind.sheet[j,4] <- i
    }
  }
}

#add replicate numbers to diffbind.sheet
diffbind.sheet <- diffbind.sheet %>% 
  group_by(Condition) %>% 
  mutate(Replicate = row_number())

#add peak caller to diffbind.sheet
diffbind.sheet$PeakCaller <- "macs"

#add ip sample bams to diffbind.sheet
for (i in ipSamples) {
  bam <- bamFiles[str_detect(bamFiles,i)]
  diffbind.sheet[diffbind.sheet$SampleID == i,7] <- bam
}

#add input sample names to diffbind.sheet
diffbind.sheet$ControlID <- inputSamples

#add input sample bam files to diffbind.sheet
for (i in inputSamples) {
  bam <- bamFiles[str_detect(bamFiles,i)]
  diffbind.sheet[diffbind.sheet$ControlID == i,9] <- bam
}

#add input sample bam files to diffbind.sheet
for (i in ipSamples){
  peak <- file.path(work.dir,"peak",genome,i,paste0(i,"_peaks.xls"))
  diffbind.sheet[diffbind.sheet$SampleID == i,10] <- peak
}

#write diffbind df to file
write.csv(diffbind.sheet,file=file.path(work.dir,"peak",genome,"diffbind_samples.csv"))






