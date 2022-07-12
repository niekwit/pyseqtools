library(DESeq2)

work.dir <- getwd()

#import sample table
sampleTable <- read.csv(file.path(work.dir,"deseq2.csv"))
sampleTable$file <- file.path(work.dir,"htseq-count", sampleTable$file)
sampleTable$condition <- factor(sampleTable$condition)
sampleTable$genotype <- factor(sampleTable$genotype)

#create DESeq2 object
directory <- file.path(work.dir,"htseq-count")
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition)

#generate size factors
sizeFactors <- estimateSizeFactors(ddsHTSeq)
sizeFactors <- sizeFactors(sizeFactors)

#prepare df to write to file
df <- as.data.frame(sizeFactors, row.names = NULL)
df$sample <- row.names(df)
df <- df[ , c(2,1)]

#write df to file
write.csv(df, file = file.path(work.dir,"sizeFactors.csv,"), row.names = FALSE)


