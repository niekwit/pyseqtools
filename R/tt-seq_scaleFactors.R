library(DESeq2)
library(dplyr)

work.dir <- getwd()

out.file <- file.path(work.dir,"scaleFactors.csv")

if (file.exists(out.file) == FALSE){
  #import sample table
  sampleTable <- read.csv(file.path(work.dir,"samples.csv"))
  sampleTable$condition <- factor(sampleTable$condition)
  sampleTable$genotype <- factor(sampleTable$genotype)
  
  #generate list of gene count file
  count.files <- Sys.glob(file.path(work.dir,"bam","R64-1-1","*","*ReadsPerGene.out.tab"))
  
  #check if number of gene count files matches sample number
  if(length(count.files) != length(rownames(sampleTable))){
    print("ERROR: number of ReadsPerGene.out.tab (STAR) files does not match sample number in samples.csv")
    exit()
  }
  
  #get gene/sample names for index names
  genes <- read.csv(count.files[1], skip=4, sep="\t", header=FALSE)
  genes <- genes$V1
  
  rownames(sampleTable) <- sampleTable$sample
  
  #generate count matrix for DESeq2
  countMatrix <- data.frame(matrix(ncol=1, nrow = 7127))
  col.names <- c("index")
  names(countMatrix) <- col.names
  countMatrix$index <- genes
  rownames(countMatrix) <- countMatrix$index
  
  
  for (i in count.files) {
    sample <- basename(i)
    sample <- sub("ReadsPerGene.out.tab","",sample)
    df <- read.csv(i, skip=4, sep="\t", header=FALSE)
    df <- subset(df, select=c(1,2))
    colnames(df) <- c("index", sample)
    countMatrix <- left_join(countMatrix, df, by="index")
  }
  #create named index
  rownames(countMatrix) <- countMatrix$index
  countMatrix$index <- NULL
  
  rownames(sampleTable) <- sampleTable$sample
  sampleTable$sample <- NULL
  
  #reorder data frame columns
  samples <- rownames(sampleTable)
  countMatrix <- countMatrix %>% select(samples)
  
  #check whether the order of the sample names between countMatrix and sampleTable is the same
  if(all(rownames(sampleTable) == colnames(countMatrix))){
    #create DESeq2 object
    if (length(unique(sampleTable$condition)) > 1){
      dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                                    colData = sampleTable,
                                    design = ~ condition)
    } else if (length(unique(sampleTable$condition)) == 1) {
      dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                                    colData = sampleTable,
                                    design = ~ genotype)
    }
    
    #calculate size factors
    sizeFactors <- estimateSizeFactors(dds)
    sizeFactors <- sizeFactors(sizeFactors)
    
    #calculate scale factors
    df <- as.data.frame(sizeFactors, row.names = NULL)
    df$sample <- row.names(df)
    df <- df[ , c(2,1)]
    df$scaleFactors <- 1 / df$sizeFactors
    df$sizeFactors <- NULL
    
    #write df to file
    write.csv(df, file = out.file, row.names = FALSE)
    
  } else{
    print("ERROR: row names of sampleTable and column names of countMatrix do not have the same order")
    exit()
  }
    
}
  

  