options(java.parameters = "-Xmx30000m") #increase available memory for java

library(DESeq2)
library(dplyr)
library(RColorBrewer)
library(pheatmap)
library(biomaRt)

library(xlsx)



#get parsed arguments
args <- commandArgs(trailingOnly=TRUE)
genome <- args[1]
cutoff <- args[2]

#get working directory
work.dir <- getwd()

#create DESeq2 output file name
deseq2.output <- file.path("work.dir","deseq2","DESeq2.csv")

if (file.exists(deseq2.output) == FALSE) {
  #import sample table
  sampleTable <- read.csv(file.path(work.dir,"samples.csv"))
  
  #set conditions for DESeq2
  sampleTable$geno_cond <- paste(sampleTable$genotype,sampleTable$condition, sep="_")
  
  #extract reference conditions
  sample_references <- sampleTable[sampleTable$ref %in% "ref",]
  sample_references <- unique(sample_references$geno_cond)
  sampleTable <- sampleTable %>% select(-one_of("genotype","condition","ref"))
  
  sampleTable$geno_cond <- factor(sampleTable$geno_cond)
  rownames(sampleTable) <- sampleTable$sample
  sampleTable$sample <- NULL
  
  #generate list of gene count files
  count.files <- Sys.glob(file.path(work.dir,"bam",genome,"*","*ReadsPerGene.out.tab"))
  
  #check if number of gene count files matches sample number
  if(length(count.files) != length(rownames(sampleTable))){
    print("ERROR: number of ReadsPerGene.out.tab (STAR) files does not match sample number in samples.csv")
    exit()
  }
  
  #get gene/sample names for index names
  genes <- read.csv(count.files[1], skip=4, sep="\t", header=FALSE)
  genes <- genes$V1
  
  #generate count matrix for DESeq2
  countMatrix <- data.frame(matrix(ncol=1, nrow = length(genes)))
  col.names <- c("index")
  names(countMatrix) <- col.names
  countMatrix$index <- genes
  rownames(countMatrix) <- countMatrix$index
  
  #add count data to countMatrix
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
  
  #reorder data frame columns
  samples <- rownames(sampleTable)
  countMatrix <- countMatrix %>% dplyr::select(samples)
  
  #load data for gene annotation
  mart <- useMart("ensembl")
  mart <- useDataset("hsapiens_gene_ensembl", mart = mart)
  
  #check whether the order of the sample names between countMatrix and sampleTable is the same
  if(all(rownames(sampleTable) == colnames(countMatrix))){
    #create DESeq2 object
    dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                                  colData = sampleTable,
                                  design = ~ geno_cond)
    
    #create list to store all data frames
    df.master <- list()
    
    #perform DESeq2 for each reference (compare to all other samples)
    for (reference in sample_references){
      #set reference level for DESeq2
      dds$geno_cond <- relevel(dds$geno_cond, reference)
      
      #apply size factors from yeast spike-in
      scale.factors <- read.csv(file.path(work.dir, "scaleFactors.csv"))
      size.factors <- scale.factors
      size.factors$sizeFactors <- 1 / size.factors$scaleFactors  
      size.factors$scaleFactors <- NULL
      sizeFactors(dds) <- size.factors$sizeFactors
      
      #run DESeq2
      dds <- DESeq(dds) 
      
      #get comparisons
      results.names <- resultsNames(dds)
      results.names <- strsplit(results.names," ")
      results.names[1] <- NULL
      
      #get results for each comparison and add to df.master
      for (results in results.names){
        #get results for comparison
        res <- results(dds, name=results)
        comparison <- sub("geno_cond_","",results)
        
        #create data frame with results
        df <- as.data.frame(res@listData)
        
        #annotate data frame
        df$ensembl_id <- res@rownames
        
        
        genes.table <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "external_gene_name","description",
                                                                       "gene_biotype", "chromosome_name","start_position",
                                                                       "end_position", "percentage_gene_gc_content"), 
                             values = df$ensembl_id, mart = mart) 
        names(genes.table)[1] <- "ensembl_id"
        
        df <- left_join(df,genes.table,by="ensembl_id")
        
        #add data to master data frame
        df$gene_length <- df$end_position - df$start_position
        df.master[[length(df.master)+1]] <- df
        names(df.master)[length(df.master)] <- comparison
        
        
      }
    }
    
    #save data frames to separate sheets of an excel file
    dir.create(file.path(work.dir,"DESeq2"), showWarnings = FALSE)
    excel <- file.path(work.dir,"DESeq2","DESeq2.xlsx")
    for (i in 1:length(df.master)){
      sheet <- names(df.master)[i]
      if (i == 1) {
        write.xlsx2(df.master[[i]], file=excel, sheetName=sheet, row.names = FALSE)
      }else{write.xlsx2(df.master[[i]], file=excel, sheetName=sheet, append = TRUE, row.names = FALSE)}
      
    }
  }
}
  





