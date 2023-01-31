#check if required packages are installed, if not install them
packages <- rownames(installed.packages())
biocmanager.packages <- c("tximport","DESeq2",
                          "GenomicFeatures","EnsDb.Mmusculus.v79",
                          "EnsDb.Hsapiens.v79","apeglm")
cran.packages <- c("BiocManager","ggplot2",
                   "readr","dplyr")

cran.packages2install <- cran.packages[! cran.packages %in% packages]
biocmanager.packages2install <- biocmanager.packages[! biocmanager.packages %in% packages]

if(length(biocmanager.packages2install) > 0){
  for (x in biocmanager.packages2install){BiocManager::install(x)} 
}

if(length(cran.packages2install) > 0){
  for (x in cran.packages2install){install.packages(x)} 
}

library(tximport)
library(readr)
library(DESeq2)
library(GenomicFeatures)
library(ggplot2)
library(dplyr)
library(biomaRt)
library(pheatmap)
library(openxlsx)
library(RColorBrewer)

#get parsed arguments
args <- commandArgs(trailingOnly=TRUE)
work.dir <- args[1]
gtf <- args[2]
script.dir <- args[3]
species <- args[4]
pvalue <- as.numeric(args[5])
genome <- args[6]


if (dir.exists(file.path(work.dir,"salmon"))){
  #Create sample files for all comparisons
  samples.master <- read.csv(file.path(work.dir,"samples.csv"), header=TRUE)
  number.exp <- ncol(samples.master)-2
  
  if (number.exp == 1){
    exp.names <- colnames(samples.master)[3]
  } else if(number.exp == 0) {
    print("ERROR: No comparisons defined for DEG analysis in samples.csv")
    exit()
  } else {
    exp.names <- colnames(samples.master[,3:ncol(samples.master)])
  }
  
  
  df.list <- list()#list for storing comparison dfs
  for (i in 1:length(exp.names)) { #create data frames for each experiment from master sample sheet
    exp <- exp.names[i]
    df.temp <- samples.master
    temp.list <- exp.names
    
    #select which columns to remove
    to.remove <- exp != temp.list
    drop.columns <- temp.list[to.remove]
    
    #create data frame with only one experiment (exp)
    df.temp <- df.temp[, !(names(df.temp) %in% drop.columns)]
    
    #remove NAs (keeps relevant experimental settings)
    df.temp <- df.temp[complete.cases(df.temp), ]
    
    #add df to df.list
    df.list[[i]] <- df.temp
  }
  
  names(df.list) <- exp.names #name dfs in list
  
  #Generate DESeq2 required variables
  #Create txdb from GTF file if it does not exist
  txdb.filename <- paste0(gsub(".gtf","",gtf, fixed=TRUE),".txdb")
  
  if(file.exists(txdb.filename) == FALSE){
    txdb <- makeTxDbFromGFF(gtf)
    saveDb(txdb, txdb.filename)
    txdb <- loadDb(txdb.filename)
  } else {txdb <- loadDb(txdb.filename)}
  
  #Create transcript to gene file
  k <- keys(txdb,keytype="TXNAME")
  tx2gene <- select(txdb,k,"GENEID","TXNAME")
  
  #create Salmon quant file list (for checking samples.csv)
  sf.list <- Sys.glob(file.path(work.dir, "salmon","*","quant.sf"))
  
  #Run DESeq2 for each sample df in df.list
  for (i in 1:length(df.list)){
    #create file list for DESeq2 input
    samples <- df.list[[i]]
    files <- file.path(work.dir,"salmon",paste0(samples$sample,"-quant"), "quant.sf")
    names(files) <- samples$sample
    
    #check if files in samples.csv match with files present:
    #files.found <- length(list.files(path=work.dir,pattern="quant.sf", recursive = TRUE))
    #if (length(files) != files.found){
    #  print("ERROR (DESeq2): number of samples in samples.csv does not match with actual samples in salmon directory")
    #  quit(save = "no")
    #}
    
    txi <- tximport(files,
                    type="salmon",
                    tx2gene=tx2gene,
                    ignoreTxVersion=TRUE)
    
    #Construct DESeq2 data set
    dds <- DESeqDataSetFromTximport(txi,
                                    colData=samples,
                                    design= ~ genotype)
    
    #Pre-filtering: remove rows with very few reads
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep,]
    
    #Set reference level: condition that is marked control in exp column
    ref <- unique(samples[samples$exp1 == "control", ]$genotype)
    dds$genotype <- relevel(dds$genotype, ref=ref)
    
    #Differential expression analysis
    dds <- DESeq(dds)
    res <- results(dds, alpha=0.01) #adjusted p-value < 0.01
    res <- res[order(res$padj),]
    
    #create directory for output
    dir.exp <- gsub("condition_","",resultsNames(dds)[2])
    dir.out <- paste0(work.dir,"/DESeq2/",dir.exp)
    if(dir.exists(dir.out) == FALSE){
      dir.create(dir.out)
    }
    
    #Generate MA plot
    resLFC <- lfcShrink(dds, 
                        coef=resultsNames(dds)[2], 
                        type="apeglm") #results for plotting (shrinkage of size effect)
    
    
    df_ma <-  plotMA(resLFC, 
                     alpha=pvalue,
                     colSig="red",
                     colLine="black",
                     ylim=c(-5,5),
                     returnData=TRUE)
    df_ma$above.limit <- abs(df_ma$lfc) > 5
    
    df_ma$lfc <- ifelse(df_ma$lfc > 5,5,df_ma$lfc)
    df_ma$lfc <- ifelse(df_ma$lfc < -5,-5,df_ma$lfc)
    
    p <- ggplot(df_ma, aes(x=`mean`,
                           y=`lfc`, 
                           color=`isDE`,
                           shape=`above.limit`)) +
      geom_point() +
      scale_x_continuous(trans='log10') +
      ylim(-5,5) + 
      scale_color_manual(values=c("black", "blue")) + 
      theme_bw(base_size = 20) + 
      guides(color="none",
             shape="none",
             size="none") +
      geom_hline(yintercept=0,
                 linetype="dashed") +
      ylab("log2(fold change)") +
      xlab("Mean of normalized counts") 
    
    ggsave(filename=file.path(dir.out,"MA-plot.pdf"),p,
           width=8,
           height=5)
    
    
    #Generate gene dispersion plot 
    pdf(file=file.path(dir.out,"dispersion-plot.pdf"))
    plotDispEsts(dds)
    dev.off()
    
    #Generate PCA plot
    vsd <- vst(dds, blind=FALSE)
    pcaData <- plotPCA(vsd, intgroup=c("genotype", "sample"), returnData=TRUE)
    percentVar <- round(100* attr(pcaData, "percentVar"))
    p <- ggplot(pcaData, aes(PC1, PC2, color=genotype,shape=sample)) +
      theme_bw(base_size = 20) +
      theme(legend.position = "right") +
      geom_point(size=5) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      coord_fixed(ratio=2.5)
    
    ggsave(paste0(dir.out,"/PCA-plot.pdf"),p)
    
    #Include normalised read counts in output table
    df <- as.data.frame(res)
    df$GENEID <- rownames(df)
    
    df.reads <- as.data.frame(counts(dds,normalized=TRUE)) #gets normalised read counts from dds
    df.reads$GENEID <- row.names(df.reads)
    df <- merge(x=df,
                y=df.reads,
                by="GENEID")
    
    #Convert Ensembl gene IDs to gene symbols
    #select ensembl data base for ensemble ID to gene symbol conversion
    
    genes <- res@rownames
    if(species == "mouse"){
      suppressMessages(library(EnsDb.Mmusculus.v79))
      gene.symbols <- ensembldb::select(EnsDb.Mmusculus.v79, 
                                        keys= genes, 
                                        keytype = "GENEID", 
                                        columns = c("SYMBOL","GENEID"))
    } else if(species == "human"){
      suppressMessages(library(EnsDb.Hsapiens.v79))
      gene.symbols <- ensembldb::select(EnsDb.Hsapiens.v79, 
                                        keys= genes, 
                                        keytype = "GENEID", 
                                        columns = c("SYMBOL","GENEID"))
    }
    
    suppressMessages(library(dplyr))
    df <- merge(x=df,
                y=gene.symbols,
                by="GENEID")
    
    #order data for padj
    df <- df[order(df$padj), ]
    
    #reorder columns
    df <- df %>%
      relocate(SYMBOL, .after=GENEID) 
    
    #write df to csv file
    write.csv(df, 
              file=paste0(dir.out,"/DESeq2-output.csv"),
              row.names=FALSE)
    
  }
} else if (dir.exists(file.path(work.dir,"bam"))){ #STAR aligned data
  #import sample table
  sampleTable <- read.csv(file.path(work.dir,"samples.csv"))
  
  #set conditions for DESeq2
  conditions <- length(unique(sampleTable$condition))
  if (conditions == 1){
    sampleTable$geno_cond <- sampleTable$genotype
  } ###to do multiple genotypes and multiple conditions at once
    
  
  #extract reference conditions
  sample_references <- sampleTable[sampleTable$ref %in% "ref",]
  sample_references <- unique(sample_references$geno_cond)
  sampleTable <- subset(sampleTable, select = c("sample","geno_cond"))
  
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
  #countMatrix <- countMatrix %>% dplyr::select(samples)
  
  #load data for gene annotation
  mart <- useMart("ensembl")
  mart <- useDataset("hsapiens_gene_ensembl", mart = mart)
  
  #check whether the order of the sample names between countMatrix and sampleTable is the same
  if(all(rownames(sampleTable) == colnames(countMatrix))){
    #output file
    excel <- file.path(work.dir,"DESeq2","DESeq2.xlsx")
    
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
      
      #apply size factors from yeast spike-in (if applicable)
      scale.factors <- file.path(work.dir, "scaleFactors.csv")
      if (file.exists(scale.factors)){
        print("Applying scale factors found in scaleFactors.csv")
        scale.factors <- read.csv(scale.factors)
        size.factors <- scale.factors
        size.factors$sizeFactors <- 1 / size.factors$scaleFactors  
        size.factors$scaleFactors <- NULL
        sizeFactors(dds) <- size.factors$sizeFactors
        excel <- file.path(work.dir,"DESeq2","DESeq2_scaled.xlsx")
      }
      
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
    
    #create output directory
    dir.create(file.path(work.dir,"DESeq2"), showWarnings = FALSE)
    
    #save data frames to separate sheets of an excel file
    #dir.create(file.path(work.dir,"DESeq2"), showWarnings = FALSE)
    write.xlsx(df.master, file = excel, asTable = TRUE)
    
    #also save data frames as separate files
    for (i in names(df.master)){
      write.csv(df.master[i],file.path(work.dir,"DESeq2",paste0("DESeq2_",i,".csv")),row.names = FALSE)
    }
    
    #generate heatmap of sample distances
    vsd <- vst(dds, blind=FALSE) #count data transformation for visualisation
    sampleDists <- dist(t(assay(vsd)))
    sampleDistMatrix <- as.matrix(sampleDists)
    replicate.number <- nrow(sampleTable) / length(levels(sampleTable$geno_cond)) 
    replicate.names <- paste0("_",1:replicate.number)
    rownames(sampleDistMatrix) <- paste0(dds$geno_cond,rep(replicate.names,length(levels(dds$geno_cond))))
    
    
    colnames(sampleDistMatrix) <- NULL
    colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
    pheatmap(sampleDistMatrix,
             clustering_distance_rows=sampleDists,
             clustering_distance_cols=sampleDists,
             col=colors,
             filename = file.path(work.dir,"DESeq2","sample_distance_heatmap.pdf"))
    
    #plot PCA
    pdf(file=file.path(work.dir,"DESeq2","PCA_plot.pdf"))
    if (conditions == 1){
      plotPCA(vsd, intgroup="geno_cond")
    } else {
      plotPCA(vsd, intgroup=c("genotype","condition"))
    }
    dev.off()
  }
  
}


