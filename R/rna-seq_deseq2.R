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

suppressMessages(library(tximport))
suppressMessages(library(readr))
suppressMessages(library(DESeq2))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(ggplot2))

#get parsed arguments
args <- commandArgs(trailingOnly=TRUE)
work.dir <- args[1]
gtf <- args[2]
script.dir <- args[3]
species <- args[4]
pvalue <- as.numeric(args[5])

#Create sample files for all comparisons
samples.master <- read.csv(file.path(work.dir,"samples.csv"), header=TRUE)
number.exp <- ncol(samples.master)-2
exp.names <- colnames(samples.master[,3:ncol(samples.master)])

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

#Run DESeq2 for each sample df in df.list
for (i in 1:length(df.list)){
  #create file list for DESeq2 input
  samples <- df.list[[i]]
  files <- file.path(work.dir,"salmon",paste0(samples$sample,"-quant"), "quant.sf")
  names(files) <- samples$sample
  
  txi <- tximport(files,
                  type="salmon",
                  tx2gene=tx2gene,
                  ignoreTxVersion=TRUE)
  
  #Construct DESeq2 data set
  dds <- DESeqDataSetFromTximport(txi,
                                  colData=samples,
                                  design= ~ condition)
  
  #Pre-filtering: remove rows with very few reads
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  
  #Set reference level: condition that is marked control in exp column
  ref <- which(samples == "control", arr.ind=TRUE) #select rows that contain control
  ref <- ref[,"row"]#get indeces of rows that contain control
  ref <- samples[ref[[1]],] #subset samples for row that contains just control condition
  ref <- ref[["condition"]] #extract reference condition
  dds$condition <- relevel(dds$condition, ref=ref)
  
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
  pcaData <- plotPCA(vsd, intgroup=c("condition", "sample"), returnData=TRUE)
  percentVar <- round(100* attr(pcaData, "percentVar"))
  p <- ggplot(pcaData, aes(PC1, PC2, color=condition,shape=sample)) +
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
            file=paste0(dir.out,"/DESeq-output.csv"),
            row.names=FALSE)
  
}
