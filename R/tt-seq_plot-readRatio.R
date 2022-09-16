library(ggplot2)
library(dplyr)
library(reshape2)
library(biomaRt)

work.dir <- getwd()

#get parsed arguments
args <- commandArgs(trailingOnly=TRUE)
genome <- args[1]

ensembl_gene_id <- scan(file="/home/niek/Documents/references/bed/hg38/Homo_sapiens.GRCh38.95_gene-only_TSS_1kb_TES_1kb_bedtools-sorted_ENSEMBL_gene-names_only.txt", what="character")

df <- as.data.frame(ensembl_gene_id)
df.ratio <- as.data.frame(ensembl_gene_id)

count.list <- Sys.glob(file.path(work.dir,"readRatio",genome,"*.txt"))

#get unique sample names
samples <- read.csv(file.path(work.dir,"samples.csv"))
samples <- samples$sample
sample_split <- function(x) {strsplit(x, split="_[1-9]")[[1]][1]}
samples <- unlist(unique(lapply(samples, sample_split)))

for (i in samples){
  TSS.count <- Sys.glob(file.path(work.dir,"readRatio",genome, paste0(i,"*_TSS.txt")))
  TES.count <- Sys.glob(file.path(work.dir,"readRatio",genome, paste0(i,"*_TES.txt")))
  df.tss <- read.csv(TSS.count, header = FALSE, sep = " ")
  names(df.tss) <- c(paste0(i,"_TSS"),"ensembl_gene_id")
  df.tes <- read.csv(TES.count, header = FALSE, sep = " ")
  names(df.tes) <- c(paste0(i,"_TES"),"ensembl_gene_id")
  
  df.temp <- left_join(df, df.tss, by="ensembl_gene_id")
  df.temp <- left_join(df.temp, df.tes, by="ensembl_gene_id")
  
  df.temp[4] <- df.temp[3] / df.temp[2]
  names(df.temp)[4] <- i
  df.temp <- na.omit(df.temp)
  df.temp[2:3] <- NULL
  
  df.ratio <- left_join(df.ratio, df.temp, by="ensembl_gene_id")
  
}

df.plot <- melt(df.ratio)

p <- ggplot(df.plot, aes(x=variable, y=value, fill=variable)) + 
  scale_fill_manual(values=c("lightskyblue", "tomato", "lightskyblue", "tomato")) +
  geom_violin(trim=FALSE) +
  geom_boxplot(width=0.1) +
  theme_bw() +
  scale_y_continuous(trans='log10') +
  labs(x = NULL,
       y = "TES / TSS ratio") +
  labs(title = "Global TES/TES ratios for all genes") + 
  theme(text=element_text(size=20)) + 
  guides(fill="none")


#load data for gene annotation
mart <- useMart("ensembl")
mart <- useDataset("hsapiens_gene_ensembl", mart = mart)

genes.table <- getBM(filters = "ensembl_gene_id", 
                     attributes = c("ensembl_gene_id", "external_gene_name","description",
                                                                 "gene_biotype", "chromosome_name","start_position",
                                                                 "end_position", "percentage_gene_gc_content"), 
                     values = df.ratio$ensembl_gene_id, 
                     mart = mart)
genes.table$gene_length <- genes.table$end_position - genes.table$start_position

df.ratio <- left_join(df.ratio, genes.table, by="ensembl_gene_id")

write.table(df.ratio, 
            file=file.path(work.dir, "readRatio", genome ,"readRatios.txt"), 
            quote=FALSE, 
            row.names=FALSE,
            sep="\t")

if (length(count.list) == 0){
  bam.list <- Sys.glob(file.path(work.dir, "bam", genome, "*bam"))
} else{
  
}
  

