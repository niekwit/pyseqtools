#!/usr/bin/env Rscript

#insprired by https://rjbioinformatics.com/2016/10/14/converting-mouse-to-human-gene-names-with-biomart-package/

#check for biomaRt package
packages <- rownames(installed.packages())

if ("biomaRt" %in% packages == FALSE){
  BiocManager::install(biomaRt)
}



library(biomaRt)

#get parsed arguments
args <- commandArgs(trailingOnly = TRUE)

conversion <- args[1]
gene.list <- args[2]
output <- args[3]


if (conversion == "hm"){
  human.genes <- read.csv(gene.list, header = FALSE)
  human.genes <- human.genes$V1
  
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", 
                   values = human.genes , mart = human, attributesL = c("mgi_symbol"), 
                   martL = mouse, uniqueRows=T)
  
  mouse.genes <- unique(genesV2[, 2])
  
  write.table(mouse.genes,
              file = output,
              col.names = FALSE,
              row.names = FALSE,
              quote = FALSE)
} else if (conversion == "mh"){
  mouse.genes <- read.csv(gene.list, header = FALSE)
  mouse.genes <- mouse.genes$V1
  
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", 
                   values = mouse.genes , mart = mouse, attributesL = c("hgnc_symbol"), 
                   martL = human, uniqueRows=T)
  human.genes <- unique(genesV2[, 2])
  
  write.table(human.genes,
              file = output,
              col.names = FALSE,
              row.names = FALSE,
              quote = FALSE)
} else {
  print("ERROR: invalid conversion selected")
  
}





