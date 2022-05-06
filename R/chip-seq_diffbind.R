#get parsed arguments
args <- commandArgs(trailingOnly = TRUE)
work.dir <- args[1]
genome <- args[2]

#check for packages
packages <- rownames(installed.packages())
biocmanager.packages <- c("tximport","DiffBind",
                          "GenomicFeatures","rtracklayer",
                          "clusterProfiler","ggupset", 
                          "ReactomePA", "TxDb.Hsapiens.UCSC.hg38.knownGene",
                          "org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg19.knownGene",
                          "ggimage", "magick"
                          )
biocmanager.packages2install <- biocmanager.packages[! biocmanager.packages %in% packages]

if(length(biocmanager.packages2install) > 0){
  for (x in biocmanager.packages2install){BiocManager::install(x)} 
}

library(DiffBind)
library(rtracklayer)

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

library(clusterProfiler)
library(ggupset)
library(ReactomePA)




#setwd("/home/niek/Documents/analyses/ChIP-Seq/COVID/bam/H3K4me3")
sample.sheet <-file.path(work.dir,"diffbind_samples.csv")

if(file.exists(sample.sheet)){
  samples <- read.csv(sample.sheet)
  print("Finding differential peaks with DiffBind")
  
  #Find differential binding sites (DBs)
  diffbind <- dba.analyze(samples)
  
  #plot correlation heatmap
  pdf(file.path(work.dir,"diffbind","sample_correlation_heatmap.pdf"))
  plot(diffbind)
  dev.off()
  
  #plot PCA
  pdf(file.path(work.dir,"diffbind","PCA.pdf"))
  dba.plotPCA(diffbind,DBA_CONDITION,label=DBA_ID)
  dev.off()
  
  #generate binding affinity heatmaps
  hmap <- colorRampPalette(c("red", "black", "green"))(n = 13)
  pdf(file.path(work.dir,"diffbind",paste0("affinity_heatmap.pdf")))
  dba.plotHeatmap(tamoxifen, correlations=FALSE,scale="row", colScheme = hmap)
  dev.off()
  
  #Plot profiles
  pdf(file.path(work.dir,"diffbind","plotProfile.pdf"))
  profiles <- dba.plotProfile(diffbind, samples=diffbind$masks$All, doPlot=TRUE)
  dev.off()
  
  #get number of conditions
  conditions <- length(diffbind$contrasts)
  
  #iterate over conditions (contrasts) and generate plots when DBs > 1
  for(i in 1:conditions){
    dbs <- dba.report(diffbind, contrast=i)
    if(sum(dbs$Fold>0) > 1 | sum(dbs$Fold<0) > 1){
      contrast <- diffbind$contrasts[[i]][[5]]
      condition <- paste0(contrast[2],"_vs_",contrast[3])
      print(paste0("Generating plots for ", contrast[2]," vs ",contrast[3]))
      
      #create dir for output
      dir.out <- file.path(work.dir,"diffbind", condition)
      if(dir.exists(dir.out) == FALSE){
        dir.create(dir.out)
      }
      
      
      #Generate MA plots
      print("MA plots")
      file <- file.path(work.dir,"diffbind",condition,paste0("MA-plot-",condition,".pdf"))
      pdf(file)
      dba.plotMA(diffbind, contrast = i)
      dev.off()
      
      #Generate volcano plots
      print("Volcano plots")
      file <- file.path(work.dir,"diffbind",condition,paste0("volcano-",condition,".pdf"))
      pdf(file)
      dba.plotVolcano(diffbind, contrast = i)
      dev.off()
      
      #Plot venn diagrams
      print("Venn diagrams")
      file <- file.path(work.dir,"diffbind",condition,paste0("venn-",condition,".pdf"))
      pdf(file)
      dba.plotVenn(tamoxifen, contrast=i, bDB=TRUE,bGain=TRUE,bLoss=TRUE,bAll=FALSE)
      dev.off()
      
      #export DBs to bed file
      print(paste0("Exporting differential binding sites to a bed file for ", contrast[2]," vs ",contrast[3]))
      bed.file <- file.path(work.dir,"diffbind",condition,paste0("DB-",condition,".bed"))
      export.bed(dbs,con=bed.file)
    
      #coverage plot
      peak <- readPeakFile(bed.file)
      file <- file.path(work.dir,"diffbind",condition,paste0("coverage-",condition,".pdf"))
      pdf(file)
      covplot(peak, weightCol="V4")
      dev.off()
      
      #peak annotation
      txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
      peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
      
      file <- file.path(work.dir,"diffbind",condition,paste0("piechart-",condition,".pdf"))
      pdf(file)
      plotAnnoPie(peakAnno)
      dev.off()
      
      file <- file.path(work.dir,"diffbind",condition,paste0("upsetplot-",condition,".pdf"))
      pdf(file)
      upsetplot(peakAnno, vennpie=TRUE)
      dev.off()
      
      #Functional enrichment analysis with ReactomePA
      pathways_all <- enrichPathway(as.data.frame(peakAnno)$geneId)
      file <- file.path(work.dir,"diffbind",condition,paste0("dotplot_all_peaks-",condition,".pdf"))
      pdf(file)
      dotplot(pathway_all)
      dev.off()
      
      gene <- seq2gene(peak, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)
      pathways_tss <- enrichPathway(gene)
      file <- file.path(work.dir,"diffbind",condition,paste0("dotplot_tss-",condition,".pdf"))
      pdf(file)
      dotplot(pathways_tss)
      dev.off()
    
      
  }
  
  }else{
    print("WARNING: Diffbind sample sheet not found, skipping differential binding site analysis")
  }
}
  
  
  