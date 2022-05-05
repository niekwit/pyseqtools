library(DiffBind)

#setwd("/home/niek/Documents/analyses/ChIP-Seq/COVID/bam/H3K4me3")
sample.sheet <-file.path(work.dir,"diffbind_samples.csv")

if(file.exists(sample.sheet)){
  samples <- read.csv(sample.sheet)
  diffbind <- dba.analyze(samples)
  
  #get number of conditions
  conditions <- length(diffbind$contrasts)
  
  #plot correlation heatmap
  pdf(file.path(work.dir,"diffbind","sample_correlation_heatmap.pdf"))
  plot(diffbind)
  dev.off()
  
  #plot PCA
  pdf(file.path(work.dir,"diffbind","PCA.pdf"))
  dba.plotPCA(diffbind,DBA_CONDITION,label=DBA_ID)
  dev.off()
  
  #Generate MA plots
  for(i in 1:conditions){
    contrast <- diffbind$contrasts[[i]][[5]]
    condition <- paste0(contrast[2],"_vs_",contrast[3])
    file <- file.path(work.dir,"diffbind",paste0("MA-plot-",condition,".pdf"))
    
    pdf(file)
    dba.plotMA(diffbind, contrast = i)
    dev.off()
  }
  
  
  
  #Generate volcano plots
  for(i in 1:conditions){
    contrast <- diffbind$contrasts[[i]][[5]]
    condition <- paste0(contrast[2],"_vs_",contrast[3])
    file <- file.path(work.dir,"diffbind",paste0("volcano-",condition,".pdf"))
    
    pdf(file)
    dba.plotVolcano(diffbind, contrast = i)
    dev.off()
  }
  
  #generate binding affinity heatmaps
  hmap <- colorRampPalette(c("red", "black", "green"))(n = 13)
  pdf(file.path(work.dir,"diffbind",paste0("affinity_heatmap.pdf")))
  dba.plotHeatmap(tamoxifen, correlations=FALSE,scale="row", colScheme = hmap)
  dev.off()
  
  #Plot profiles
  
  pdf(file.path(work.dir,"diffbind","plotProfile.pdf"))
  profiles <- dba.plotProfile(diffbind, doPlot=TRUE)
  dev.off()
  
  
  #get DBs
  diffbind.DB <- dba.report(diffbind)
  
}else{
  print("Diffbind sample sheet not found")
}








