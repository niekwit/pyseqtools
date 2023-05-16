library(DiffBind)
library(rtracklayer)
library(ChIPseeker)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggupset)
library(ReactomePA)
library(stringr)
library(dplyr)
library(glue)

library(BiocParallel)
register(SerialParam())

#get parsed arguments
args <- commandArgs(trailingOnly = TRUE)
work.dir <- args[1]
genome <- args[2]

#create output dir
out.dir <- file.path(work.dir,"diffbind")
dir.create(out.dir,showWarnings = FALSE)

#load appropriate genome db
if (genome == "hg38"){
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  } else if (genome == "hg19"){
      library(TxDb.Hsapiens.UCSC.hg19.knownGene)
      txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    }

####prepare sample sheet####

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

#use sample_group for Condition (simplifies Diffbind experiment design)
conditions <- unique(sampleInfo$sample_group)
for (i in conditions){
  geno <- str_split(i,"_")[[1]][1]
  treatment <- str_split(i,"_")[[1]][2]
  
  for (j in 1: nrow(diffbind.sheet)){
    if (str_detect(diffbind.sheet[j,"SampleID"], geno)){
      if (str_detect(diffbind.sheet[j,"Treatment"], treatment)){
        diffbind.sheet[j,"Condition"] <- i
      }
    }
  }
}

#add replicate numbers to diffbind.sheet
diffbind.sheet$temp <- paste0(diffbind.sheet$Condition,"_",diffbind.sheet$Treatment)
diffbind.sheet <- diffbind.sheet %>% 
  group_by(temp) %>% 
  mutate(Replicate = row_number())
diffbind.sheet$temp <- NULL

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
  peak <- file.path(work.dir,"peaks",genome,i,paste0(i,"_peaks.xls"))
  diffbind.sheet[diffbind.sheet$SampleID == i,10] <- peak
}

#write diffbind df to file
sample.sheet <- file.path(out.dir,"diffbind_samples.csv")
write.csv(diffbind.sheet,
          file=sample.sheet,
          row.names = FALSE)


####differential binding analysis####

diffbind <- dba.analyze(diffbind.sheet,
                        bGreylist=FALSE,
                        design="~Condition")

#diffbind <-dba(sampleSheet=diffbind.sheet) %>%
#            dba.blacklist() %>%
#            dba.count(bParallel=FALSE) %>%
#            dba.normalize() %>%
#            dba.contrast(minMembers = 2) %>%
#            dba.analyze(design="~Condition",
#                        bGreylist=FALSE)

#save diffbind object to file (e.g. for tailored analyses)
save(diffbind, file = file.path(work.dir,"diffbind",glue("diffbind.Rdata")))


#if (control.condition == "None"){
#  diffbind <- dba.analyze(diffbind.sheet,
#                        bGreylist=FALSE)
#} else {
#    if (control.condition %in% sampleInfo$sample_group ==TRUE){
#      diffbind <- dba.analyze(diffbind.sheet,
#                        bGreylist=FALSE,
#                        design="~Condition")
#
#    #set contrasts
#    diffbind <- dba.contrast(diffbind, reorderMeta=list(Condition=control.condition))
#    } else {
#      stop("ERROR: parsed control condition not found in samples.csv")
#    }
#  }

#plot correlation heatmap
pdf(file.path(out.dir,"sample_correlation_heatmap.pdf"))
plot(diffbind)
dev.off()

#plot PCA
pdf(file.path(out.dir,"PCA.pdf"))
dba.plotPCA(diffbind,DBA_CONDITION,label=DBA_ID)
dev.off()

#generate binding affinity heatmaps
hmap <- colorRampPalette(c("red", "black", "forestgreen"))(n = 13)
pdf(file.path(out.dir,paste0("affinity_heatmap.pdf")))
dba.plotHeatmap(diffbind, correlations=FALSE,scale="row", colScheme = hmap)
dev.off()

#Plot profiles
tryCatch({
  pdf(file.path(out.dir,"plotProfile.pdf"))
  dba.plotProfile(diffbind, samples=diffbind$masks$All, doPlot=TRUE)
  dev.off()
}, error = function(e){
  warning("WARNING: failed to profile plot")
})

#for each condition perform differential binding analysis as control
for (condition in conditions){
  #set control condition
  diffbind <- dba.contrast(diffbind, reorderMeta=list(Condition=condition))
  
  #get number of contrasts
  contrasts <- length(diffbind$contrasts)

  #iterate over contrasts and generate plots when DBs > 1
  for(i in 1:contrasts){
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
      file <- file.path(out.dir,condition,paste0("MA-plot-",condition,".pdf"))
      pdf(file)
      dba.plotMA(diffbind, contrast = i)
      dev.off()
      
      #Generate volcano plots
      print("Volcano plots")
      file <- file.path(out.dir,condition,paste0("volcano-",condition,".pdf"))
      pdf(file)
      dba.plotVolcano(diffbind, contrast = i)
      dev.off()
      
      #Plot venn diagrams
      tryCatch({
        print("Creating Venn diagram")
        file <- file.path(out.dir,condition,paste0("venn-",condition,".pdf"))
        pdf(file)
        dba.plotVenn(diffbind, contrast=i, bDB=TRUE,bGain=TRUE,bLoss=TRUE,bAll=FALSE)
        dev.off()
      }, error = function(e){
        print("ERROR: failed to create Venn diagram")
      }, warning = function(w){
        print("WARNING: failed to create Venn diagram")
      })

      #export DBs to bed file
      print(paste0("Exporting differential binding sites to a bed file for ", contrast[2]," vs ",contrast[3]))
      bed.file <- file.path(out.dir,condition,paste0("DB-",condition,".bed"))
      export.bed(dbs,con=bed.file)
      
      #prepend "chr" to chromosome names in bed file
      system(glue("sed -i 's/^/chr/' {bed.file}"))

      #coverage plot
      peak <- readPeakFile(bed.file)
      file <- file.path(out.dir,condition,paste0("coverage-",condition,".pdf"))
      pdf(file)
      covplot(peak, weightCol="V4")
      dev.off()
      
            peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
      out <- as.data.frame(peakAnno)
      write.csv(out,
              file = file.path(out.dir,condition,glue("DB-{condition}_annotated.csv")),
              col.names = FALSE)

      file <- file.path(out.dir,condition,paste0("piechart-",condition,".pdf"))
      pdf(file)
      plotAnnoPie(peakAnno)
      dev.off()
      
      file <- file.path(out.dir,condition,paste0("upsetplot-",condition,".pdf"))
      pdf(file)
      upsetplot(peakAnno, vennpie=TRUE)
      dev.off()
      
      ###Functional enrichment analysis with ReactomePA
      #all peaks
      pathways_all <- enrichPathway(as.data.frame(peakAnno)$geneId)
      
      tryCatch({
        file <- file.path(out.dir,condition,paste0("dotplot_all_peaks-",condition,".pdf"))
        pdf(file)
        dotplot(pathways_all)
        dev.off()
      }, error = function(e){
        print("ERROR: failed to create dotplot_all_peaks")
      }, warning = function(w){
        print("WARNING: failed to create dotplot_all_peaks")
      })
            
      go.all <- pathways_all@result
      go.all$geneSymbol <- NA
      suppressMessages(for(j in 1:nrow(go.all)){
        entrez.genes <- unlist(strsplit(go.all[j,8], "/"))
        entrez.genes <- mapIds(org.Hs.eg.db, entrez.genes,'SYMBOL','ENTREZID')
        entrez.genes <- paste(entrez.genes, collapse = "/")
        go.all[j,10] <- entrez.genes
      })
      write.csv(go.all, file = file.path(out.dir,condition,paste0("GO_all-peaks_",condition,".csv")))
      
      #peaks around tss
      gene <- seq2gene(peak, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)
      pathways_tss <- enrichPathway(gene)
      
      tryCatch({
        file <- file.path(out.dir,condition,paste0("dotplot_tss-",condition,".pdf"))
        pdf(file)
        dotplot(pathways_tss)
        dev.off()
      }, error = function(e){
        print("ERROR: failed to create dotplot_tss")
      }, warning = function(w){
        print("WARNING: failed to create dotplot_tss")
      })


      go.peak <- pathways_all@result
      go.peak$geneSymbol <- NA
      suppressMessages(for(j in 1:nrow(go.peak)){
        entrez.genes <- unlist(strsplit(go.peak[j,8], "/"))
        entrez.genes <- mapIds(org.Hs.eg.db, entrez.genes,'SYMBOL','ENTREZID')
        entrez.genes <- paste(entrez.genes, collapse = "/")
        go.peak[j,10] <- entrez.genes
      })
      write.csv(go.all, file = file.path(out.dir,condition,paste0("GO_tss-peaks_",condition,".csv")))
    }
  }

}




