#!/usr/bin/env Rscript

#check packages
cran.packages <- c("ggplot2","ggrepel","viridis","dplyr")
installed.packages <- installed.packages()[,1]

for (i in cran.packages){
  if(!i %in% installed.packages){
    cat("R: ",i,"package missing, installing now\n")
    install.packages(i,repos='http://cran.us.r-project.org')
  }
}

library(ggplot2)
library(ggrepel)
library(viridis)
suppressWarnings(suppressMessages(library("dplyr")))

plot.bagel2 <- function(work.dir,df.file,save.path,title){
  
  df <- read.csv(file=df.file, sep="\t")
  df.top.ten <- df[1:10, ]
  df <- arrange(df, df$Gene)
  
  #plot hits
  p <- ggplot(df, aes(x=`Gene`,y=`BF`)) +
    theme_bw() +
    theme(axis.text = element_text(size=16),
          axis.title = element_text(size=16),
          plot.title = element_text(hjust = 0.5),
          axis.line = element_line(colour = "black"),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_blank(),
          panel.grid.major = element_blank(),
          axis.ticks.x = element_blank()) +
    xlab("Genes") +
    ylab("log2 Bayes Factor") +
    geom_point(df, mapping=aes(size=`Precision`, fill=`Recall`)
               ,shape=21) +
    ggtitle(title) +
    scale_fill_viridis(discrete=FALSE,
                       guide = guide_colorbar(frame.colour = "black", 
                                              ticks.colour = "black"), 
                       limits=c(0,1))+
    geom_hline(yintercept= 0, 
               linetype="dashed", 
               color = "red") +
    geom_label_repel(data = df.top.ten, 
                     aes(x = `Gene`, y = `BF`, label = `Gene`))
  
  file.plot <- paste0("bagel2-hits-",title, ".pdf")
  
  ggsave(plot=p,
         file=file.path(save.path,file.plot),
         width = 7,
         height = 5,
         useDingbats = FALSE)
  
  df <- arrange(df, df$Recall)
  
  #plot Precision-Recall plot
  pp <- ggplot(df, aes(x=`Recall`,y=`Precision`)) +
    theme_bw() +
    theme(axis.text = element_text(size=16),
          axis.title = element_text(size=16),
          plot.title = element_text(hjust = 0.5),
          axis.line = element_line(colour = "black"),
          panel.border = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_blank()) +
    xlab("Recall") +
    ylab("Precision (1-FDR)") +
    geom_line() +
    ggtitle(paste0("Precision-Recall plot ",title))
  
  plot.file <- paste0("PR-",title,".pdf")
  file.plot <- file.path(save.path,plot.file)
  
  ggsave(plot=pp,
         file=file.plot,
         width = 7,
         height = 5,
         useDingbats = FALSE)
}


plot.mageck <- function(work.dir,df.file,save.path,fdr){
  
  df <- read.csv(file=df.file, sep="\t")
  df$fdr.cutoff <- fdr
  
  #function to plot top 10 hits
  plot.hits <- function(input){
    
    if(input == "neg"){
      x <- "neg.rank"
      y <- "neg.score"
      z <- "neg.fdr"
      title <- "Drop out gene ranking"
      file <- "dropout.pdf"
    } else if(input == "pos") {
      x <- "pos.rank"
      y <- "pos.score"
      z <- "pos.fdr"
      title <- "Enrichment gene ranking"
      file <- "enrichment.pdf"
    }
    
    df$log.score <- -log10(df[[y]]) #for plotting
    
    #calculates score value for fdr cut off line in plot
    df$fdr.diff <- as.numeric(df$fdr.cutoff) - as.numeric(df[[z]])
    df$fdr.diff.abs <- abs(df$fdr.diff)
    fdr.min <- min(df$fdr.diff.abs)
    df.temp <- df[df$fdr.diff.abs == fdr.min,]
    fdr.cut.off.df <- df.temp[df.temp$log.score == max(df.temp$log.score), ]
    fdr.cut.off <- fdr.cut.off.df$log.score
    
    #determines top 10 hits for ggrepel labels
    df <- arrange(df, x)
    df.label <- df[df[[x]] %in% 1:10, ]
    
    options(ggrepel.max.overlaps = Inf)
    
    #generates plot
    p <- ggplot(df, aes_string(x = x, y = df$log.score)) +
      theme_bw() +
      theme(axis.text = element_text(size=16),
            axis.title = element_text(size=16),
            plot.title = element_text(hjust = 0.5),
            axis.line = element_line(colour = "black"),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.text.x = element_blank(),
            panel.grid.major = element_blank(),
            axis.ticks.x = element_blank()) +
      xlab("Genes") +
      ylab("-log(MAGeCK score)") +
      guides(color = "none",
             shape = "none") +
      ggtitle(title) +
      geom_point(df, 
                 alpha = 0.6,
                 shape = 1,
                 size = 5,
                 colour = "black",
                 mapping = aes_string(x = x,
                                      y = df$log.score)) + 
      geom_hline(yintercept= fdr.cut.off, 
                 linetype="dashed", 
                 color = "red") +
      annotate("text", 
               x = nrow(df)*0.95, 
               y = fdr.cut.off*1.05, 
               label = paste0("FDR < ",fdr),
               size = 5,
               colour="red") 
    
    if(input == "neg"){
      p <- p + geom_label_repel(data = df.label, 
                                aes(x = `neg.rank`, y = `log.score`, label = `id`))
    } else if(input == "pos") {
      p <- p + geom_label_repel(data = df.label, 
                                aes(x = `pos.rank`, y = `log.score`, label = `id`))
    }
    
    ggsave(plot=p,
           file=file.path(save.path,file),
           width = 7,
           height = 5,
           useDingbats = FALSE)
    
    #logFC plot
    
    #order x axis from low to high logFC
    df <- arrange(df, neg.lfc)
    level.order <- df$id
    level.order <- unique(level.order)
    
    df.top.neg <- df[1:5,]
    df.nrow <- nrow(df)
    df.nrow.minus <- nrow(df) - 4
    df.top.pos <- df[df.nrow.minus:df.nrow,]
    
    pp <- ggplot(df, aes(x = factor(id, level=level.order) , y = `neg.lfc`)) +
      theme_bw() +
      geom_point(size=5,
                 shape=21,
                 fill="grey",
                 alpha=0.5) +
      ylab("log2FC") +
      xlab("Gene") +
      theme(axis.text=element_text(size=16),
            axis.title=element_text(size=16),
            plot.title = element_text(hjust = 0.5),
            axis.line = element_line(colour = "black"),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.text.x = element_blank(),
            panel.grid.major = element_blank(),
            axis.ticks.x = element_blank())+ 
      geom_hline(yintercept= 0, 
                 color = "black") +
      geom_label_repel(size=4,
                       aes(x = `id`,
                           y = `neg.lfc`,
                           label = id), 
                       data = df.top.neg,
                       colour = "black",
                       nudge_x = 2,
                       nudge_y = 0.075)+
      geom_label_repel(size=4,
                       aes(x = `id`,
                           y = `neg.lfc`,
                           label = id), 
                       data = df.top.pos,
                       colour = "black",
                       nudge_x = 2,
                       nudge_y = 0.075)
    
    save.file <- file.path(save.path,"logFC-plot.pdf")
    ggsave(plot=pp,
           file=save.file,
           width = 7,
           height = 5,
           useDingbats = FALSE)
  }
  
  #Plots top 10 hits for dropout and enrichment
  input <- c("neg","pos")
  lapply(input, plot.hits)
}

args = commandArgs(trailingOnly=TRUE)

work.dir <- args[1]
df.file <- args[2]
test <- args[3]
save.path <- args[4]
title <- args[5]
script.dir <- args[6]
fdr <- args[7]

if(test == "mageck"){plot.mageck(work.dir,df.file,save.path,fdr)
} else if(test == "bagel2"){plot.bagel2(work.dir,df.file,save.path,title)
}