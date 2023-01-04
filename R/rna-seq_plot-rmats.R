library(tidyverse)
library(reshape2)

#get parsed arguments
args <- commandArgs(trailingOnly = TRUE)
work.dir <- args[1]

#splice events
se <- c("A3SS","A5SS","MXE","RI","SE")
base.file <- ".MATS.JCEC.txt"

#get rMATS data dirs
rmats.dirs <- Sys.glob(file.path(work.dir,"rmats","*","*"))
rmats.dirs <-rmats.dirs[dir.exists(rmats.dirs)] #subset for only dirs

#plot data from each dir
for (i in rmats.dirs){
  #get sample names
  control <- str_split_fixed(basename(i),"_vs_",2)[1]
  test <- str_split_fixed(basename(i),"_vs_",2)[2]
  
  #set colours for plotting
  colours <- unname(palette.colors(palette = "Classic Tableau",n=5))
  
  #plot bar graph with summary data
  summary <- read.csv(file.path(i,"summary.txt"), sep="\t")
  summary <- subset(summary, select=c("EventType","SigEventsJCECSample1HigherInclusion","SigEventsJCECSample2HigherInclusion"))
  colnames(summary) <- c("EventType",control,test)
  
  summary <- melt(summary)
  
  p <- ggplot(summary, aes(fill=EventType,y=value,x=variable)) +
    geom_bar(position="stack", stat="identity",colour="black") +
    theme_bw() +
    xlab(NULL) +
    ylab("Number of significant events") +
    scale_fill_manual(values = colours) +
    theme(text=element_text(size=18))
  
  ggsave(file.path(i,"summary-plot.pdf"),p,width=6,height=6)
  
  #create pie chart with differential splice events
  df.se <- data.frame(matrix(ncol = 2, nrow = 5))
  names(df.se) <- c("EventType","Number")
  df.se$EventType <- se
  
  for (j in se){
    se.file = file.path(i,paste0(j,base.file))
    df <- read.csv(se.file,sep="\t")
    df <- subset(df, FDR < 0.1 & IncLevelDifference >= 0.3 | FDR < 0.1 & IncLevelDifference <= -0.3)
    df.se[df.se$EventType == j,2] = nrow(df)
    
    write_csv(df,file=file.path(i,paste0(j,".MATS.JCEC_filtered.csv"))) #save significant diff events to csv
    }
  
  pp <- ggplot(df.se, aes(x = "", y = Number, fill = EventType)) +
    geom_col(color = "black") +
    geom_text(aes(label = Number),
              position = position_stack(vjust = 0.5),
              size = 14) +
    coord_polar(theta = "y") +
    scale_fill_manual(values = colours) +
    theme_void() +
    theme(text=element_text(size=26),
          legend.key.height=unit(1.75,"cm"),
          legend.key.width=unit(1.75,"cm"),
          plot.title = element_text(hjust = 0.5)) +
    labs(fill="Number of\ndifferential\nsplice events",
         title=basename(i)) 
  
  ggsave(file.path(i,paste0(basename(i),"_pie-chart_splice-events.pdf")))
}






