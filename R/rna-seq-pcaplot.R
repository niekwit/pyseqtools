#check if required packages are installed, if not install them
packages <- rownames(installed.packages())
cran.packages <- c("ggfortify","ggplot2",
                   "stringr","data.table")

cran.packages2install <- cran.packages[! cran.packages %in% packages]

if(length(cran.packages2install) > 0){
  for (x in cran.packages2install){install.packages(x)} 
}

library(ggfortify)
library(ggplot2)
library(stringr)
library(data.table)

#get parsed arguments
args <- commandArgs(trailingOnly = TRUE)

#get all Salmon quant files in experimental directory
work.dir <- args[1]
file.list <- Sys.glob(file.path(work.dir,"salmon","*","quant.sf"))

#get number of transcripts
nrows <- nrow(read.csv(file.list[1]))

#create empty data frame to store TPM values of all samples
df.all <- data.frame(matrix(nrow = nrows, ncol = length(file.list)))
column.names <- basename(dirname(file.list))
column.names <- str_replace_all(column.names, "-quant","")
colnames(df.all) <- column.names

#add TPM values to empty data frame
for (i in 1:length(file.list)){
  df <- read.csv(file.list[i], sep="\t")
  df.all[i] <- df$TPM
}

#perform PCA analysis
df.all.t <- as.data.frame(t(df.all)) #transposes data
pca_res <- prcomp(df.all.t)

#create sample and condition names
df.all.t$samples <- row.names(df.all.t)
df.all.t$Condition <- gsub('.{2}$', 
                           '', 
                           x = df.all.t$samples) #removes -1/-2 from sample names (last two characters)

#check number of replicates for shapes plots
sample.number <- length(df.all.t$samples)
unique.samples <- length(unique(df.all.t$Condition))
number.replicates <- sample.number / unique.samples

shape.list <- c(15,16,17,18,4,3,8,11,10,13,21,22,23,24,25)

shape.value.vector <- shape.list[1:number.replicates]

#plot PCA
p <- autoplot(pca_res, 
         data = df.all.t, 
         colour = "Condition",
         shape = "samples",
         size = 6,
         ) +
  theme_set(theme_bw()) +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size=12)
        ) +
  scale_shape_manual(name = "Samples",
                     values = c(rep(shape.value.vector, 
                                  length(unique(df.all.t$Condition))))) +
  coord_fixed(ratio = 0.75)

ggsave(file.path(file = work.dir,"salmon","PCAplot.pdf"), 
       plot = p)
