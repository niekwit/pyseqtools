#check packages
cran.packages <- c("ggplot2","DescTools", "stringr", "data.table")
installed.packages <- installed.packages()[,1]

for (i in cran.packages){
  if(!i %in% installed.packages){
    cat("R: ",i,"package missing, installing now\n")
    install.packages(i, repos = 'http://cran.us.r-project.org')
  }
}


library("DescTools")
library("ggplot2")
library("stringr")
library("data.table")

#get parsed arguments
args <- commandArgs(trailingOnly = TRUE)
work.dir <- args[1]
fasta <- args [2]

#read count file
df <- read.csv(file.path(work.dir, "count","counts-aggregated.tsv"), sep = "\t")

###LORENZ CURVE###
#compute Lorenz curves and GINI indeces
lorenz.pre <- Lc(df$pre)
lorenz.post <- Lc(df$post)

#create df for plotting
df.plot <- data.frame(matrix(NA, ncol = 3, nrow = length(lorenz.pre$p)))
column.names <- c("pre-amplification", "post-amplification", "p")
names(df.plot) <- column.names

df.plot$pre.L <- lorenz.pre$L #y value
df.plot$p <- lorenz.pre$p #x value

df.plot$post.L <- lorenz.post$L #y value

#create plot
p <- ggplot() +
  theme_bw(base_size = 14) +
  geom_line(data = df.plot, aes(x = `p`, y = `pre.L`, color = "pre-amplification")) +
  geom_line(data = df.plot, aes(x = `p`, y = `post.L`, color = "post-amplification")) +
  geom_abline(aes(color = "Ideal library"), intercept = 0, slope = 1 ) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_color_manual(values = c("pre-amplification" = "forestgreen", "post-amplification" = "red", "Ideal library" = "black")) +
  labs(color = "",
       x = "sgRNAs ranked by abundance",
       y = "Cumulative fraction of reads represented") +
  theme(legend.background=element_blank(),
        legend.key=element_blank(),
        legend.position = c(0.175, 0.9)) +
  annotate("text", x = 0.045, y = 0.7, hjust = 0, label = "Gini index:") + 
  annotate("text", x = 0.045, y = 0.65, hjust = 0, label = paste0("pre-amplification: ", sprintf(lorenz.pre$Gini, fmt = '%#.3f'))) + 
  annotate("text", x = 0.045, y = 0.6, hjust = 0, label = paste0("post-amplification: ", sprintf(lorenz.post$Gini, fmt = '%#.3f')))
  
#save plot
output.dir <- file.path(work.dir, "library-analysis")
if (!dir.exists(output.dir)){
  dir.create(output.dir)
}

ggsave(file.path(output.dir, "lorenz-curve.pdf"), p)


###NORMALISED GUIDE FREQUENCY###

#calculate total read count
pre.lib.sum <- sum(df$pre)
post.lib.sum <- sum(df$post)

#calculate normalised guide frequencies
df$pre.norm <- df$pre / pre.lib.sum
df$post.norm <- df$post / post.lib.sum

#calculate ratio
df$ratio <- df$pre.norm / df$post.norm
df <- df[order(df$ratio), ]
df$x <- 1:nrow(df) #x value for plotting

#remove all NaN and inf
df <- df[!is.infinite(df$ratio) & !is.na(df$ratio) , ]

#create plot
p <- ggplot(data = df, aes(x = `x`, y = `ratio`)) +
  theme_bw(base_size = 14) +
  geom_line() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_log10(expand = c(0,0)) +
  labs(color = "",
       x = "sgRNA",
       y = "Ratio normalised sgRNA frequency (pre/post)") 

ggsave(file.path(output.dir, "norm-guide-freq.pdf"), p)

###GC BIAS###

#create df with sgRNA names and sequences
df.temp <- read.csv(fasta, header = FALSE)
df.fasta <- data.frame(matrix(nrow = nrow(df.temp)/2))
df.fasta$sgRNA <- df.temp[grep(">", df.temp$V1), ]
df.fasta$sgRNA <- gsub(">", "",df.fasta$sgRNA)
df.fasta$sequence <- df.temp[grep(">", df.temp$V1, invert = TRUE), ]
df.fasta <- df.fasta[ , -1]
df.fasta$sgRNA <- sapply(str_split(df.fasta$sgRNA,"_sg", 2), `[`, 2)
df.fasta$sgRNA <- paste0("sg", df.fasta$sgRNA)

df.bias <- df[c("sgRNA", "pre", "post")]
df.bias <- merge(df.bias, df.fasta, by = "sgRNA")

#calculate GC content for each sgRNA
df.bias$`G-count` <- str_count(df.bias$sequence, pattern = "G")
df.bias$`C-count` <- str_count(df.bias$sequence, pattern = "C")
df.bias$`length` <- nchar(df.bias$sequence)

df.bias$`%GC` <- (df.bias$`G-count` + df.bias$`C-count`) / df.bias$`length` * 100


#df.melt <- melt(df.bias, 
#                measure.vars = c("pre", "post"),
#                variable.name = "sample",
#                value.name = "Count")

gc.content <- unique(df.bias$`%GC`)
df.plot <- data.frame(matrix(NA, ncol = 3, nrow = length(gc.content)))
names(df.plot) <- c("pre", "post", "%GC")
df.plot$`%GC` <- gc.content



total.reads.pre <- sum(df.plot$pre)
total.reads.post <- sum(df.plot$post)

df.plot$pre <- df.plot$pre / total.reads.pre
df.plot$post <- df.plot$post / total.reads.post

df.melt <- melt(df.plot,
                measure.vars = c("pre", "post"),
                variable.name = "sample",
                value.name = "Count")

df.melt <- df.melt[order(df.melt$`%GC`), ]

ggplot(df.melt, aes(x = as.character(`%GC`), y = `Count`, fill = `sample`)) +
  theme_bw(base_size = 14) +
  geom_col(position = "dodge",
           size = 0.5,
           color = "black",
           width = 0.80) +
  labs(y = "Normalised read count",
       x = "%GC") +
  scale_fill_manual(values = c("pre" = "royalblue", "post" = "orange")) +
  theme(legend.position = c(0.9, 0.9),
        legend.background=element_blank(),
        legend.key=element_blank())


