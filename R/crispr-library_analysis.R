library("DescTools")
library("ggplot2")

#get parsed arguments
args <- commandArgs(trailingOnly = TRUE)
work.dir <- args[1]

#read count file
df <- read.csv(file.path(work.dir, "count","counts-aggregated.tsv"), sep = "\t")

#compute Lorenz curves and GINI indeces
lorenz.pre <- Lc(df$pre)
lorenz.post <- Lc(df$post)

#create df for plotting
df.plot <- data.frame(matrix(NA, ncol = 3, nrow = length(lorenz.pre$p)))
column.names <- c("pre.L", "post.L", "p")
names(df.plot) <- column.names

df.plot$pre.L <- lorentz.pre$L #y value
df.plot$p <- lorentz.pre$p #x value

df.plot$post.L <- lorentz.post$L #y value

#create plot
p <- ggplot() +
  theme_bw(base_size = 14) +
  geom_line(data = df.plot, aes(x = `p`, y = `pre.L`, color = "pre.L")) +
  geom_line(data = df.plot, aes(x = `p`, y = `post.L`, color = "post.L")) +
  geom_abline(aes(color = "Ideal library"), intercept = 0, slope = 1 ) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_color_manual(values = c("pre.L" = "forestgreen", "post.L" = "red", "Ideal library" = "black")) +
  labs(color = "",
       x = "sgRNAs ranked by abundance",
       y = "Cumulative fraction of reads represented") +
  theme(legend.background=element_blank(),
        legend.key=element_blank(),
        legend.position = c(0.1, 0.9)) +
  annotate("text", x = 0.080, y = 0.7, label = "Gini index:") + 
  annotate("text", x = 0.138, y = 0.65, label = paste0("pre-amplification: ", sprintf(lorentz.pre$Gini, fmt = '%#.3f'))) + 
  annotate("text", x = 0.141, y = 0.6, label = paste0("post-amplification: ", sprintf(lorentz.post$Gini, fmt = '%#.3f')))
  
#save plot
output.dir <- file.path(work.dir, "library-analysis")
if (!dir.exists(output.dir)){
  dir.create(output.dir)
}


ggsave(file.path(output.dir, "lorenz-curve.pdf"), p)


