library(DiffBind)

setwd("/home/niek/Documents/analyses/ChIP-Seq/COVID/bam/H3K4me3")
samples <- read.csv("diffbind_samples.csv")
diffbind <- dba.analyze(samples)



setwd(system.file('extra',package='DiffBind'))
tmpdir <- tempdir()
url <- 'https://content.cruk.cam.ac.uk/bioinformatics/software/DiffBind/DiffBind_vignette_data.tar.gz'
file <- basename(url)
options(timeout=600)
download.file(url, file.path(tmpdir,file))
untar(file.path(tmpdir,file), exdir = tmpdir )
setwd(file.path(tmpdir,"DiffBind_Vignette"))

tamoxifen <- dba.analyze("tamoxifen.csv")
