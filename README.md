# pyseqtools
Tool box for a variety of NGS data analyses

Available analyses:

- CRISPR screens
- RNA-Seq
- ChIP-Seq
- CUT & RUN

## Installation

> `git clone https://github.com/niekwit/pyseqtools.git`

When running an analysis, `pyseqtools.py` will first check whether the relevant dependencies are available in `$PATH`. If not, it will then look for it in any `$HOME` directories, and finally if it does not find them there it will install them.
## Usage:

```
pyseqtools.py [-h] {crispr,rna-seq,chip-seq,cutrun} ...

Analysis of a variety of NGS data

positional arguments:
  {crispr,rna-seq,chip-seq,cutrun}
                        Available analyses:
    crispr              CRISPR screens
    rna-seq             RNA-Seq
    chip-seq            ChIP-Seq
    cutrun              CUT & RUN

optional arguments:
  -h, --help            show this help message and exit

```

### CRISPR screen analysis

```
pyseqtools.py crispr [-h] -l
                            {CRISPR library}
                            [-t <int>] [-r] [-m N] [-a {mageck,bagel2}] [-f <FDR value>] [--cnv <CCLE cell line>] [--go]
                            [--gene-sets <GO gene set>] [--essential-genes <Custom essential gene list>] [--csv2fasta <CSV FILE>] [--skip-fastqc]
                            [--skip-stats]

optional arguments:
  -h, --help            show this help message and exit
  -l {CRISPR library}, --library {CRISPR library}
                        CRISPR library
  -t <int>, --threads <int>
                        Number of CPU threads to use (default is 1). Use max to apply all available CPU threads
  -r, --rename          Rename fastq files according to rename.config
  -m N, --mismatch N    Number of mismatches (0 or 1) allowed during alignment
  -a {mageck,bagel2}, --analysis {mageck,bagel2}
                        Statistical analysis with MAGeCK or BAGEL2 (default is MAGeCK)
  -f <FDR value>, --fdr <FDR value>
                        Set FDR cut off for MAGeCK hits (default is 0.25)
  --cnv <CCLE cell line>
                        Activate CNV correction for MAGeCK/BAGEL2 with given cell line
  --go                  Gene set enrichment analysis with enrichR
  --gene-sets <GO gene set>
                        Gene sets used for GO analysis (default is GO_Molecular_Function_2021, GO_Cellular_Component_2021, and
                        GO_Biological_Process_2021). Gene sets can be found on https://maayanlab.cloud/Enrichr/#stats
  --essential-genes <Custom essential gene list>
                        Essential gene list (default is Hart et al 2015 Cell)
  --csv2fasta <CSV FILE>
                        Convert CSV file with sgRNA names and sequences to fasta format. The first column should contain sgRNA names and the
                        second sgRNA sequences (headers can be anything).
  --skip-fastqc         Skip FastQC/MultiQC
  --skip-stats          Skip MAGeCK/BAGEL2

```

### RNA-Seq analysis

```
pyseqtools.py rna-seq [-h] [-t <int>] [-r {gencode-v35,hg19,hg38,mm10,mm9}] -s {mouse,human} [-a {salmon,hisat2}] [-p <P value>] [--go]
                             [--gene-sets <GO gene set>] [--skip-fastqc]

optional arguments:
  -h, --help            show this help message and exit
  -t <int>, --threads <int>
                        Number of CPU threads to use (default is 1). Use max to apply all available CPU threads. For Salmon 8-12 threads are
                        optimal
  -r {gencode-v35,hg19,hg38,mm10,mm9}, --reference {gencode-v35,hg19,hg38,mm10,mm9}
                        Reference genome
  -s {mouse,human}, --species {mouse,human}
                        Set species.
  -a {salmon,hisat2}, --align {salmon,hisat2}
                        Choose aligner. Default is Salmon.
  -p <P value>, --pvalue <P value>
                        Set P value cut off (default is 0.001)
  --go                  Gene set enrichment analysis with Enrichr
  --gene-sets <GO gene set>
                        Gene sets used for GO analysis (default is GO_Molecular_Function_2021, GO_Cellular_Component_2021, and
                        GO_Biological_Process_2021). Gene sets can be found on https://maayanlab.cloud/Enrichr/#stats
  --skip-fastqc         Skip FastQC/MultiQC

```

### ChIP-Seq analysis

```
pyseqtools.py chip-seq [-h] [-t THREADS] [-r] [-f] [-g GENOME] [-a {hisat2,bwa}] [-d] [-s] [-b] [-q] [-p] [-n]

optional arguments:
  -h, --help            show this help message and exit
  -t THREADS, --threads THREADS
                        <INT> number of CPU threads to use (default is 1). Use max to apply all available CPU threads
  -r, --rename          Rename fq files
  -f, --fastqc          Perform FASTQC
  -g GENOME, --genome GENOME
                        Choose reference genome (default is hg19)
  -a {hisat2,bwa}, --align {hisat2,bwa}
                        Trim and align raw data to index using HISAT2 or BWA
  -d, --deduplication   Perform deduplication of BAM files
  -s, --downsample      Perform downsampling of BAM files
  -b, --bigwig          Create BigWig files
  -q, --qc              Perform QC analysis of BAM files
  -p, --peaks           Call and annotate peaks
  -n, --ngsplot         Generate metageneplots and heatmaps with ngs.plot

```

### CUT&RUN analysis

under construction
