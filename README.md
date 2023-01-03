# pyseqtools
Tool box for a variety of NGS data analyses

Available pipelines:

- [CRISPR screens](https://github.com/niekwit/pyseqtools#crispr-screen-analysis)
- [RNA-Seq](https://github.com/niekwit/pyseqtools#rna-seq-analysis)
- [ChIP-Seq](https://github.com/niekwit/pyseqtools#chip-seq-analysis)
- [CUT & RUN](https://github.com/niekwit/pyseqtools#cutrun-analysis)
- [DamID-Seq](https://github.com/niekwit/pyseqtools#damid-seq-analysis)
- [TT-Seq](https://github.com/niekwit/pyseqtools#tt-seq-analysis)
- [3end-Seq](https://github.com/niekwit/pyseqtools#3end-seq-analysis)


## Installation

Perform the following command in the directory of your choice:
```
git clone https://github.com/niekwit/pyseqtools.git
```
For your convenience you can add the `pyseqtools` directory to your `$PATH` by adding the following line to your `~/.bashrc` file:
```
export PATH=/home/path/to/pyseqtools:$PATH
```

Optional: enabe auto-completion of the different `pyseqtools` pipelines, genomes, and CRISPR libraries by adding the following line to your `~/.bashrc` file:
```
source /home/path/to/pyseqtools/bash/auto-complete.sh
```
## Software requirements

- Python >= 3.7
- Perl
- R

When running an analysis locally, `pyseqtools.py` will first check whether the relevant dependencies are available in `$PATH`. If not, it will then look for it in any `$HOME` directories, and finally if it does not find them there it will install them. When using the `SLURM` option however, each dependency is expected to be in `$PATH`

## Before analysis
1. Create analysis directory (parent directory). `pyseqtools` should be run from this parent directory. 
2. All fastq files should be located in a sub directory `raw-data`.
3. The configuration files `stats.config` (CRISPR screen analysis only), `rename.config` and `samples.csv` should be located in the parent directory.

## Usage:

```
usage: pyseqtools.py [-h] {crispr,rna-seq,chip-seq,cutrun,damid,tt-seq,3end-seq,genesymconv,subsetgtf} ...

Analysis pipelines for a variety of NGS data, and related tools

positional arguments:
  {crispr,rna-seq,chip-seq,cutrun,damid,tt-seq,3end-seq,genesymconv,subsetgtf}
                        Available pipelines/tools:
    crispr              CRISPR screen analysis
    rna-seq             RNA-Seq analysis
    chip-seq            ChIP-Seq analysis
    cutrun              CUT & RUN analysis
    damid               DamID analysis
    tt-seq              TT-Seq analysis according to https://github.com/crickbabs/DRB_TT-seq
    3end-seq            3end-Seq analysis
    genesymconv         Inter-species gene symbol conversion
    subsetgtf           Subset GTF files for selected genes

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
Before running `pyseqtools` prepare the following:
1. `stats.config` file.  
This contains the sample comparisons that will be fed into `MAGeCK` or `BAGEL2`:
```
t;c
D8SORT;D8LIB
D21SORT;D21LIB
D8SORT,D21SORT;D8LIB,D21LIB
```
	The first line is the header line (test;sample). Each following line contains a unique sample comparison: test sample name(s);control sample name(s)
	The sample names should be the same as the fastq files, but wihout the fastq.gz/fq.gz file extension. Multiple samples can be grouped when seperated with a comma (see example).
	In the `MAGeCK` output file scores for depletion and enrichment which will correspond to depletion/enrichment of sgRNAs against genes in the test sample compared to the control sample.

2. `rename.config` file.  
This contains information for renaming fastq files before analyis (optional):
```
SLX-20801.DNAA002.H55VLDRX2.s_2.r_1.fq.gz;D8LIB.fq.gz
SLX-20801.DNAA004.H55VLDRX2.s_2.r_1.fq.gz;D8SORT.fq.gz
SLX-20801.DNAA007.H55VLDRX2.s_2.r_1.fq.gz;D21LIB.fq.gz
SLX-20801.DNAA016.H55VLDRX2.s_2.r_1.fq.gz;D21SORT.fq.gz
```
Each line contains: original file name;new file name  

IMPORTANT NOTE: if your data is paired-end, only use the fastq files of mate 1 (read1).  

### RNA-Seq analysis

```
usage: pyseqtools.py rna-seq [-h] [-t <int>] [-g {hg38_76,mm10,hg38_50,hg19,hg38_100,mm9,hg38,gencode-v35}] [--trim]
                             [--peTags PETAGS] [-a {salmon,star}] [--rsemIndex RSEMINDEX [RSEMINDEX ...]]
                             [--isoformAnalysis {None,miso,rmats}] [--indexBAM] [--sortBAM] [-f] [--deseq2] [-b] [--TE]
                             [--go] [--gene-sets <GO gene set>] [-p <P value>] [--fastqc] [--slurm]

Analysis pipeline for RNA-Seq experiments

optional arguments:
  -h, --help            show this help message and exit
  -t <int>, --threads <int>
                        Number of CPU threads to use (default is 1). Use max to apply all available CPU threads. For
                        Salmon 8-12 threads are optimal
  -g {hg38_76,mm10,hg38_50,hg19,hg38_100,mm9,hg38,gencode-v35}, --genome {hg38_76,mm10,hg38_50,hg19,hg38_100,mm9,hg38,gencode-v35}
                        Reference genome
  --trim                Quality trimming of data using Trim_galore!
  --peTags PETAGS       Comma-separated paired-end file tags (e.g. _R1_001.fq.gz,_R2_001.fq.gz). Only required when
                        --trim argument is called.
  -a {salmon,star}, --align {salmon,star}
                        Program to align fastq files
  --rsemIndex RSEMINDEX [RSEMINDEX ...]
                        Generate index for STAR alignment via RSEM (--rsemIndex genome read-length out-dir
  --isoformAnalysis {None,miso,rmats}
                        Alternative isoform analysis with RSEM/MISO or rMATS
  --indexBAM            Index BAM files with samtools
  --sortBAM             Sort BAM files with samtools
  -f, --scaleFactors    Calculate scale factors from yeast spike-in RNA with DESeq2
  --deseq2              Calculate differential genes using DESeq2
  -b, --bigwig          Create BigWig files using bamCoverage (requires indexed BAM files
  --TE                  Transposable element expression analysis (Requires STAR alignment)
  --go                  Gene set enrichment analysis with Enrichr
  --gene-sets <GO gene set>
                        Gene sets used for GO analysis (default is GO_Molecular_Function_2021,
                        GO_Cellular_Component_2021, and GO_Biological_Process_2021). Gene sets can be found on
                        https://maayanlab.cloud/Enrichr/#stats
  -p <P value>, --pvalue <P value>
                        Set P value cut off (default is 0.001)
  --fastqc              Run FastQC/MultiQC
  --slurm               Submit jobs to Cambridge HPC using SLURM

```

### ChIP-Seq analysis

```
usage: pyseqtools.py chip-seq [-h] [-t THREADS] [-r] [-g GENOME] [--indexBAM] [--trim] [--peTags PETAGS]
                              [--align [{hisat2,bwa-mem,bwa-aln}]] [-d] [--downsample] [-s] [-b] [--qc] [-p] [--metagene]
                              [--slurm] [--fastqc]

Analysis pipeline for ChIP-Seq experiments

optional arguments:
  -h, --help            show this help message and exit
  -t THREADS, --threads THREADS
                        <INT> number of CPU threads to use (default is 1). Use max to apply all available CPU threads
  -r, --rename          Rename fq files
  -g GENOME, --genome GENOME
                        Choose reference genome (default is hg19)
  --indexBAM            Index BAM files with samtools
  --trim                Quality trimming of data using Trim_galore!
  --peTags PETAGS       Comma-separated paired-end file tags (e.g. _R1_001.fq.gz,_R2_001.fq.gz). Only required when
                        --trim argument is called.
  --align [{hisat2,bwa-mem,bwa-aln}]
                        Create BAM files using HISAT2 or BWA (mem or aln)
  -d, --deduplication   Perform deduplication of BAM files
  --downsample          Perform downsampling of BAM files to lowest read count
  -s, --spike-in        Perform normalisation based on Drosophila chromatin spike-in
  -b, --bigwig          Create BigWig files
  --qc                  Perform QC analysis of BAM files
  -p, --peaks           Call and annotate peaks with MACS3/HOMER
  --metagene            Generate metageneplots and heatmaps with plotProfile (deepTools)
  --slurm               Submit jobs to Cambridge HPC using SLURM
  --fastqc              Run FastQC/MultiQC

```

### CUT&RUN analysis

```
pyseqtools.py cutrun [-h] [-t THREADS] [-g GENOME] [-a [{bowtie,bowtie2,hisat2,bwa-mem,bwa-aln}]]
                            [-d] [-b] [--qc] [-p] [--metagene] [--skip-fastqc]

Analysis pipeline for CUT & RUN experiments

optional arguments:
  -h, --help            show this help message and exit
  -t THREADS, --threads THREADS
                        <INT> number of CPU threads to use (default is 1). Use max to apply all available
                        CPU threads
  -g GENOME, --genome GENOME
                        Choose reference genome (default is hg19)
  -a [{bowtie,bowtie2,hisat2,bwa-mem,bwa-aln}], --align [{bowtie,bowtie2,hisat2,bwa-mem,bwa-aln}]
                        Create BAM files using Bowtie2, HISAT2 or BWA (mem or aln
  -d, --deduplication   Perform deduplication of BAM files
  -b, --bigwig          Create BigWig files
  --qc                  Perform QC analysis of BAM files
  -p, --peaks           Call and annotate peaks with MACS3/HOMER
  --metagene            Generate metageneplots and heatmaps with plotProfile (deepTools)
  --skip-fastqc         Skip FastQC/MultiQC

```

### DamID-Seq analysis
```
pyseqtools.py damid [-h] [-t THREADS] [-g GENOME] [--skip-fastqc]

Analysis pipeline for DamID, based on https://github.com/owenjm/damidseq_pipeline

optional arguments:
  -h, --help            show this help message and exit
  -t THREADS, --threads THREADS
                        <INT> number of CPU threads to use (default is 1). Use max to apply all available
                        CPU threads
  -g GENOME, --genome GENOME
                        Choose reference genome (default is hg19)
  --skip-fastqc         Skip FastQC/MultiQC

```

### TT-Seq analysis
```
usage: pyseqtools.py tt-seq [-h] [-t THREADS] [-r] [-a] [-g [GENOME]] [--splitBAM] [-f] [-b] [--deseq2] [--metagene]
                            [--readRatio] [--slurm]

Analysis pipeline for TT-Seq experiments

optional arguments:
  -h, --help            show this help message and exit
  -t THREADS, --threads THREADS
                        <INT> number of CPU threads to use (default is 1). Use max to apply all available CPU threads
  -r, --rename          Rename fastq files
  -a, --align           Alignment of fastq files using STAR
  -g [GENOME], --genome [GENOME]
                        Choose reference genome (default is hg38)
  --splitBAM            Generate forward and reverse strand specific BAM files
  -f, --sizeFactors     Calculate size factors from yeast spike-in RNA with DESeq2
  -b, --bigwig          Create BigWig files
  --deseq2              Calculate differential genes using DESeq2
  --metagene            Generate metagene plots and heatmaps with ngs.plot
  --readRatio           Calculate ratios of read numbers in areas around TSS vs TES
  --slurm               Submit jobs to Cambridge HPC using SLURM
```

### 3end-Seq analysis

Under construction

### Gene symbol conversion
```
pyseqtools.py genesymconv [-h] -c {hm,mh} -i INPUT -o OUTPUT

Convert human gene symbols to mouse gene symbols, or vice versa

optional arguments:
  -h, --help            show this help message and exit
  -c {hm,mh}, --conversion {hm,mh}
                        Conversion type: human to mouse (hm) or mouse to human (mh)
  -i INPUT, --input INPUT
                        Input gene list file. Each gene symbol should be on a new line
  -o OUTPUT, --output OUTPUT
                        Output file name

```

### Generate subset of GTF based on selected genes
```
pyseqtools.py subsetgtf [-h] -l LIST -i INPUT -o OUTPUT

Tool for subsetting GTF files according to input gene list

optional arguments:
  -h, --help            show this help message and exit
  -l LIST, --list LIST  Input gene list file. Each gene symbol should be on a new line
  -i INPUT, --input INPUT
                        GTF file to be subsetted (must be specified with genome name, loaded from chip-
                        seq.yaml)
  -o OUTPUT, --output OUTPUT
                        Output file name
```
