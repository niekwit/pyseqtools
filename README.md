# pyseqtools
Tool box for a variety of NGS data analyses

Available analyses/tools:

- CRISPR screens
- RNA-Seq
- ChIP-Seq
- CUT & RUN
- DamID (work in progress)


- Interspecies gene symbol conversion
- Subsetting of GTF file based on selected genes

## Installation

Perform the following command in the directory of your choice:
```
git clone https://github.com/niekwit/pyseqtools.git
```

## Software requirements

- Python >= 3.7
- Perl
- R

When running an analysis, `pyseqtools.py` will first check whether the relevant dependencies are available in `$PATH`. If not, it will then look for it in any `$HOME` directories, and finally if it does not find them there it will install them.

## Usage:

```
usage: pyseqtools.py [-h] {crispr,rna-seq,chip-seq,cutrun,damid,genesymconv,subsetgtf} ...

Analysis pipelines for a variety of NGS data, and related tools

positional arguments:
  {crispr,rna-seq,chip-seq,cutrun,damid,genesymconv,subsetgtf}
                        Available pipelines/tools:
    crispr              CRISPR screen analysis
    rna-seq             RNA-Seq analysis
    chip-seq            ChIP-Seq analysis
    cutrun              CUT & RUN analysis
    damid               DamID analysis
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
pyseqtools.py chip-seq [-h] [-t THREADS] [-r] [-g GENOME] [-a [{hisat2,bwa-mem,bwa-aln}]] [-d]
                              [--downsample] [-b] [--qc] [-p] [--metagene] [--skip-fastqc]

optional arguments:
  -h, --help            show this help message and exit
  -t THREADS, --threads THREADS
                        <INT> number of CPU threads to use (default is 1). Use max to apply all available
                        CPU threads
  -r, --rename          Rename fq files
  -g GENOME, --genome GENOME
                        Choose reference genome (default is hg19)
  -a [{hisat2,bwa-mem,bwa-aln}], --align [{hisat2,bwa-mem,bwa-aln}]
                        Create BAM files using HISAT2 or BWA (mem or aln
  -d, --deduplication   Perform deduplication of BAM files
  --downsample          Perform downsampling of BAM files
  -b, --bigwig          Create BigWig files
  --qc                  Perform QC analysis of BAM files
  -p, --peaks           Call and annotate peaks with MACS3/HOMER
  --metagene            Generate metageneplots and heatmaps with plotProfile (deepTools)
  --skip-fastqc         Skip FastQC/MultiQC


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
