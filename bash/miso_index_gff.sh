#!/bin/bash

#SBATCH -A JNATHAN-SL3-CPU
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH -p cclake
#SBATCH -D /home/nw416/rds/hpc-work/gff3/hg38/
#SBATCH -o gff_index.log
#SBATCH -c 12
#SBATCH -t 01:00:00
#SBATCH --mem=20G
#SBATCH -J miso-index

#open conda environment that contains MISO installation (Python=2.7.6+)
source ~/.bashrc
conda activate miso

#prepare GFF3 file for hg38
wget -o /home/nw416/rds/hpc-work/gff3/hg38/wgEncodeGencodeCompV28.txt.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/wgEncodeGencodeCompV28.txt.gz
gunzip /home/nw416/rds/hpc-work/gff3/hg38/wgEncodeGencodeCompV28.txt.gz
genePredToGtf 


index_gff --index /home/nw416/rds/hpc-work/gff3/hg38/Homo_sapiens.GRCh38.107.gff3 /home/nw416/rds/hpc-work/gff3/hg38/miso-index

