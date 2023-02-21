#!/bin/bash

##### pyseqtools Conda environment ####
conda create -n pyseqtools python=3.10
conda activate pyseqtools

#add Bioconda channel
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

#install conda pakcages
conda install samtools
conda install fastqc
conda install picard
conda install trim-galore
conda install mageck
conda install r-ngsplot
conda install salmon

#Python packages
pip install numpy matplotlib seaborn pandas tqdm clint pysam gseapy PyYAML multiqc deeptools



