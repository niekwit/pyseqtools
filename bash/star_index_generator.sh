#!/usr/bin/env bash



STAR --runThreadN 44 --runMode genomeGenerate --genomeDir $PWD --genomeFastaFiles /home/niek/Documents/scripts/pyseqtools/fasta/hg38/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile /home/niek/Documents/scripts/scripts_jonnie/Homo_sapiens.GRCh38.107.gtf --sjdbOverhang 49
