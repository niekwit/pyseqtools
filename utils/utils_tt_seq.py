#!/usr/bin/env python3

import glob
import os
import sys
import subprocess
import shutil
import pandas as pd

from clint.textui import colored, puts
import pysam

script_dir = os.path.abspath(os.path.dirname(__file__))
script_dir = os.path.dirname(script_dir)
sys.path.append(os.path.join(script_dir, "utils"))
import utils_general as utils


        
def STAR(work_dir, threads, script_dir, tt_seq_settings, genome):
    '''
    based on https://github.com/crickbabs/DRB_TT-seq/
    '''

    file_list = glob.glob(os.path.join(work_dir,"trim","*_val_1.fq.gz"))

    
    #function for alignment with STAR
    def align(work_dir,file_list, index, threads, genome):
        for read1 in file_list:
            read2 = read1.replace("_R1_001_val_1.fq.gz","_R2_001_val_2.fq.gz")
            
            #create sample name
            sample = os.path.basename(read1).replace("_R1_001_val_1.fq.gz","")
            
            #create temp dir path for STAR and make sure it does not exist
            temp_dir = os.path.join(work_dir,"tmp")
            if os.path.exists(temp_dir):
                shutil.rmtree(temp_dir)
            
            #create output dir
            os.makedirs(os.path.join(work_dir,"bam",genome,sample), exist_ok = True)
            
            bam = os.path.join(work_dir,"bam",genome,sample,sample+"Aligned.out.bam")
            sorted_bam= bam.replace("Aligned.out.bam","_sorted.bam")
            
            #create STAR command
            star = ["STAR", "--runThreadN", threads,"--runMode", "alignReads", "--genomeDir", index,
                    "--readFilesIn", read1, read2, "--readFilesCommand", "zcat", "--quantMode",
                    "TranscriptomeSAM", "GeneCounts", "--twopassMode", "Basic", "--outSAMunmapped",
                    "None", "--outSAMattrRGline","ID:"+sample,"PU:"+sample,"SM:"+sample,"LB:unknown",
                    "PL:illumina", "--outSAMtype","BAM", "Unsorted", "--outTmpDir", temp_dir,
                    "--outFileNamePrefix", os.path.join(work_dir,"bam",genome,sample,sample)]
            
            #run STAR
            puts(colored.green(f"Aligning {sample} to {genome}"))
            if genome == "hg38": 
                if not utils.file_exists(sorted_bam):
                    utils.write2log(work_dir, " ".join(star), "" )
                    subprocess.call(star)
            elif genome == "R64-1-1": #it is better for HTSeq to use unsorted BAM files
                if not utils.file_exists(bam):
                    utils.write2log(work_dir, " ".join(star), "" )
                    subprocess.call(star)

            #only sort hg38 bam files
            if genome == "hg38":
                puts(colored.green("Sorting " + os.path.basename(bam)))
                
                if not utils.file_exists(sorted_bam):
                    pysam.sort("--threads", threads,"-o",sorted_bam,bam)
                
                #remove unsorted bam file
                if os.path.exists(sorted_bam):
                    if os.path.exists(bam):
                        os.remove(bam)
             
    #align trimmed reads to selected genome    
    index = tt_seq_settings["index"][genome]
    align(work_dir,file_list, index, threads,genome)
    
    #align trimmed reads to yeast genome (spike-in)
    yeast_index = tt_seq_settings["index-yeast"]
    align(work_dir,file_list, yeast_index, threads,"R64-1-1")
    
    #index all bam files
    utils.indexBam(work_dir, threads, genome)
    utils.indexBam(work_dir, threads, "R64-1-1")
        

def splitBam(threads, work_dir, genome):
    '''
    based on https://www.biostars.org/p/92935/
    
    '''
    puts(colored.green("Creating forward and reverse strand-specific BAM files"))

    file_list = glob.glob(os.path.join(work_dir, "bam", genome, "*", "**_sorted.bam"))
    
        
    for bam in file_list:
        print(os.path.basename(bam))
        ##forward strand
        fwd1 = bam.replace("*_sorted.bam","_fwd1.bam")
        fwd2 = bam.replace("*_sorted.bam","_fwd2.bam")
        
        #alignments of the second in pair if they map to the forward strand
        pysam.view("-@",threads,"-b","-f","128","-F","16",bam,"-o",fwd1, catch_stdout=False)
        pysam.index(fwd1)
        
        #alignments of the first in pair if they map to the reverse  strand
        pysam.view("-@",threads,"-b","-f","80",bam,"-o",fwd2, catch_stdout=False)
        pysam.index(fwd2)
        
        #merge all forward reads
        fwd = bam.replace("*_sorted.bam","_fwd.bam")
        pysam.merge("-f",fwd,fwd1,fwd2, catch_stdout=False)
        pysam.index(fwd)
        
        ##reverse strand
        rev1 = bam.replace("*_sorted.bam","_rev1.bam")
        rev2 = bam.replace("*_sorted.bam","_rev2.bam")
        
        #alignments of the second in pair if they map to the reverse strand
        pysam.view("-b","-f","144",bam,"-o", rev1, catch_stdout=False)
        pysam.index(rev1)
        
        #alignments of the first in pair if they map to the forward strand
        pysam.view("-@",threads,"-b","-f","64","-F","16",bam,"-o",rev2, catch_stdout=False)
        pysam.index(rev2)
        
        #merge all reverse reads
        rev = bam.replace("*_sorted.bam","_rev.bam")
        pysam.merge("-f",rev,rev1,rev2, catch_stdout=False)
        pysam.index(fwd)
        
        #remove all non-merged bam files
        remove = [fwd1,fwd2,rev1,rev2]
        for file in remove:
            os.remove(file)
            os.remove(file.replace(".bam",".bam.bai"))
            
            
def sizeFactors(work_dir, tt_seq_settings):
    """
    Based on https://github.com/crickbabs/DRB_TT-seq/blob/master/bigwig.md with modifications

    """
    puts(colored.green("Generating size factors for normalisation using DESeq2"))
    #first prepare htseq-count input for DESeq2
    print("Quantifying reads with HTSeq for DESeq2 input")
    file_list = glob.glob(os.path.join(work_dir,"bam","R64-1-1","*","*_sorted.bam"))
    rc_file = [os.path.join(work_dir, "htseq-count", os.path.basename(x.replace("_sorted.bam","_count.txt"))) for x in file_list] 
    
    #create output dir
    os.makedirs(os.path.join(work_dir,"htseq-count"), exist_ok = True)
    
    #load yeast gtf file
    gtf = tt_seq_settings["gtf-yeast"]
    
    for bam,count_file in zip(file_list,rc_file):
                
        if not utils.file_exists(count_file):
           htseq_count = ["htseq-count", "--quiet", "--format", "bam", "--stranded", "yes", bam, gtf, ">", count_file]
           #utils.write2log(work_dir, " ".join(htseq_count), "" )
           print("Quantifying " + os.path.basename(bam) + " with htseq-count")
           subprocess.call(htseq_count)
   
    #run DESeq2 to obtain size factors for normalisation
    print("Creating size factors with DESeq2 with HTSeq files")
    deseq2 = ["Rscript", os.path.join(script_dir, "R", "tt-seq_sizefactors.R")]
    subprocess.call(deseq2)
    
    
def ttSeqBigWig(work_dir, threads, tt_seq_settings, genome):
    """
    Create BigWig files for TT-Seq with scaling factors derived from DESeq2

    """
    
    
    #load scaling factors
    try:
        scaling_factors = pd.read_csv(os.path.join(work_dir, "sizeFactors.csv"))
    except FileNotFoundError:
        return(puts(colored.red("ERROR: sizeFactors.csv not found")))
    
    #create BigWig directory
    os.makedirs(os.path.join(work_dir,"bigwig"), exist_ok = True)
    
    #create scaled BigWig files
    for index,row in scaling_factors.iterrows():
        bam = os.path.join(work_dir,"bam", genome, 
                           row["sample"].replace("_count.txt","")
                           ,row["sample"].replace("_count.txt","_sorted.bam"))
        scaling_factor = row["sizeFactors"]
        
        #BigWig function
        def bigWig(work_dir, threads, bam, strand):
            
            fw_output = os.path.join(work_dir, "bigwig",os.path.basename(bam).replace("_sorted.bam","_" + strand +".bigwig"))
            bigwig = ["bamCoverage","--scaleFactor",scaling_factor, "-p", threads, "-b", bam, "-o", fw_output]
            puts(colored.green(f"Generating {strand} BigWig file for {os.path.basename(bam)}"))
            if not utils.file_exists(fw_output):
                utils.write2log(work_dir, " ".join(bigwig))
                subprocess.call(bigwig)
        
        strands = ["fwd", "rev"]
        for i in strands:
            bigWig(work_dir, threads, bam, i)
    
    
    #create mean Wig files for all technical replicates with wiggletools
    samples = pd.read_csv(os.path.join(work_dir,"samples.csv"))
    genotypes = set(samples["genotype"])
    conditions = set(samples["condition"])
    
    for genotype in genotypes:
        for condition in conditions:
            for strand in strands:
                sub_samples = samples[samples["genotype"] == genotype]
                sub_samples = sub_samples[samples["condition"] == condition]
                sub_samples = list(sub_samples["sample"])
                in_bigwigs = [os.path.join(work_dir,"bigwig",x + "_" + strand +".bigwig") for x in sub_samples]
                
                wig_mean = os.path.join(work_dir,"bigwig",genotype+"_"+condition+"_mean.wig")
                wiggletools = ["wiggletools", "write", wig_mean, "mean", " ".join(in_bigwigs)]
                utils.write2log(work_dir, " ".join(wiggletools))
                subprocess.call(wiggletools)
    
    #generate chrom.sizes file needed for converting Wig to BigWig
    fasta = tt_seq_settings["fasta"][genome]
    
    chrom_sizes = os.path.join(script_dir,"chrom.sizes",os.path.basename(fasta))
    chrom_sizes = chrom_sizes.rsplit(".",1)[0] + ".chrom.sizes"
    
    if not utils.file_exists(chrom_sizes):
        fasta_index = fasta + ".fai"
        if not utils.file_exists(fasta_index):
            pysam.faidx(fasta)
            bash = ["cut", "-f1,2", fasta_index, ">", chrom_sizes]
            utils.write2log(work_dir, " ".join(bash),"Generate chrom.sizes file: ")
            subprocess.call(bash)
    
    #convert Wig to BigWig
    file_list = glob.glob(os.path.join(work_dir,"bigwig","*_mean.wig"))
    
    
    for wig in file_list:
        bw = wig.replace("*_mean.wig","_mean.bigwig")
        if not utils.file_exists(bw):
            wig2bigwig = ["wigToBigWig", wig, chrom_sizes, bw]
            utils.write2log(work_dir,wig2bigwig)
            subprocess.call(wig2bigwig)


def bwQC(work_dir, threads):
    """
    Quality control for BigWig files using deepTools
    """
    #samples = pd.read_csv(os.path.join(work_dir,"samples.csv"))
    bigwig_all = glob.glob(os.path.join(work_dir,"bigwig","*[!_mean].bigwig"))
    bigwig_mean = glob.glob(os.path.join(work_dir,"bigwig","*_mean.bigwig"))
    
    

def metaProfiles(work_dir, threads, tt_seq_settings, genome):
    #subset GTF for protein coding genes
    gtf = tt_seq_settings["gtf"][genome]
    gtf_pc = gtf.replace(".gtf", "_pc.gtf")
    
    if not utils.file_exists(gtf_pc):
        grep = ["grep", "protein_coding", gtf, ">", gtf_pc]
        utils.write2log(work_dir, " ".join(grep), "Subset GTF file for only protein coding genes: ")
        subprocess.call(grep)
        
    #create compute matrix with deepTools
    matrix = os.path.join(work_dir, "bigwig", )


def trxReadThrough(work_dir, threads):
    pass