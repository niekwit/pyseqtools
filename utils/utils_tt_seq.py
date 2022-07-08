#!/usr/bin/env python3

import glob
import os
import sys
import subprocess
import shutil
import collections


from clint.textui import colored, puts
import pysam
import HTSeq

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
            
            puts(colored.green("Aligning " +sample + " to " + genome))
            if not utils.file_exists(sorted_bam):
                star = ["STAR", "--runThreadN", threads,"--runMode", "alignReads", "--genomeDir", index,
                        "--readFilesIn", read1, read2, "--readFilesCommand", "zcat", "--quantMode",
                        "TranscriptomeSAM", "GeneCounts", "--twopassMode", "Basic", "--outSAMunmapped",
                        "None", "--outSAMattrRGline","ID:"+sample,"PU:"+sample,"SM:"+sample,"LB:unknown",
                        "PL:illumina", "--outSAMtype","BAM", "Unsorted", "--outTmpDir", temp_dir,
                        "--outFileNamePrefix", os.path.join(work_dir,"bam",genome,sample,sample)]
                
                utils.write2log(work_dir, " ".join(star), "" )
                subprocess.call(star)

            #sort bam file
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
        
               
def splitBam(threads, work_dir, genome):
    '''
    https://www.biostars.org/p/92935/
    
    '''
    puts(colored.green("Creating forward and reverse strand BAM files"))

    file_list = glob.glob(os.path.join(work_dir, "bam", genome, "*", "*_sorted.bam"))
    
        
    for bam in file_list:
        ##forward strand
        fwd1 = bam.replace(".bam","_fwd1.bam")
        fwd2 = bam.replace(".bam","_fwd2.bam")
        
        #alignments of the second in pair if they map to the forward strand
        pysam.view("-@",threads,"-b","-f","128","-F","16",bam,"-o",fwd1, catch_stdout=False)
        pysam.index(fwd1)
        
        #alignments of the first in pair if they map to the reverse  strand
        pysam.view("-@",threads,"-b","-f","80",bam,"-o",fwd2, catch_stdout=False)
        pysam.index(fwd2)
        
        #merge all forward reads
        fwd = bam.replace(".bam","_fwd.bam")
        pysam.merge("-f",fwd,fwd1,fwd2, catch_stdout=False)
        pysam.index(fwd)
        
        ##reverse strand
        rev1 = bam.replace(".bam","_rev1.bam")
        rev2 = bam.replace(".bam","_rev2.bam")
        
        #alignments of the second in pair if they map to the reverse strand
        pysam.view("-b","-f","144",bam,"-o", rev1, catch_stdout=False)
        pysam.index(rev1)
        
        #alignments of the first in pair if they map to the forward strand
        pysam.view("-@",threads,"-b","-f","64","-F","16",bam,"-o",rev2, catch_stdout=False)
        pysam.index(rev2)
        
        #merge all reverse reads
        rev = bam.replace(".bam","_rev.bam")
        pysam.merge("-f",rev,rev1,rev2, catch_stdout=False)
        pysam.index(fwd)
        
        #remove all non-merged bam files
        remove = [fwd1,fwd2,rev1,rev2]
        for file in remove:
            os.remove(file)
            os.remove(file.replace(".bam",".bam.bai"))
            
            
def sizeFactors(work_dir, gtf):
    #first prepare htseq-count input for DESeq2
    file_list = os.path.join(work_dir,"bam","R64-1-1","*","*_sorted.bam")
        
    for bam in file_list:
        rc_file = bam.replace("_sorted.bam","") + "_count.txt"
        
        
        if not utils.file_exists(rc_file):
           htseq_count = ["htseq-count", "--format", "bam", "--stranded", "yes", bam, gtf]
           utils.write2log(work_dir, " ".join(htseq_count), "" )
           subprocess.call(htseq_count)
    
    
    
    
    
    
    