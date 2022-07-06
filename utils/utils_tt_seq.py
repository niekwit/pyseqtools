#!/usr/bin/env python3

import glob
import os
import sys

import pysam

script_dir = os.path.abspath(os.path.dirname(__file__))
script_dir = os.path.dirname(script_dir)
sys.path.append(os.path.join(script_dir, "utils"))
import utils_general as utils


def splitBam(threads, work_dir):
    '''
    https://www.biostars.org/p/92935/
    
    '''
    print("Splitting BAM files into forward and reverse")
    
    file_list = glob.glob(os.path.join(work_dir, "bam", "*.bam"))
    
    
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
        
        