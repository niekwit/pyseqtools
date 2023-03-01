#!/usr/bin/env python3

'''Merges replicate fq.gz files:
    
'''

import os
import subprocess
import sys
import glob


args = sys.argv

#get parsed arguments
work_dir = args[1]

#get bam files for read1/2
samples_r1 = glob.glob(os.path.join(work_dir,"raw-data","*.r_1.fq.gz"))
samples_r1 = [x.replace(".r_1.fq.gz","") for x in samples_r1]

samples_r2 = glob.glob(os.path.join(work_dir,"raw-data","*.r_2.fq.gz"))
samples_r2 = [x.replace(".r_2.fq.gz","") for x in samples_r2]

#get base sample names
base_samples = list(set([x.rsplit(".",1)[0] for x in samples_r1]))

#function to merge fq files
def merge(base,extension):
    fq = glob.glob(os.path.join(work_dir,"raw-data",base + "*" + extension))
    fq.sort()
    fq = " ".join(fq)
    merged_fq = base + extension
    
    #create cat command
    command = f"cat {fq} > {merged_fq}"
    
    #run command
    subprocess.call(command,shell=True)

#merge replicate read1/read2 samples
for i in base_samples:
        #merge read1 samples:
        merge(i,".r_1.fq.gz")
            
        #merge read2 samples:
        merge(i,".r_2.fq.gz")
