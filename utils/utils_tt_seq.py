#!/usr/bin/env python3

import glob
import os
import sys
import subprocess
import shutil
import pandas as pd
import yaml
from pathlib import Path
import tempfile

from clint.textui import colored, puts
import pysam
#import pybedtools
#import biomart

script_dir = os.path.abspath(os.path.dirname(__file__))
script_dir = os.path.dirname(script_dir)
sys.path.append(os.path.join(script_dir, "utils"))
import utils_general as utils

        
def STAR(work_dir, threads, script_dir, tt_seq_settings, genome, slurm=False, job_id_trim=None):
    '''
    Alignment of trimmed fastq files with STAR.
    Based on https://github.com/crickbabs/DRB_TT-seq/
    If running on the HPC, this function will generate the SLURM script for parallel alignments
    '''
    
    #function for alignment with STAR
    def align(work_dir, index, threads, genome, slurm=False):
        #get list of trimmed fastq files
        if slurm == False:
            file_list = glob.glob(os.path.join(work_dir,"trim","*_val_1.fq.gz"))
        else:
            #create trim output file list
            extension = utils.get_extension(work_dir)
            file_list = glob.glob(os.path.join(work_dir, "raw-data","*R1_001." + extension))
            file_list = [x.split(".",1)[0] + "_val_1.fq.gz" for x in file_list]
            file_list = [x.replace("raw-data", "trim") for x in file_list]
        
            #create empty csv file for commands
            Path(os.path.join(work_dir,"slurm",f"slurm_STAR_{genome}.csv")).touch()
        
        for read1 in file_list:
            read2 = read1.replace("_R1_001_val_1.fq.gz","_R2_001_val_2.fq.gz")
            
            #create sample name
            sample = os.path.basename(read1).replace("_R1_001_val_1.fq.gz","")
            
            #create temp dir path for STAR and make sure it does not exist
            if slurm == False:
                temp_dir = os.path.join(work_dir,"tmp")
                if os.path.exists(temp_dir):
                    shutil.rmtree(temp_dir)
            else: 
                #each dir should have a unique name otherwise parallel alignments cannot be run
                temp_dir = tempfile.mkdtemp()
                shutil.rmtree(temp_dir)
            
            #create output dir
            os.makedirs(os.path.join(work_dir,"bam",genome,sample), exist_ok = True)
            
            bam = os.path.join(work_dir,"bam",genome,sample,sample+"Aligned.out.bam")
            sorted_bam= bam.replace("Aligned.out.bam","_sorted.bam")
            
            #if running on cluster, get thread count from cluster.yaml
            if slurm == True:
                with open(os.path.join(script_dir,"yaml","slurm.yaml")) as file:
                    slurm_settings = yaml.full_load(file)
                threads = str(slurm_settings["TT-Seq"]["STAR_CPU"])
            
            #create STAR command
            star = ["STAR", "--runThreadN", threads,"--runMode", "alignReads", "--genomeDir", index,
                    "--readFilesIn", read1, read2, "--readFilesCommand", "zcat", "--quantMode",
                    "TranscriptomeSAM", "GeneCounts", "--twopassMode", "Basic", "--outSAMunmapped",
                    "None", "--outSAMattrRGline","ID:"+sample,"PU:"+sample,"SM:"+sample,"LB:unknown",
                    "PL:illumina", "--outSAMtype","BAM", "Unsorted", "--outTmpDir", temp_dir,
                    "--outFileNamePrefix", os.path.join(work_dir,"bam",genome,sample,sample)]
            
            if slurm == False:
                #run STAR locally
                puts(colored.green(f"Aligning {sample} to {genome}"))
                if not utils.file_exists(sorted_bam):
                    utils.write2log(work_dir, " ".join(star), "" )
                    subprocess.call(star)
                    #only sort non-R64-1-1 bam files (it is better for HTSeq to use unsorted BAM files) 
                    if genome != "R64-1-1":
                        puts(colored.green("Sorting " + os.path.basename(bam)))
                        if not utils.file_exists(sorted_bam):
                            pysam.sort("--threads", threads,"-o",sorted_bam,bam)
                        
                        #remove unsorted bam file
                        if os.path.exists(sorted_bam):
                            if os.path.exists(bam):
                                os.remove(bam)
            else:
                #create csv files with STAR commands for SLURM job
                if not utils.file_exists(bam):
                    csv = open(os.path.join(work_dir,"slurm",f"slurm_STAR_{genome}.csv"), "a")  
                    csv.write(" ".join(star) +"\n")
                    csv.close()   
                    
        #index all bam files
        if slurm == False:
            utils.indexBam(work_dir, threads, genome)
            utils.indexBam(work_dir, threads, "R64-1-1")
        
        if slurm == True:
            #load slurm settings  
            mem = str(slurm_settings["TT-Seq"]["STAR_mem"])
            slurm_time = str(slurm_settings["TT-Seq"]["STAR_time"])
            account = slurm_settings["groupname"]
            partition = slurm_settings["TT-Seq"]["partition"]
            
            #create slurm bash script for splitting bam files
            print(f"Generating slurm_STAR_{genome}.sh")
            csv = os.path.join(work_dir,"slurm",f"slurm_STAR_{genome}.csv")
            commands = int(subprocess.check_output(f"cat {csv} | wc -l", shell = True).decode("utf-8"))
            script_ = os.path.join(work_dir,"slurm",f"slurm_STAR_{genome}.sh")
            script = open(script_, "w")  
            script.write("#!/bin/bash" + "\n")
            script.write("\n")
            script.write("#SBATCH -A " + account + "\n")
            script.write("#SBATCH --mail-type=BEGIN,FAIL,END" + "\n")
            script.write("#SBATCH -p " + partition + "\n")
            script.write("#SBATCH -D " + work_dir + "\n")
            script.write(f"#SBATCH -o slurm/slurm_STAR_%a_{genome}.log" + "\n")
            script.write("#SBATCH -c " + threads + "\n")
            script.write("#SBATCH -t " + slurm_time + "\n")
            script.write("#SBATCH --mem=" + mem + "\n")
            script.write("#SBATCH -J " + "STAR_"+genome + "\n")
            script.write("#SBATCH -a " + "1-" + str(commands) + "\n")
            script.write("\n")
            script.write("sed -n ${SLURM_ARRAY_TASK_ID}p " + csv +" | bash\n")
            script.close()
                    
            #run slurm script
            if job_id_trim is None:
                script = os.path.join(work_dir,"slurm",f"slurm_STAR_{genome}.sh")
                print("Submitting slurm script to cluster")
                job_id_align = subprocess.check_output(f"sbatch {script} | cut -d ' ' -f 4", shell = True)
                return(job_id_align)
            else:
                print("Submitting slurm script to cluster")
                job_id_align = subprocess.check_output(f"sbatch --dependency=afterok:{job_id_trim} {script} | cut -d ' ' -f 4", shell = True)  
                return(job_id_align)
            
            #sort BAM files
            bamSortSLURM(work_dir, job_id_align, genome)
              
    
    def alignPerformedCluster(work_dir,genome):
        '''
        To check if alignment has already been done on cluster
        '''
        #get number of output bam files
        if genome != "R64-1-1":
            bam_number = len(glob.glob(os.path.join(work_dir,"bam",genome, "*", "*_sorted.bam")))
        else:
            bam_number = len(glob.glob(os.path.join(work_dir,"bam","R64-1-1", "*", "*Aligned.out.bam")))
        
        #get number of align commands
        csv = os.path.join(work_dir,"slurm",f"slurm_STAR_{genome}.csv")
        if not os.path.isfile(csv): #if the slurm csv file does not exist, then no alignments have previously been run at all 
            return(False)
        
        commands = int(subprocess.check_output(f"cat {csv} | wc -l", shell = True).decode("utf-8"))
        
        #check if number of bam files is the same as align commands number:
        if bam_number == commands:
            return(True)
        else:
            return(False)
                
    #align trimmed reads to selected genome
    index = tt_seq_settings["STAR"][genome]
    puts(colored.green(f"Aligning fastq files to {genome} with STAR"))
    
    if slurm == True:
        if alignPerformedCluster(work_dir, genome) == False:
            align(work_dir, index, threads, genome, slurm)
        else:
            print(f"Skipping STAR alignment, all output BAM files for {genome} already detected")
    else:
        align(work_dir, index, threads, genome)
        
       
    #align trimmed reads to yeast genome (spike-in)
    puts(colored.green("Aligning fastq files to R64-1-1 (spike-in) with STAR"))
    yeast_index = tt_seq_settings["STAR"]["yeast"]
    
    if slurm == True:
        if alignPerformedCluster(work_dir, "R64-1-1") == False:
            job_id_align = align(work_dir, yeast_index, threads,"R64-1-1", slurm)
            return(job_id_align)
        else:
            print("Skipping STAR alignment, all output BAM files for R64-1-1 already detected")
    else:
        align(work_dir, yeast_index, threads, "R64-1-1")
    
    
    

def bamSortSLURM(work_dir, job_id_align, genome="hg38"):
    puts(colored.green("Sorting BAM files on cluster"))
    #create BAM file list (files might not exist yet)
    extension = utils.get_extension(work_dir)
    file_list = glob.glob(os.path.join(work_dir, "raw-data","*R1_001." + extension))
    file_list = [os.path.basename(x.split(".",1)[0].replace("_R1_001","")) for x in file_list]
    file_list = [os.path.join(work_dir,"bam",genome,x, x  + "Aligned.out.bam") for x in file_list]
    file_list_index = [x + ".bai" for x in file_list]
    
    #load SLURM settings
    with open(os.path.join(script_dir,"yaml","slurm.yaml")) as file:
        slurm_settings = yaml.full_load(file)
    threads = str(slurm_settings["samtools-sort_CPU"])
    
    mem = str(slurm_settings["samtools-sort_mem"])
    slurm_time = str(slurm_settings["samtools-sort_time"])
    account = slurm_settings["groupname"]
    partition = slurm_settings["partition"]
    
    #create csv file with samtools index commands
    csv = os.path.join(work_dir,"slurm","slurm_sortBAM.csv")
    if os.path.exists(csv):
        os.remove(csv)
    
    for bam, index in zip(file_list,file_list_index):
        if not utils.file_exists(index):
            samtools = ["samtools","sort","-@", threads, "-o", index, bam]
            csv = open(os.path.join(work_dir,"slurm","slurm_sortBAM.csv"), "a")  
            csv.write(" ".join(samtools) +"\n")
            csv.close() 
    
    #create slurm bash script 
    print("Generating slurm_sortBAM.sh")
    csv = os.path.join(work_dir,"slurm","slurm_sortBAM.csv")
    commands = int(subprocess.check_output(f"cat {csv} | wc -l", shell = True).decode("utf-8"))
    script_ = os.path.join(work_dir,"slurm","slurm_sortBAM.sh")
    script = open(script_, "w")  
    script.write("#!/bin/bash" + "\n")
    script.write("\n")
    script.write("#SBATCH -A " + account + "\n")
    script.write("#SBATCH --mail-type=BEGIN,FAIL,END" + "\n")
    script.write("#SBATCH -p " + partition + "\n")
    script.write("#SBATCH -D " + work_dir + "\n")
    script.write("#SBATCH -o slurm/slurm_sortBAM_%a.log" + "\n")
    script.write("#SBATCH -c " + threads + "\n")
    script.write("#SBATCH -t " + slurm_time + "\n")
    script.write("#SBATCH --mem=" + mem + "\n")
    script.write("#SBATCH -J " + "sortBAM" +"\n")
    script.write("#SBATCH -a " + "1-" + str(commands) + "\n")
    script.write("\n")
    script.write("sed -n ${SLURM_ARRAY_TASK_ID}p " + csv +" | bash\n")
    script.close()
   
    #submit SLURM bash script to cluster
    print("Submitting slurm_sortBAM.sh to cluster")
    script_ = os.path.join(work_dir,"slurm","slurm_sortBAM.sh")
    subprocess.call(["sbatch", f"--dependency=afterok:{job_id_align}", script])
    
    
def hisat2(work_dir, threads, tt_seq_settings, genome, slurm=False, job_id_trim=None):
    """
    Align TT-Seq quality trimmed paired-end files to genome with HISAT2
    """
    #function for aligning
    def align(work_dir,threads, tt_seq_settings,genome, slurm=False, job_id_trim=None):
        if slurm is False:
            read1_list = glob.glob(os.path.join(work_dir,"trim","*_R1_001_val_1.fq.gz"))
            hisat2_index = tt_seq_settings["hisat2"][genome]
            os.makedirs(os.path.join(work_dir,"bam", genome), exist_ok=True)
            
            for read1 in read1_list:
                read2 = read1.replace("_R1_001_val_1.fq.gz","_R2_001_val_2.fq.gz")
                bam = os.path.join(work_dir,"bam", genome,os.path.basename(read1.replace("_R1_001_val_1.fq.gz","_sort.bam")))
                sample = os.path.basename(read1.replace("_R1_001_val_1.fq.gz",""))
                puts(colored.green(f"Aligning {sample} to {genome} with HISAT2"))
                
                if not utils.file_exists(bam):
                    hisat2 = ["hisat2", "-p", str(threads), "-x", hisat2_index,"-1", read1, "-2", read2,
                              "2>>", os.path.join(work_dir,"align.log"), "|", "samtools", "view", "-q", "15", "-F", 
                              "260", "-b", "-@", str(threads), "-"]
                    if genome != "R64-1-1":
                        extend_hisat2 = ["|", "samtools", "sort", "-@", str(threads), "-", ">", bam]
                        hisat2.extend(extend_hisat2)
                    else:
                        #don't sort yeast bam files
                        bam = bam.replace("_sort.bam",".bam")
                        extend_hisat2 = [">", bam]
                        hisat2.extend(extend_hisat2)
                    utils.write2log(work_dir, " ".join(hisat2))
                    subprocess.call(hisat2)
        else:
            #if running on cluster, get thread count from cluster.yaml
            if slurm == True:
                with open(os.path.join(script_dir,"yaml","slurm.yaml")) as file:
                    slurm_settings = yaml.full_load(file)
                threads = str(slurm_settings["TT-Seq"]["STAR_CPU"])
            
            #create file list of trimming output
            extension = utils.get_extension(work_dir)
            file_list = glob.glob(os.path.join(work_dir, "raw-data","*R1_001." + extension))
            file_list = [x.split(".",1)[0] + "_val_1.fq.gz" for x in file_list]
            file_list = [x.replace("raw-data", "trim") for x in file_list]
            
            #create empty csv file for commands
            Path(os.path.join(work_dir,"slurm",f"slurm_STAR_{genome}.csv")).touch()
            
            for read1 in file_list:
                read1 = read1.replace("_R1_001_val_1.fq.gz","_R2_001_val_2.fq.gz")
                ###to do
      
    #align to selected genome
    align(work_dir,threads, tt_seq_settings,genome, slurm, job_id_trim)
    
    #align to yeast genome (spike-in)
    genome = "R64-1-1"
    align(work_dir,threads, tt_seq_settings,genome, slurm, job_id_trim)


def splitBam(threads, work_dir, genome):
    
    '''
    based on https://www.biostars.org/p/92935/
    
    '''
    puts(colored.green("Generating forward and reverse strand-specific BAM files with samtools"))

    file_list = glob.glob(os.path.join(work_dir, "bam", genome, "*", "*_sorted.bam"))
    
        
    for bam in file_list:
        print(os.path.basename(bam))
        ###forward strand
        fwd1 = bam.replace("_sorted.bam","_fwd1.bam")
        fwd2 = bam.replace("_sorted.bam","_fwd2.bam")
        
        print("\tGenerating forward strand-specific BAM file")
        #alignments of the second in pair if they map to the forward strand
        if not utils.file_exists(fwd1):
            pysam.view("-@",threads,"-b","-f","128","-F","16",bam,"-o",fwd1, catch_stdout=False)
            pysam.index(fwd1)
        
        #alignments of the first in pair if they map to the reverse  strand
        if not utils.file_exists(fwd2):
            pysam.view("-@",threads,"-b","-f","80",bam,"-o",fwd2, catch_stdout=False)
            pysam.index(fwd2)
        
        #merge all forward reads
        fwd = bam.replace("_sorted.bam","_fwd.bam")
        if not utils.file_exists(fwd):
            pysam.merge("-@",threads,fwd,fwd1,fwd2, catch_stdout=False)
            pysam.index(fwd)
        
        ###reverse strand
        rev1 = bam.replace("_sorted.bam","_rev1.bam")
        rev2 = bam.replace("_sorted.bam","_rev2.bam")
        
        print("\tGenerating reverse strand-specific BAM file")
        #alignments of the second in pair if they map to the reverse strand
        if not utils.file_exists(rev1):
            pysam.view("-b","-f","144",bam,"-o", rev1, catch_stdout=False)
            pysam.index(rev1)
        
        #alignments of the first in pair if they map to the forward strand
        if not utils.file_exists(rev2):
            pysam.view("-@",threads,"-b","-f","64","-F","16",bam,"-o",rev2, catch_stdout=False)
            pysam.index(rev2)
        
        #merge all reverse reads
        rev = bam.replace("_sorted.bam","_rev.bam")
        if not utils.file_exists(rev):
            pysam.merge("-@",threads,rev,rev1,rev2, catch_stdout=False)
            pysam.index(rev)
        
        #remove all non-merged bam files
        remove = [fwd1,fwd2,rev1,rev2]
        for file in remove:
            os.remove(file)
            os.remove(file.replace(".bam",".bam.bai"))
            
'''

    puts(colored.green("Generating strand-specific BAM files with samtools"))
    
    #get sorted bam files
    file_list = glob.glob(os.path.join(work_dir, "bam", genome, "*", "*_sorted.bam"))
    
    if len(file_list) == 0:
        return(puts(colored.red("ERROR: No BAM files found to split")))
    
    #load slurm settings  
    with open(os.path.join(script_dir,
                               "yaml",
                               "slurm.yaml")) as file:
            slurm_settings = yaml.full_load(file)
    
    threads = str(slurm_settings["TT-Seq"]["splitBAM_CPU"])
    mem = str(slurm_settings["TT-Seq"]["splitBAM_mem"])
    time = str(slurm_settings["TT-Seq"]["splitBAM_time"])
    account = slurm_settings["groupname"]
    partition = slurm_settings["TT-Seq"]["partition"]
    #email = slurm_settings["email"]
    
    #create csv file with samtools view commands for slurm bash script
    for bam in file_list:
        ##forward strand
        fwd1 = bam.replace("*_sorted.bam","_fwd1.bam")
        fwd2 = bam.replace("*_sorted.bam","_fwd2.bam")
        
        #alignments of the second in pair if they map to the forward strand
        samtools_fwd1 = ["samtools", "view", "-@", threads, "-b", "-f", "128", "-F", "16", bam, "-o", fwd1]
        samtools_fwd1 = " ".join(samtools_fwd1)
        #alignments of the first in pair if they map to the reverse  strand
        samtools_fwd2 = ["samtools", "view", "-@", threads, "-b", "-f", "80", bam, "-o", fwd2]
        samtools_fwd2 = " ".join(samtools_fwd2)
        
        ##reverse strand
        rev1 = bam.replace("*_sorted.bam","_rev1.bam")
        rev2 = bam.replace("*_sorted.bam","_rev2.bam")
        
        #alignments of the second in pair if they map to the reverse strand
        samtools_rev1 = ["samtools", "view", "-b", "-f", "144", bam, "-o", rev1]
        samtools_rev1 = " ".join(samtools_rev1)
        #alignments of the first in pair if they map to the forward strand
        samtools_rev2 = ["samtools", "view", "-@", threads, "-b", "-f", "64", "-F", "16", bam, "-o", rev2]
        samtools_rev2 = " ".join(samtools_rev2)
        
        csv = open(os.path.join(work_dir,"slurm","slurm_splitBAM.csv"), "a")  
        csv.write(samtools_fwd1)
        csv.write(samtools_fwd2)
        csv.write(samtools_rev1)
        csv.write(samtools_rev2)
        csv.close()
    
    #create slurm bash script for splitting bam files
    print("Generating slurm_splitBAM.sh")
    script = os.path.join(work_dir,"slurm","slurm_splitBAM.sh")
    script = open(script, "w")  
    script.write("#!/bin/bash" + "\n")
    script.write("\n")
    script.write("#SBATCH -A " + account + "\n")
    script.write("#SBATCH --mail-type=FAIL" + "\n")
    script.write("#SBATCH --mail-type=END" + "\n")
    script.write("#SBATCH -p " + partition + "\n")
    script.write("#SBATCH -D " + work_dir + "\n")
    script.write("#SBATCH -o slurm/slurm_split_BAM_%a.log" + "\n")
    script.write("#SBATCH -c " + threads + "\n")
    script.write("#SBATCH -t " + time + "\n")
    script.write("#SBATCH --mem=" + mem + "\n")
    script.write("#SBATCH -J " + "split_bam" + "\n")
    script.write("#SBATCH -a " + "1-" + str(len(file_list) * 4) + "\n")
    script.write("\n")
    script.write("sed -n ${SLURM_ARRAY_TASK_ID}p slurm/slurm_splitBAM.csv | bash")
    script.close()
    
    #run slurm script and get job id
    print("Submitting slurm_splitBAM.sh to cluster")
    job_id = subprocess.check_output("sbatch slurm/slurm_splitBAM.sh | cut -d " " -f 4", shell = True)
    job_id = job_id.decode("utf-8")
    
    #load slurm settings for samtools index
    threads = str(slurm_settings["TT-Seq"]["samtools-index_CPU"])
    mem = str(slurm_settings["TT-Seq"]["samtools-index_mem"])
    time = str(slurm_settings["TT-Seq"]["samtools-index_time"])
    
    #create csv file with bam files to be indexed
    extensions = ["*_fwd1.bam", "*_fwd2.bam", "*_rev1.bam", "*_rev2.bam"]
    csv = os.path.join(work_dir,"slurm","slurm_indexBAM.csv")
    os.remove(csv) #delete a preexisting file                   
    
    for extension in extensions:
        file_list = glob.glob(os.path.join(work_dir, "bam", genome, "*", extension))
        for bam in file_list:
            csv = open(csv, "a")
            index = ["samtools","index", "-@", threads]
            csv.write(index)
            csv.close()
    
    #create bash script indexing splitted bam files
    script = os.path.join(work_dir,"slurm","slurm_indexBAM.sh")
    script = open(script, "w")  
    script.write("#!/bin/bash" + "\n")
    script.write("\n")
    script.write("#SBATCH -A " + account + "\n")
    #script.write("#SBATCH ---mail-user=" + email + "\n")
    script.write("#SBATCH --mail-type=FAIL" + "\n")
    script.write("#SBATCH --mail-type=END" + "\n")
    script.write("#SBATCH -p " + partition + "\n")
    script.write("#SBATCH -D " + work_dir + "\n")
    script.write("#SBATCH -o slurm/slurm_index_BAM_%a.log" + "\n")
    script.write("#SBATCH -c " + threads + "\n")
    script.write("#SBATCH -t " + time + "\n")
    script.write("#SBATCH --mem=" + mem + "\n")
    script.write("#SBATCH -J " + "samtools_index" + "\n")
    script.write("#SBATCH -a " + "1-" + str(len(file_list) * 4) + "\n")
    script.write("#SBATCH --dependency=afterok:" + job_id)
    script.write("\n")
    script.write("conda activate ttseq")
    script.write("module load samtools/1.10\n")
    script.write("\n")
    script.write("sed -n %ap slurm/slurm_indexBAM.csv | bash")
    script.close()
    
    #run slurm script and get job id
    job_id = subprocess.check_output("sbatch slurm/slurm_indexBAM.sh | cut -d " " -f 4", shell = True)
    job_id = job_id.decode("utf-8")
    
    #merge all forward and reverse reads
   ''' 
   

def sizeFactors(script_dir, work_dir, slurm=False):
    """
    Creates size factors based on yeast RNA spike-in for generating BigWig files
    Based on https://github.com/crickbabs/DRB_TT-seq/blob/master/bigwig.md with modifications

    """
    puts(colored.green("Generating size factors for normalisation using DESeq2"))
            
    #run DESeq2 to obtain size factors for normalisation
    deseq2 = ["Rscript", os.path.join(script_dir, "R", "tt-seq_sizeFactors.R.R")]
    if slurm == False:
        subprocess.call(deseq2)
    else:
        #load SLURM settings
        with open(os.path.join(script_dir,"yaml","slurm.yaml")) as file:
            slurm_settings = yaml.full_load(file)        

        threads = slurm_settings["RNA-Seq"]["deseq2_CPU"]
        mem = slurm_settings["RNA-Seq"]["deseq2_mem"]
        time = slurm_settings["RNA-Seq"]["deseq2_time"]
        account = slurm_settings["groupname"]
        partition = slurm_settings["RNA-Seq"]["partition"]
        
        #generate SLURM script
        deseq2 = " ".join(deseq2)
        print("Generating SLURM script for caculating scale factors with DESeq2")
        script_ = os.path.join(work_dir,"slurm","slurm_scaleFactors.sh")
        script = open(script_, "w")  
        script.write("#!/bin/bash" + "\n")
        script.write("\n")
        script.write(f"#SBATCH -A {account}\n")
        script.write("#SBATCH --mail-type=BEGIN,FAIL,END" + "\n")
        script.write(f"#SBATCH -p {partition}\n")
        script.write(f"#SBATCH -D {work_dir}\n")
        script.write("#SBATCH -o slurm/slurm_scaleFactors.log" + "\n")
        script.write(f"#SBATCH -c {threads}\n")
        script.write(f"#SBATCH -t {time}\n")
        script.write(f"#SBATCH --mem={mem}\n")
        script.write("#SBATCH -J scaleFactors\n")
        script.write("\n")
        script.write(f"{deseq2}\n")
        script.write("\n")
        script.close()
        
        #send job to cluster
        job_id = subprocess.check_output(f"sbatch {script_} | cut -d ' ' -f 4", shell = True)
        job_id = job_id.decode("UTF-8").replace("\n","")
        print(f"Submitted SLURM script to cluster (job ID {job_id})")
    
    
def ttSeqBigWig(work_dir, threads, tt_seq_settings, genome, slurm):
    """
    Create BigWig files for TT-Seq with scaling factors derived from DESeq2

    """
        
    #load scaling factors
    try:
        scaling_factors = pd.read_csv(os.path.join(work_dir, "scaleFactors.csv"))
    except FileNotFoundError:
        return(puts(colored.red("ERROR: scaleFactors.csv not found")))
    
    #create BigWig directory
    os.makedirs(os.path.join(work_dir,"bigwig",genome), exist_ok = True)
    
    #BigWig function
    def bigWig(work_dir, threads, base, strand, scaling_factor, slurm):
        
        bw_output = f"{base}_{strand}.bigwig"
        bam = f"{base}_{strand}.bam".replace("bigwig", "bam")
        bigwig = ["bamCoverage","--scaleFactor", str(scaling_factor), "-p", threads, "-b", bam, "-o", bw_output]
        puts(colored.green(f"Generating {strand} BigWig file for {os.path.basename(base)}"))
        
        if not utils.file_exists(bw_output):
            if slurm == False:
                utils.write2log(work_dir, " ".join(bigwig))
                subprocess.call(bigwig)
            else:
                #create csv file with bamCoverage commands
                os.makedirs(os.path.join(work_dir,"slurm"), exist_ok = True)
                csv = os.path.join(work_dir,"slurm","slurm_bamCoverage.csv")
                csv = open(csv, "a")  
                csv.write(" ".join(bigwig) +"\n")
                csv.close()  
                
                #load SLURM settings
                with open(os.path.join(script_dir,"yaml","slurm.yaml")) as file:
                    slurm_settings = yaml.full_load(file)
                threads = str(slurm_settings["TT-Seq"]["bamCoverage_CPU"])
                
                mem = str(slurm_settings["TT-Seq"]["bamCoverage_mem"])
                slurm_time = str(slurm_settings["TT-Seq"]["bamCoverage_time"])
                account = slurm_settings["groupname"]
                partition = slurm_settings["TT-Seq"]["partition"]
                
                #create slurm bash script 
                print("Generating slurm_bamCoverage.sh")
                csv = os.path.join(work_dir,"slurm","slurm_bamCoverage.csv")
                commands = int(subprocess.check_output(f"cat {csv} | wc -l", shell = True).decode("utf-8"))
                script_ = os.path.join(work_dir,"slurm","slurm_bamCoverage.sh")
                script = open(script_, "w")  
                script.write("#!/bin/bash" + "\n")
                script.write("\n")
                script.write("#SBATCH -A " + account + "\n")
                script.write("#SBATCH --mail-type=BEGIN,FAIL,END" + "\n")
                script.write("#SBATCH -p " + partition + "\n")
                script.write("#SBATCH -D " + work_dir + "\n")
                script.write("#SBATCH -o slurm/slurm_bamCoverage_%a.log" + "\n")
                script.write("#SBATCH -c " + threads + "\n")
                script.write("#SBATCH -t " + slurm_time + "\n")
                script.write("#SBATCH --mem=" + mem + "\n")
                script.write("#SBATCH -J " + "bamCoverage" +"\n")
                script.write("#SBATCH -a " + "1-" + str(commands) + "\n")
                script.write("\n")
                script.write("sed -n ${SLURM_ARRAY_TASK_ID}p " + csv +" | bash\n")
                script.close()
                script_ = os.path.join(work_dir,"slurm","slurm_bamCoverage.sh")
                
                print("Submitting SLURM script to cluster")
                job_id_bigwig = subprocess.check_output(f"sbatch {script} | cut -d ' ' -f 4", shell = True)
                return(job_id_bigwig)
        
    
    #create scaled BigWig files for both fwd and rev strands
    for index,row in scaling_factors.iterrows():
        os.makedirs(os.path.join(work_dir,"bigwig", genome, row["sample"] ), exist_ok = True)
        base = os.path.join(work_dir,"bigwig", genome, row["sample"] ,row["sample"]) 
        scaling_factor = row["scaleFactors"]
        strands = ["fwd", "rev"]
        for i in strands:
            bigWig(work_dir, threads, base, i, scaling_factor, slurm)
    
    
    #create mean Wig files for all technical replicates with wiggletools
    print("Generating mean Wig files from technical replicates using Wiggletools")
    samples = pd.read_csv(os.path.join(work_dir,"samples.csv"))
    genotypes = set(samples["genotype"])
    conditions = set(samples["condition"])
    
    for genotype in genotypes: 
        for condition in conditions:
            for strand in strands:
                sub_samples = samples[samples["genotype"] == genotype]
                sub_samples = sub_samples[samples["condition"] == condition]
                sub_samples = list(sub_samples["sample"])
                
                in_bigwigs = [os.path.join(work_dir,"bigwig",genome,x,x + "_" + strand +".bigwig") for x in sub_samples]
                
                wig_mean = os.path.join(work_dir,"bigwig", genome, genotype+"_"+condition+f"_{strand}_mean.wig")
                wiggletools = ["wiggletools", "write", wig_mean, "mean", " ".join(in_bigwigs)]
                utils.write2log(work_dir, " ".join(wiggletools))
                subprocess.call(wiggletools)
    
    #generate chrom.sizes file needed for converting Wig to BigWig
    fasta = tt_seq_settings["fasta"][genome]
    print(f"Checking if chrom.sizes file exists for {os.path.basename(fasta)}")
    chrom_sizes = os.path.join(script_dir,"chrom.sizes",os.path.basename(fasta))
    chrom_sizes = chrom_sizes.rsplit(".",1)[0] + ".chrom.sizes"
    os.makedirs(os.path.join(script_dir,"chrom.sizes"), exist_ok = True)
    
    if not utils.file_exists(chrom_sizes):
        fasta_index = fasta + ".fai"
        if not utils.file_exists(fasta_index):
            pysam.faidx(fasta, catch_stdout=False)
            bash = ["cut", "-f","1,2", fasta_index, ">", chrom_sizes]
            utils.write2log(work_dir, " ".join(bash),"Generate chrom.sizes file: ")
            subprocess.call(bash)
    
    #convert Wig to BigWig
    print("Converting Wig to BigWig")
    file_list = glob.glob(os.path.join(work_dir,"bigwig", genome,"*_mean.wig"))
        
    for wig in file_list:
        bw = wig.replace("_mean.wig","_mean.bigwig")
        if not utils.file_exists(bw):
            wig2bigwig = ["wigToBigWig", wig, chrom_sizes, bw]
            utils.write2log(work_dir," ".join(wig2bigwig))
            subprocess.call(wig2bigwig)
            os.remove(wig)
    

def bwQC(work_dir, threads, genome):
    """Quality control for BigWig files using deepTools
    """
    #get BigWig files
    bigwig_all_fwd = glob.glob(os.path.join(work_dir,"bigwig", genome, "*","*fwd.bigwig"))
    bigwig_all_rev = glob.glob(os.path.join(work_dir,"bigwig", genome, "*","*rev.bigwig"))
    bigwig_mean_fwd = glob.glob(os.path.join(work_dir,"bigwig", genome,"*fwd_mean.bigwig"))
    bigwig_mean_rev = glob.glob(os.path.join(work_dir,"bigwig", genome,"*rev_mean.bigwig"))
    
    files = [bigwig_all_fwd, bigwig_all_rev, bigwig_mean_fwd, bigwig_mean_rev]
    names = ["all_fwd","all_rev","mean_fwd","mean_rev"]
    
    #function to create BigWig files
    def multiBigWigSummary(work_dir, threads, file_list, name,genome):
        file_list = os.path.join(work_dir,"bigwig"," ".join(file_list))
        matrix = os.athp.join(work_dir, "bigwig",genome, f"{name}_multiBigWigSummary.npz")
        if not utils.file_exists(matrix):
            command = ["multiBigwigSummary", "bins", "-p", threads, "--smartLabels","-b", file_list, "-o", matrix]
            utils.write2log(work_dir, " ".join(command))
            subprocess.call(command)
    
    def plotPCA(work_dir, threads, genome, matrix):
        plot_file = matrix.replace(".npz","_PCA.pdf")
        if not utils.file_exists(plot_file):
            command = ["plotPCA", "-in", matrix, "--outFileNameData", matrix.replace("npz","tab"), "-o", plot_file]
            utils.write2log(work_dir, " ".join(command))
            subprocess.call(command)
    
    
    #create BigWig Summary for all sample groups
    for list_,name in zip(files,names):
        multiBigWigSummary(work_dir, threads, list_, name, genome)
    
    #create PCA plot for each sample group
    matrix_list = glob.glob(os.path.join(work_dir, "bigwig", genome, "*.npz"))
    for matrix in matrix_list:
        plotPCA(work_dir, threads, genome, matrix)


def metaProfiles(work_dir, threads, tt_seq_settings, genome):
    '''Generate meta plots using DeepTools for TT-Seq data
    '''
    bigwig_mean_fwd = glob.glob(os.path.join(work_dir,"bigwig", genome,"*fwd_mean.bigwig"))
    bigwig_mean_rev = glob.glob(os.path.join(work_dir,"bigwig", genome,"*rev_mean.bigwig")) 
    
    #subset GTF for protein coding genes
    gtf = tt_seq_settings["gtf"][genome]
    gtf_pc = gtf.replace(".gtf", "_pc.gtf")
    
    if not utils.file_exists(gtf_pc):
        grep = ["grep", "protein_coding", gtf, ">", gtf_pc]
        utils.write2log(work_dir, " ".join(grep), "Subset GTF file for only protein coding genes: ")
        subprocess.call(grep)
    
    os.makedirs(os.path.join(work_dir,"deeptools"), exist_ok=True)
    
    def computeMatrix(work_dir, threads, gtf, list_, strand):
        #create compute matrix with deepTools
        sample_names = [os.path.basename(x.replace(f"_{strand}_mean.bigwig","")) for x in list_]
        sample_names = " ".join(sample_names)
        bigwig_mean = " ".join(list_)
        matrix = os.path.join(work_dir, "deeptools",  f"{strand}_plotProfile_matrix.mat.gz")
        if not utils.file_exists(matrix):
            computematrix = ["computeMatrix", "scale-regions","-m", "6000", "-b", "3000","-a", "20000", 
                             "-S", bigwig_mean, "-R", gtf, "--samplesLabel",
                             sample_names, "-p", threads, "-o", matrix]
            utils.write2log(work_dir, " ".join(computematrix))
            subprocess.call(computematrix)
        
    def plotProfile(work_dir, threads, strand):
        #plot meta profiles
        matrix = os.path.join(work_dir, "deeptools",  f"{strand}_plotProfile_matrix.mat.gz")
        out_file = os.path.join(work_dir,"deeptools",f"{strand}_meta_profile_whole_genome.pdf")
        plotProfile = ["plotProfile","-m", matrix, "-out", out_file, "--plotType", "se",
                       "-z", "Protein_coding_genes", "--perGroup"]
        if not utils.file_exists(out_file):
            utils.write2log(work_dir, " ".join(plotProfile))
            subprocess.call(plotProfile)
    
    #generate matrices
    computeMatrix(work_dir, threads, gtf_pc, bigwig_mean_fwd, "fwd")
    computeMatrix(work_dir, threads, gtf_pc, bigwig_mean_rev, "rev")
    
    #generate meta plots
    plotProfile(work_dir, threads, "fwd")
    plotProfile(work_dir, threads, "rev")


def bedBamOverlap(work_dir, genome, slurm, bed):
    """
    Quantify reads from BAM files that overlap with specified BED file

    """
    #merge replicate bam files
    samples = pd.read_csv(os.path.join(work_dir,"samples.csv"))
    genotypes = list(set(samples["genotype"]))
    conditions = list(set(samples["condition"]))
    
    for genotype in genotypes:
        for condition in conditions:
            sub_samples = samples[samples["genotype"] == genotype]
            sub_samples = sub_samples[samples["condition"] == condition]
            sub_samples = list(sub_samples["sample"])
    

def DESeq2(script_dir, genome):
    pass



def txReadThrough(work_dir, threads):
    pass