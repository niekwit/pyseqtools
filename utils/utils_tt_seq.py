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
            bamSortSLURM(work_dir, genome, job_id_align)
              
    
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
    

def bamSortSLURM(work_dir, genome="hg38", job_id_align=None):
    puts(colored.green("Sorting BAM files using samtools"))
    #create BAM file list (files might not exist yet)
    file_list = glob.glob(os.path.join(work_dir,"bam","*.bam"))
    if len(file_list) == 0:
        file_list = glob.glob(os.path.join(work_dir, "bam", genome, "*", "*Aligned.out.bam"))
        extension = "Aligned.out.bam"
        if len(file_list) == 0:
            extension = utils.get_extension(work_dir)
            file_list = glob.glob(os.path.join(work_dir, "raw-data","*R1_001." + extension))
            file_list = [os.path.basename(x.split(".",1)[0].replace("_R1_001","")) for x in file_list]
            file_list = [os.path.join(work_dir,"bam",genome,x, x  + "Aligned.out.bam") for x in file_list]
    
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
    
    csv = open(os.path.join(work_dir,"slurm","slurm_sortBAM.csv"), "a")  
    for bam in file_list:
        bam_sorted = bam.replace("Aligned.out.bam","Aligned.out_sorted.bam")
        if not utils.file_exists(bam_sorted):
            samtools = ["samtools","sort","-@", threads, bam,"-o", bam_sorted]
            
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
    script_ = os.path.join(work_dir,"slurm","slurm_sortBAM.sh")
    if job_id_align == None:
        job_id = subprocess.check_output(f"sbatch {script_} | cut -d ' ' -f 4", shell = True)
    else:
        job_id = subprocess.check_output(f"sbatch --dependency=afterok:{job_id_align} {script} | cut -d ' ' -f 4", shell = True)
    job_id = job_id.decode("UTF-8").replace("\n","")
    print(f"Script submitted to cluster (job id {job_id})")
        
    
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


def splitBam(work_dir, genome, slurm,threads='1'):
    
    '''
    based on https://www.biostars.org/p/92935/
    
    '''
    puts(colored.green("Generating forward and reverse strand-specific BAM files with samtools"))

    file_list = glob.glob(os.path.join(work_dir, "bam", genome, "*", "*_sorted.bam"))
    
    if slurm == False:    
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
    else:
        
        #load SLURM settings
        with open(os.path.join(script_dir,"yaml","slurm.yaml")) as file:
            slurm_settings = yaml.full_load(file)        

        threads = slurm_settings["TT-Seq"]["samtools"]["cpu"]
        mem = slurm_settings["TT-Seq"]["samtools"]["mem"]
        time = slurm_settings["TT-Seq"]["samtools"]["time"]
        account = slurm_settings["groupname"]
        partition = slurm_settings["TT-Seq"]["partition"]
        
        slurm = {"threads": threads, 
                 "mem": mem,
                 "time": time,
                 "account": account,
                 "partition": partition
                 }
        
        #csv files for commands
        csv_view_fwd1 = os.path.join(work_dir,"slurm","view_fwd1.csv")
        csv_index_fwd1 = os.path.join(work_dir,"slurm","index_fwd1.csv")
        
        csv_view_fwd2 = os.path.join(work_dir,"slurm","view_fwd2.csv")
        csv_index_fwd2 = os.path.join(work_dir,"slurm","index_fwd2.csv")
        
        csv_merge_fwd = os.path.join(work_dir,"slurm","merge_fwd.csv")
        csv_index_fwd = os.path.join(work_dir,"slurm","index_fwd.csv")
        
        csv_view_rev1 = os.path.join(work_dir,"slurm","view_rev1.csv")
        csv_index_rev1 = os.path.join(work_dir,"slurm","index_rev1.csv")
        
        csv_view_rev2 = os.path.join(work_dir,"slurm","view_rev2.csv")
        csv_index_rev2 = os.path.join(work_dir,"slurm","index_rev2.csv")
        
        csv_merge_rev = os.path.join(work_dir,"slurm","merge_rev.csv")
        csv_index_rev = os.path.join(work_dir,"slurm","index_rev.csv")
        
        csv_list = [csv_view_fwd1,csv_index_fwd1,csv_view_fwd2,csv_index_fwd2,
                    csv_merge_fwd,csv_index_fwd,csv_view_rev1,csv_index_rev1,
                    csv_view_rev2,csv_index_rev2,csv_merge_rev,csv_index_rev]
        
        utils.removeFiles(csv_list) #make sure they do not exist already
        
        #create list of bam + index files to be removed post analysis
        remove = []
        
        for bam in file_list:
            ###forward strand
            fwd1 = bam.replace("_sorted.bam","_fwd1.bam")
            fwd2 = bam.replace("_sorted.bam","_fwd2.bam")
            
            #create commands for alignments of the second in pair if they map to the forward strand
            samtools_view_fwd1 = ["samtools","view","-@",threads,"-b","-f","128","-F","16",bam,"-o",fwd1]
            utils.appendCSV(csv_view_fwd1,samtools_view_fwd1)
            
            index_fwd1 = ["samtools","index","-@",threads,fwd1]
            utils.appendCSV(csv_index_fwd1,index_fwd1)
            
            #create commands for alignments of the first in pair if they map to the reverse  strand
            samtools_view_fwd2 = ["samtools","view","-@",threads,"-b","-f","80",bam,"-o",fwd2]
            utils.appendCSV(csv_view_fwd2,samtools_view_fwd2) 
            
            index_fwd2 = ["samtools","index","-@",threads,fwd2]
            utils.appendCSV(csv_index_fwd2,index_fwd2)
            
            #merge all forward reads
            fwd = bam.replace("_sorted.bam","_fwd.bam")
            merge_fwd = ["samtools","merge","-@",threads,fwd,fwd1,fwd2]
            utils.appendCSV(csv_merge_fwd,merge_fwd)
            
            index_fwd = ["samtools","index","-@",threads,fwd]
            utils.appendCSV(csv_index_fwd,index_fwd)
                
            ###reverse strand
            rev1 = bam.replace("_sorted.bam","_rev1.bam")
            rev2 = bam.replace("_sorted.bam","_rev2.bam")
            
            #alignments of the second in pair if they map to the reverse strand
            samtools_view_rev1 = ["samtools","view","-@",threads,"-b","-f","144",bam,"-o",rev1]
            utils.appendCSV(csv_view_rev1,samtools_view_rev1)
            
            index_rev1 = ["samtools","index","-@",threads,rev1]
            utils.appendCSV(csv_index_rev1,index_rev1)
            
            #alignments of the first in pair if they map to the forward strand
            samtools_view_rev2 = ["samtools","view","-@",threads,"-b","-f","64","-F","16",bam,"-o",rev2]
            utils.appendCSV(csv_view_rev2,samtools_view_rev2)   
                        
            index_rev2 = ["samtools","index","-@",threads,rev2]
            utils.appendCSV(csv_index_rev2,index_rev2)
            
            #merge all reverse reads
            rev = bam.replace("_sorted.bam","_rev.bam")
            merge_rev = ["samtools","merge","-@",threads,rev,rev1,rev2]
            utils.appendCSV(csv_merge_rev,merge_rev)
            
            index_rev = ["samtools","index","-@",threads,rev]
            utils.appendCSV(csv_index_rev,index_rev) 
            
            bam_remove = [fwd1,fwd2,rev1,rev2]
            remove.extend(bam_remove)
            remove.extend([x+".bai" for x in bam_remove])
            
        #generate slurm script
        slurm_file = os.path.join(work_dir, "slurm", f"splitBAM_{genome}.sh")
        utils.slurmTemplateScript(work_dir,"splitBAM",slurm_file,slurm,None,True,csv_list,None)
        
        #submit slurm script to cluster
        job_id_split = utils.runSLURM(work_dir, slurm_file, "splitBAM")
        
        #remove all non-merged fw/rev bam files
        remove_command = ["rm"]
        remove_command.extend(remove)
        remove_command = [" ".join(remove_command)]
        
        #generate slurm script
        slurm_file = os.path.join(work_dir, "slurm", f"removeBAM_{genome}.sh")
        utils.slurmTemplateScript(work_dir,"removeBAM",slurm_file,slurm,remove_command,False,None,job_id_split)
        
        #submit slurm script to cluster
        job_id_remove = utils.runSLURM(work_dir, slurm_file, "removeBAM")
        
        #merge replicate fwd/rev strand split BAM files
        csv_merge_replicates = os.path.join(work_dir,"slurm","merge_replicates.csv")
        utils.removeFiles(csv_merge_replicates)
        
        samples = utils.getSampleNames(work_dir)
        strand = ["fwd","rev"]
        
        for i in samples:
            for j in strand:
                bam_files = sorted(glob.glob(os.path.join(work_dir,"bam",genome,f"{i}*",f"*_{j}.bam")))
                out_bam = os.path.join(work_dir,"bam",genome,f"{i}_{j}_merged.bam")
                merge = ["samtools","merge","-f","-@",threads,out_bam]
                merge.extend(bam_files)
                
                #write command to csv
                utils.appendCSV(csv_merge_replicates,merge)
        
        #generate slurm script
        slurm_file = os.path.join(work_dir, "slurm", f"mergeBAM_same-strands_{genome}.sh")
        utils.slurmTemplateScript(work_dir,"merge",slurm_file,slurm,None,True,csv_merge_replicates,job_id_remove)
                                 
        #run slurm script
        job_id_merge = utils.runSLURM(work_dir, slurm_file, "mergeBAM")
        

def scaleFactors(script_dir, work_dir, slurm=False):
    """
    Creates scale factors based on yeast RNA spike-in
    Based on https://github.com/crickbabs/DRB_TT-seq/blob/master/bigwig.md with modifications

    NOTE: Multiply by scale factors and divide by size factors for correction
    """
    puts(colored.green("Generating scale factors for normalisation using DESeq2"))
            
    #run DESeq2 to obtain scale factors for normalisation
    if slurm == False:
        deseq2 = ["Rscript", os.path.join(script_dir, "R", "tt-seq_scaleFactors.R")]
        subprocess.call(deseq2)
    else:
        scale_script = os.path.join(script_dir, "R", "tt-seq_scaleFactors.R")
        deseq2 = f"Rscript {scale_script}"
        
        #load SLURM settings
        with open(os.path.join(script_dir,"yaml","slurm.yaml")) as file:
            slurm_settings = yaml.full_load(file)        

        threads = slurm_settings["RNA-Seq"]["deseq2"]["cpu"]
        mem = slurm_settings["RNA-Seq"]["deseq2"]["mem"]
        time = slurm_settings["RNA-Seq"]["deseq2"]["time"]
        account = slurm_settings["groupname"]
        partition = slurm_settings["RNA-Seq"]["partition"]
        
        slurm = {"threads": threads, 
                 "mem": mem,
                 "time": time,
                 "account": account,
                 "partition": partition
                 }
        
        
        #generate slurm script
        slurm_file = os.path.join(work_dir,"slurm","scaleFactors.sh")
        utils.slurmTemplateScript(work_dir,"scaleFactors",slurm_file,slurm,deseq2)
                                 
        #submit slurm script to HPC
        job_id = utils.runSLURM(work_dir, slurm_file, "scaleFactors")
        
                
        
def scaleBAM(work_dir,slurm,genome):
    ''' Generate scaled BAM files based on scale factors
    
    '''
    if slurm == True:
        #load scale factors
        scale_factors = pd.read_csv(os.path.join(work_dir,"scaleFactors.csv"))
        
        #calculate maximum scale factor
        max_scale_factor = max(scale_factors["scaleFactors"])
        
        #convert each scale factor to 0.0 ≤ FLOAT ≤ 1.0 required for samtools
        scale_factors["subsample_float"] = scale_factors["scaleFactors"] / max_scale_factor
        
        #set seed for subsampling
        seed = "123"
        
        #create samtools float
        scale_factors["subsample_float"] = scale_factors["subsample_float"].astype(str)
        scale_factors["subsample_float"] = seed + "." + scale_factors["subsample_float"].str.rsplit(".",expand=True).iloc[:,1]
        
        #create csv files
        csv_subsample = os.path.join(work_dir,"slurm","subsample.csv")
        csv_subsample_index = os.path.join(work_dir,"slurm","subsample_index.csv")
            
        csv_list = [csv_subsample,csv_subsample_index]
        
        utils.removeFiles(csv_list) #make sure they do not exist already
        
        #load slurm settings
        with open(os.path.join(script_dir,"yaml","slurm.yaml")) as file:
            slurm_settings = yaml.full_load(file)        

        threads = slurm_settings["samtools"]["cpu"]
        mem = slurm_settings["samtools"]["mem"]
        time = slurm_settings["samtools"]["time"]
        account = slurm_settings["groupname"]
        partition = slurm_settings["partition"]
        
        slurm = {"threads": threads, 
                 "mem": mem,
                 "time": time,
                 "account": account,
                 "partition": partition
                 }
        
        #generate samtools commands for scaling and indexing
        for index,row in scale_factors.iterrows():
            sample = row[0]
            subsample_float = row[2]
            bam_in = os.path.join(work_dir,"bam",genome,sample,f"{sample}_sorted.bam")
            bam_subsample_sorted = os.path.join(work_dir,"bam",genome,sample,f"{sample}_subsample_sorted.bam")
            
            #scaling
            if subsample_float == seed + ".0": #sample that does not require scaling will just be copied
                copy = f"cp {bam_in} {bam_subsample_sorted}"
                utils.appendCSV(csv_subsample, copy)
            else:
                samtools = f"samtools view -h -@ {threads} -s {subsample_float} {bam_in} | samtools sort -@ {threads} - > {bam_subsample_sorted}"
                utils.appendCSV(csv_subsample, samtools)
            
            #index
            index = f"samtools index -@ {threads} {bam_subsample_sorted}"
            utils.appendCSV(csv_subsample_index, index)
            
        #generate slurm script
        slurm_file = os.path.join(work_dir, "slurm", f"scaleBAM_{genome}.sh")
        utils.slurmTemplateScript(work_dir,"scale",slurm_file,slurm,None,True,csv_list)                   
        
        #run slurm script
        job_id_scaleBAM = utils.runSLURM(work_dir, slurm_file, "mergeBAM")
        
        return(job_id_scaleBAM)
    
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


def mergeBAM(work_dir, script_dir,genome,extension,slurm=False,dependency=None,threads='1',):
    '''
    Merge replicate BAM files using samtools
    '''
    
    
    if slurm == False:
        
        '''
        puts(colored.green("Merging replicate BAM files using samtools"))
        #load samples.csv
        df = pd.read_csv(os.path.join(work_dir,"samples.csv"))
        df.drop("ref", axis=1, inplace = True)
    
        #merge genotype and condition
        df["merge"] = df["genotype"] + "_" + df["condition"]
        df.drop(["genotype","condition"], axis=1,inplace = True)
        
        #generate sample names for merged BAM files
        sample_names = list(set(df["merge"]))
        
        #merge BAM files
        
        for sample_name in sample_names:
            #get BAM files to be merged
            samples = df[df["merge"] == sample_name]
            samples = samples["sample"]
            out_bam = os.path.join(work_dir, "bam", genome, f"{sample_name}_merged.bam")
            bam_files = [os.path.join(work_dir, "bam", genome, x, f"{x}_sorted.bam") for x in samples]
            print_ = " ".join([os.path.basename(x) for x in bam_files])
            print(f"Merging {print_}")
            
            if not utils.file_exists(out_bam):
                command = ["samtools", "merge", "-@", threads, "-"]
                for i in bam_files:
                    command.append(i)
                command.extend(["|", "samtools", "sort", "-@", threads, "-" ,">", out_bam])
                
                if slurm == False:
                    utils.write2log(work_dir, " ".join(command))
                    subprocess.call(command)
                else:
                    pass
        
        #deduplicate merged BAM files
        file_list = glob.glob(os.path.join(work_dir,"bam", genome,"*_merged.bam"))
        
        print("Removing duplicates from merged BAM files using PICARD")
        for bam in file_list:
            if slurm == False:
                print(os.path.basename(bam))
                picard = utils.checkPicard(script_dir)
                dedup_bam = bam.replace(".bam","_dedup.bam")
                log = bam.replace(".bam","_deduplication.log")
                command = ["java", "-jar", picard, "MarkDuplicates", f"INPUT={bam}" , 
                           f"OUTPUT={dedup_bam}", "REMOVE_DUPLICATES=TRUE", f"METRICS_FILE={log}"]
                utils.write2log(work_dir, " ".join(command))
                if utils.file_exists(dedup_bam):
                    subprocess.call(command)
        '''
    else:
        #get sample names
        sample_names = utils.getSampleNames(work_dir)
        
        #load SLURM settings
        with open(os.path.join(script_dir,"yaml","slurm.yaml")) as file:
            slurm_settings = yaml.full_load(file)        

        threads = slurm_settings["samtools"]["cpu"]
        mem = slurm_settings["samtools"]["mem"]
        time = slurm_settings["samtools"]["time"]
        account = slurm_settings["groupname"]
        partition = slurm_settings["partition"]
        
        slurm = {"threads": threads, 
                 "mem": mem,
                 "time": time,
                 "account": account,
                 "partition": partition
                 }
        
        #csv files for commands
        csv_merge = os.path.join(work_dir,"slurm","merge.csv")
        csv_index = os.path.join(work_dir,"slurm","merge_index.csv")
        
        csv_list = [csv_merge,csv_index]
        
        utils.removeFiles(csv_list)
        
        #create commands for merging strand specific bam files
        for sample in sample_names:
            bams = sorted(glob.glob(os.path.join(work_dir,"bam",genome,f"{sample}*",f"{sample}*{extension}")))
            new_extension = extension.replace(".bam","_merged.bam")            
            
            #create commands to merge and sort replicate BAM files
            out_bam = os.path.join(work_dir,"bam",genome,f"{sample}{new_extension}")
            merge = ["samtools", "merge", "-f","-@", threads, "-"]
            merge.extend(bams)
            merge.extend(["|", "samtools", "sort", "-@", threads, "-" ,">", out_bam])
            utils.appendCSV(csv_merge,merge)
                         
            #create commands to index BAM files
            index_bam = ["samtools", "index", "-@",threads,out_bam]
            utils.appendCSV(csv_index,index_bam)

        #generate slurm script
        slurm_file = os.path.join(work_dir, "slurm", f"merge-bam_{genome}.sh")
        utils.slurmTemplateScript(work_dir,"merge",slurm_file,slurm,None,True,csv_list,dependency)
        
        #run slurm script
        job_id = utils.runSLURM(work_dir, slurm_file, "merge-bam")
        return(job_id)


def readRatio(work_dir, script_dir, genome, slurm=False, threads = "1"):
    """
    Calculate ratios of reads around TSS and TES on replicate merged BAM files.
    Note: use sorted BED file to reduce memory usage

    """
    puts(colored.green("Calculating ratios of reads at TSS vs TES in replicate merged BAM files")) 
    #get merged BAM files
    file_list = glob.glob(os.path.join(work_dir,"bam", genome,"*_merged_dedup.bam"))
  
    if slurm == False:
        #first merge replicate BAM files
        if len(file_list) == 0:
            puts(colored.orange("WARNING: no merged deduplicated BAM files found"))
        mergeBAM(work_dir, genome, threads, slurm)
           
    else:
        if len(file_list) == 0:
            puts(colored.red("ERROR: no merged deduplicated BAM files found"))
            return()
        
        #load SLURM settings
        with open(os.path.join(script_dir,"yaml","slurm.yaml")) as file:
            slurm_settings = yaml.full_load(file)
        threads = slurm_settings["TT-Seq"]["readRatio_CPU"]
        
        mem = slurm_settings["TT-Seq"]["readRatio_mem"]
        slurm_time = slurm_settings["TT-Seq"]["readRatio_time"]
        account = slurm_settings["groupname"]
        partition = slurm_settings["partition"]
        
        #load bed file location
        with open(os.path.join(script_dir,"yaml","tt-seq.yaml")) as file:
            ttseq_settings = yaml.full_load(file)
        bed = ttseq_settings["readRatio"]["bed"]
        genome_file = ttseq_settings["genome_file"][genome]
        
        #create csv file with commands
        os.makedirs(os.path.join(work_dir,"slurm"), exist_ok = True)
        
        csv = os.path.join(work_dir,"slurm",f"slurm_bedtools-intersect_{genome}.csv")
        if os.path.exists(csv):
            os.remove(csv)
        
        csv = open(csv, "a")  
        for bam in file_list:
            out_bed = os.path.join(work_dir, "readRatio", os.path.basename(bam).replace("_merged_dedup.bam",".bed"))
            if not utils.file_exists(out_bed):
                #base bedtools command
                command = ["bedtools", "intersect", "-sorted", "-g", genome_file, "-s", "-a", bed, "-b" ]
                #create final command
                command.extend([bam, ">", out_bed])
                csv.write(" ".join(command) +"\n")
            else:
                continue
        csv.close()
        
        #create output directory
        os.makedirs(os.path.join(work_dir, "readRatio"), exist_ok = True)
        
        #create SLURM script
        print("Generating SLURM script for bedtools intersect")
        csv = os.path.join(work_dir,"slurm",f"slurm_bedtools-intersect_{genome}.csv")
        commands = str(subprocess.check_output(f"cat {csv} | wc -l", shell = True).decode("utf-8"))
        
        if commands != "0":
            script_ = os.path.join(work_dir,"slurm",f"slurm_bedtools-intersect_{genome}.sh")
            script = open(script_, "w")  
            script.write("#!/bin/bash" + "\n")
            script.write("\n")
            script.write(f"#SBATCH -A {account}\n")
            script.write("#SBATCH --mail-type=BEGIN,FAIL,END\n")
            script.write(f"#SBATCH -p {partition}\n")
            script.write(f"#SBATCH -D {work_dir}\n")
            script.write("#SBATCH -o slurm/slurm_bedtools-intersect_%a.log" + "\n")
            script.write(f"#SBATCH -c {threads}\n")
            script.write(f"#SBATCH -t {slurm_time}\n")
            script.write(f"#SBATCH --mem={mem}\n")
            script.write("#SBATCH -J bedtools-intersect\n")
            script.write(f"#SBATCH -a 1-{commands}\n")
            script.write("\n")
            script.write("sed -n ${SLURM_ARRAY_TASK_ID}p " + csv +" | bash\n")
            script.close()
            
            print("Submitting SLURM script to cluster")
            job_id_intersect = subprocess.check_output(f"sbatch {script_} | cut -d ' ' -f 4", shell = True)
            job_id_intersect = job_id_intersect.decode("utf-8").replace("\n","")
            print(f"Submitted SLURM script to cluster (job ID {job_id_intersect})")
        else:
            print("Bedtools intersect has already been run for all samples")
            job_id_intersect = None
        
        #get read count for each TSS/TES
        bed_list = [os.path.join(work_dir, "readRatio", os.path.basename(x).replace("_merged_dedup.bam",".bed")) for x in file_list]
        
        csv = os.path.join(work_dir,"slurm",f"slurm_countreads-bedtools_{genome}.csv")
        if os.path.exists(csv):
            os.remove(csv)
        
        csv = open(csv, "a")  
        for bed in bed_list:
            out_bed_tss = bed.replace(".bed","_count_TSS.txt")
            out_bed_tes = bed.replace(".bed","_count_TES.txt")
            
            if not utils.file_exists(out_bed_tss):
                command = [os.path.join(script_dir, "bash", "readRatio.sh"), bed, out_bed_tss, out_bed_tes]
                csv.write(" ".join(command) +"\n")
            else:
                continue
        csv.close()
        
        #create SLURM script
        print("Generating SLURM script for counting reads at TES/TSS")
        csv = os.path.join(work_dir,"slurm",f"slurm_countreads-bedtools_{genome}.csv")
        commands = str(subprocess.check_output(f"cat {csv} | wc -l", shell = True).decode("utf-8"))
        script_ = os.path.join(work_dir,"slurm",f"slurm_countreads-bedtools_{genome}.sh")
        script = open(script_, "w")  
        script.write("#!/bin/bash" + "\n")
        script.write("\n")
        script.write(f"#SBATCH -A {account}\n")
        script.write("#SBATCH --mail-type=BEGIN,FAIL,END\n")
        script.write(f"#SBATCH -p {partition}\n")
        script.write(f"#SBATCH -D {work_dir}\n")
        script.write(f"#SBATCH -o slurm/slurm_countreads-bedtools_{genome}_%a.log" + "\n")
        script.write(f"#SBATCH -c {threads}\n")
        script.write(f"#SBATCH -t {slurm_time}\n")
        script.write(f"#SBATCH --mem={mem}\n")
        script.write("#SBATCH -J countreads-bedtools\n")
        script.write(f"#SBATCH -a 1-{commands}\n")
        script.write("\n")
        script.write("sed -n ${SLURM_ARRAY_TASK_ID}p " + csv +" | bash\n")
        script.close()
        
        #submit to cluster
        script = os.path.join(work_dir,"slurm",f"slurm_countreads-bedtools_{genome}.sh")
        if job_id_intersect is None:
            print("Submitting SLURM script to cluster")
            job_id_count = subprocess.check_output(f"sbatch {script} | cut -d ' ' -f 4", shell = True)
        else:
            print("Submitting SLURM script to cluster")
            job_id_count = subprocess.check_output(f"sbatch --dependency=afterok:{job_id_intersect} {script} | cut -d ' ' -f 4", shell = True)
            print(f"Submitted SLURM script to cluster (job ID {job_id_count})")
            
        #plot read count with R
        


def DESeq2(work_dir, script_dir, genome, slurm=False):
    '''
    Differential transcript analysis for TT-Seq using DESeq2
    '''
    puts(colored.green("Differential expression analysis using DESeq2"))
    
    samples = os.path.join(work_dir,"samples.csv")
    if not os.path.exists(samples):
        puts(colored.red(f"ERROR: {samples} not found"))
        return()
    
    if slurm == False:
        pass
    else:
        pass


def txReadThrough(work_dir, threads):
    pass





def ngsPlotSlurm(work_dir,genome):
    '''
    Creation of metagene profiles using ngs.plot (https://github.com/crickbabs/DRB_TT-seq/blob/master/metaprofiles.md)
    '''
    puts(colored.green("Generating metagene profiles using ngs.plot"))
    
    #create ngs.plot output dir
    os.makedirs(os.path.join(work_dir,"ngsplot"),exist_ok=True)
    
    #load slurm settings
    with open(os.path.join(script_dir,"yaml","slurm.yaml")) as file:
        slurm_settings = yaml.full_load(file)        

    threads = slurm_settings["samtools"]["cpu"]
    mem = slurm_settings["samtools"]["mem"]
    time = slurm_settings["samtools"]["time"]
    account = slurm_settings["groupname"]
    partition = slurm_settings["partition"]
    
    slurm = {"threads": threads, 
             "mem": mem,
             "time": time,
             "account": account,
             "partition": partition
             }
    
    #### create scaled bam files
    job_id_scaleBAM = scaleBAM(work_dir,slurm,genome)
    
    #### merge replicate scaled bam files
    job_id_merge = mergeBAM(work_dir, script_dir,genome,"_subsample_sorted.bam",True,threads,job_id_scaleBAM)
    
    #### create mate1 bam files
       
    #create csv file names for commands
    csv_mate1 = os.path.join(work_dir,"slurm","mate1.csv")
    csv_reheader = os.path.join(work_dir,"slurm","reheader.csv")
    csv_reheader_index = os.path.join(work_dir,"slurm","index_reheader.csv")
    
    csv_ngsplot = os.path.join(work_dir,"slurm","ngsplot.csv")
    
    csv_list = [csv_mate1,csv_reheader,csv_reheader_index, csv_ngsplot]
    
    utils.removeFiles(csv_list) #make sure they do not exist already
    
    #create first mate only BAM files and reheader for hg38 (required for ngs.plot)
    if genome == "hg38":
        
        bam_files = glob.glob(os.path.join(work_dir,"bam",genome,"*[!fwdrev]_merged.bam"))
        
        for bam in bam_files:
           bam_mate1 = bam.replace(".bam","_mate1.bam")
           bam_mate1_reheader = bam.replace(".bam","_reheader_mate1.bam")
            
           mate1 = f"samtools view -@ {threads} -h -b -f 64 {bam} -o {bam_mate1}"
           utils.appendCSV(csv_mate1,mate1)
           
           reheader = f"samtools view -@ $THREADS -H ${bam_mate1} | sed -e 's/SN:\([0-9XY]*\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' | samtools reheader - {bam_mate1} > {bam_mate1_reheader}"
           utils.appendCSV(csv_reheader,reheader)
           
           index = f"samtools index -@ {threads} {bam_mate1_reheader}"
           utils.appendCSV(csv_reheader_index,index)
        
        #generate slurm script
        slurm_file = os.path.join(work_dir, "slurm", f"mate1_{genome}.sh")
        utils.slurmTemplateScript(work_dir,"mate1",slurm_file,slurm,None,True,csv_list,job_id_merge)
        
        #submit slurm script to HPC
        job_id_mate1 = utils.runSLURM(work_dir, slurm_file, "mate1")
        
        #get sample names
        samples = utils.getSampleNames(work_dir)
        
        #### run ngsplot
        #create ngs.plot config file
        #make sure config file does not exist pre-plotting
        config = os.path.join(work_dir,"ngsplot","config.txt")
        utils.removeFiles(config)
            
        #add sample info to config file
        for sample in samples:
            bam = os.path.join(work_dir,"bam",genome,f"{sample}_merged.bam")
            line = f'{bam}\t"-1"\t"{sample}"'
            utils.appendCSV(config, line)
        
        #load slurm settings
        threads = slurm_settings["TT-Seq"]["ngsplot"]["cpu"]
        mem = slurm_settings["TT-Seq"]["ngsplot"]["mem"]
        time = slurm_settings["TT-Seq"]["ngsplot"]["time"]
        account = slurm_settings["groupname"]
        partition = slurm_settings["partition"]
        
        slurm = {"threads": threads, 
                 "mem": mem,
                 "time": time,
                 "account": account,
                 "partition": partition
                 }
        
        #prepare ngsplot commands
        regions = ["genebody","tss","tes"]
        strands = ["both","same","opposite"]
        
        for region in regions:
            for strand in strands:
                output = os.path.join(work_dir,"ngsplot",f"{region}_{strand}")
                ngsplot = f"ngs.plot.r -G {genome} -R {region} -C {config} -O {output} -P {threads} -SS {strand} -SE 1 -L 5000 -F chipseq -D ensembl"
                utils.appendCSV(csv_ngsplot, ngsplot)
                
        #generate slurm script
        slurm_file = os.path.join(work_dir, "slurm", f"ngsplot_{genome}.sh")
        utils.slurmTemplateScript(work_dir,"ngsplot",slurm_file,slurm,None,True,csv_ngsplot,job_id_mate1)
        
        #submit slurm script to HPC
        job_id_ngsplot = utils.runSLURM(work_dir, slurm_file, "ngsplot1")
        
        
        
        

        
        
        
        