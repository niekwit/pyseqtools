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
    def align(work_dir, index, threads, genome, slurm):
        #get list of trimmed fastq files
        if slurm == False:
            file_list = glob.glob(os.path.join(work_dir,"trim","*_val_1.fq.gz"))
        else:
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
            else:
                #create csv files with STAR commands for slurm job
                if not utils.file_exists(bam):
                    print(star)
                    csv = open(os.path.join(work_dir,"slurm",f"slurm_STAR_{genome}.csv"), "a")  
                    csv.write(" ".join(star) +"\n")
                    csv.close()    
        
        '''
        #if alignment has already been done return none
        csv = os.path.join(work_dir,"slurm",f"slurm_STAR_{genome}.csv")
        command_number = subprocess.check_output(f'cat {csv} | wc -l', shell = True).decode("utf-8").replace("\n","")

        if len(command_number) == 0:
            print("Skipping STAR alignment (already performed for all files)")
            return(None)        
        '''
        if slurm == True:
            #load slurm settings  
            mem = str(slurm_settings["TT-Seq"]["STAR_mem"])
            slurm_time = str(slurm_settings["TT-Seq"]["STAR_time"])
            account = slurm_settings["groupname"]
            partition = slurm_settings["TT-Seq"]["partition"]
            #email = slurm_settings["email"]
            
            #create slurm bash script for splitting bam files
                            
            print(f"Generating slurm_STAR_{genome}.sh")
            commands = int(subprocess.check_output(f"cat {csv} | wc -l", shell = True).decode("utf-8"))
            script_ = os.path.join(work_dir,"slurm",f"slurm_STAR_{genome}.sh")
            script = open(script_, "w")  
            script.write("#!/bin/bash" + "\n")
            script.write("\n")
            script.write("#SBATCH -A " + account + "\n")
            script.write("#SBATCH --mail-type=FAIL" + "\n")
            script.write("#SBATCH --mail-type=END" + "\n")
            script.write("#SBATCH -p " + partition + "\n")
            script.write("#SBATCH -D " + work_dir + "\n")
            script.write("#SBATCH -o slurm/slurm_STAR_%a.log" + "\n")
            script.write("#SBATCH -c " + threads + "\n")
            script.write("#SBATCH -t " + slurm_time + "\n")
            script.write("#SBATCH --mem=" + mem + "\n")
            script.write("#SBATCH -J " + "STAR_"+genome + "\n")
            script.write("#SBATCH -a " + "1-" + str(commands) + "\n")
            script.write("\n")
            script.write("sed -n ${SLURM_ARRAY_TASK_ID}p slurm/slurm_STAR.csv | bash")
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
            '''
            CREATE SEPARATE FUNCTION FOR SORTING BAM FILES
            #only sort hg38 bam files (it is better for HTSeq to use unsorted BAM files)
            if slurm == False:
                if genome == "hg38":
                    puts(colored.green("Sorting " + os.path.basename(bam)))
                    
                    if not utils.file_exists(sorted_bam):
                        pysam.sort("--threads", threads,"-o",sorted_bam,bam)
                    
                    #remove unsorted bam file
                    if os.path.exists(sorted_bam):
                        if os.path.exists(bam):
                            os.remove(bam)
            else:
                pass #to do
            '''     
    
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
    if alignPerformedCluster(work_dir, genome) == False:
        print("test")
        align(work_dir, index, threads, genome, slurm)
    else:
        print(f"Skipping STAR alignment, output BAM files for {genome} already detected")
       
    #align trimmed reads to yeast genome (spike-in)
    puts(colored.green("Aligning fastq files to R64-1-1 (spike-in) with STAR"))
    yeast_index = tt_seq_settings["STAR"]["yeast"]
    if alignPerformedCluster(work_dir, "R64-1-1") == False:
        align(work_dir, yeast_index, threads,"R64-1-1", slurm)
    else:
        print("Skipping STAR alignment, output BAM files for R64-1-1 already detected")
    
    #index all bam files
    if slurm == False:
        utils.indexBam(work_dir, threads, genome)
        utils.indexBam(work_dir, threads, "R64-1-1")


def bamSortSTAR():
    pass


def hisat2(work_dir, threads, tt_seq_settings, genome, slurm, job_id_trim):
    """
    Align TT-Seq quality trimmed paired-end files to genome with HISAT2
    """
    if slurm is False:
        read1_list = glob.glob(os.path.join(work_dir,"trim","*_R1_001_val_1.fq.gz"))
        hisat2_index = tt_seq_settings["hisat2"][genome]
        os.path.makedirs(os.path.join(work_dir,"bam", genome), exist_ok=True)
        
        for read1 in read1_list:
            read2 = read1.replace("_R1_001_val_1.fq.gz","_R2_001_val_1.fq.gz")
            bam = os.path.join(work_dir,"bam", genome,os.path.basename(read1.replace("_R1_001_val_1.fq.gz","_sort.bam")))
            sample = os.path.basename(read1.replace("_R1_001_val_1.fq.gz",""))
            puts(colored.green(f"Aligning {sample} to {genome} with HISAT2"))
            
            if not utils.file_exists(bam):
                hisat2 = ["hisat2", "-p", str(threads), "-x", hisat2_index,"-1", read1, "-2", read2,
                          "2>>", os.path.join(work_dir,"align.log"), "|", "samtools", "view", "-q", "15", "-F", 
                          "260", "-b", "-@", str(threads), "-", "|", "samtools", "sort", "-@", str(threads), "-",
                          ">", bam]
                utils.write2log(work_dir, " ".join(hisat2))
                subprocess.call(hisat2)
    else:
        pass #to do


def splitBam(threads, work_dir, genome):
    
    '''
    based on https://www.biostars.org/p/92935/
    
    '''
    puts(colored.green("Creating forward and reverse strand-specific BAM files with samtools"))

    file_list = glob.glob(os.path.join(work_dir, "bam", genome, "*", "*_sorted.bam"))
    
        
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
        
def sizeFactors(work_dir, tt_seq_settings, slurm):
    """
    Based on https://github.com/crickbabs/DRB_TT-seq/blob/master/bigwig.md with modifications

    """
    puts(colored.green("Generating size factors for normalisation using DESeq2"))
    #first prepare htseq-count input for DESeq2
    print("Quantifying reads with HTSeq for DESeq2 input")
    file_list = glob.glob(os.path.join(work_dir,"bam","R64-1-1","*","*Aligned.out.bam"))
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
    file_list = glob.glob(os.path.join(work_dir,"bigwig","*_mean.bigwig"))
    matrix = os.path.join(work_dir, "deeptools",  "matrix.npz")


def trxReadThrough(work_dir, threads):
    pass