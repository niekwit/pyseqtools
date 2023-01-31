#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import glob
import os
import subprocess
import sys
import pkg_resources
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import math
import urllib.request
import shutil
from pathlib import Path
import tempfile
import warnings

import yaml
from clint.textui import colored, puts
import pysam
import gseapy as gp
from gseapy.plot import gseaplot

script_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.join(script_dir, "utils"))
import utils_general as utils
import utils_tt_seq as ttseq


def install_packages(): #check for required python packages; installs if absent
    required = {"pyyaml, cutadapt, multiqc"}
    installed = {pkg.key for pkg in pkg_resources.working_set}
    missing = required - installed
    if missing:
        python = sys.executable
        print("Installing missing required Python3 packages")
        subprocess.check_call([python, '-m', 'pip3', 'install', 
                               *missing], stdout = subprocess.DEVNULL)


def salmonSLURM(work_dir):
    pass


def salmon(salmon_index, threads, work_dir, gtf, fasta, script_dir, settings, reference):
    
    #check for salmon
    path = os.environ["PATH"].lower()
    
    if "salmon" not in path:
        #Check for Salmon elsewhere
        salmon = [line[0:] for line in subprocess.check_output("find $HOME -name salmon ! -path '*/multiqc*'", 
                                                               shell = True).splitlines()]
        
        salmon_file = None
        for i in salmon:
            i = i.decode("utf-8")
            if "bin/salmon" in i:
                salmon_file = i
            
        if salmon_file == None:
            print("WARNING: Salmon was not found\nInstalling Salmon now")
            url = "https://github.com/COMBINE-lab/salmon/releases/download/v1.5.2/salmon-1.5.2_linux_x86_64.tar.gz"
            download_file = os.path.join(script_dir,"salmon-1.5.2_linux_x86_64.tar.gz")
            urllib.request.urlretrieve(url, download_file)
            #untar Salmon download file
            tar_command = "tar -xzf " + download_file + " --directory " + script_dir
            subprocess.run(tar_command,
                           shell = True)
            
            #remove download file
            os.remove(download_file)
            
            salmon_file = os.path.join(script_dir, 
                                       "salmon-1.5.2_linux_x86_64", 
                                       "bin", 
                                       "salmon")
    else:
        salmon_file = "salmon"
    
    
    if salmon_index == "": #Salmon index not found, make it
        print("No Salmon index found: generating Salmon index")
        if os.path.isfile(fasta):
            index_dir = os.path.join(script_dir,"index", "salmon", reference)
            os.makedirs(index_dir, exist_ok = True)
            salmon_index_command = [salmon_file, "index", "-t", fasta, "-i", index_dir, "--gencode"]
            #log commands
            with open(os.path.join(work_dir, "commands.log"), "a") as file:
                file.write("Salmon index: ")
                print(*salmon_index_command, 
                      sep = " ", 
                      file = file)

            subprocess.run(salmon_index_command) #run Salmon index

            #Write salmon index file location to rna-seq.yaml
            with open(os.path.join(script_dir, "yaml", "rna-seq.yaml")) as f:
                doc = yaml.safe_load(f)
            doc["salmon_index"][reference] = index_dir
            with open(os.path.join(script_dir,"yaml" ,"rna-seq.yaml"), "w") as f:
                yaml.dump(doc,f)
        else:
            print("ERROR: no FASTA file specified in rna-seq.yaml")
            sys.exit()

    
    
    
    def salmonSE(work_dir, threads, gtf, salmon_index, reference):
        salmon_output_dir = os.path.join(work_dir,"salmon")
        os.makedirs(salmon_output_dir,exist_ok = True)
        trim_list = glob.glob(os.path.join(work_dir,"trim/*_trimmed.fq.gz"))
        for read1 in trim_list:
            base_read1 = os.path.basename(read1).replace("_trimmed.fq.gz", "") + "-quant"
            salmon_folder_test = os.path.join(salmon_output_dir, base_read1)
            if not utils.file_exists(salmon_folder_test):
                print("Mapping sample " + read1.replace("_trimmed.fq.gz", ""))
                
                out_file = os.path.basename(read1.replace("_trimmed.fq.gz",""))
                salmon_output_file = os.path.join(work_dir, "salmon", out_file)+"-quant"
                salmon_index = settings["salmon_index"][reference] #reload index
                salmon_command = [salmon_file, "quant", "--index", salmon_index, "-l", "A",
                "-g", gtf, "-p", threads, "-r", read1, "--validateMappings",
                "--gcBias", "-o", salmon_output_file]
                with open(os.path.join(work_dir,"commands.log"), "a") as file:
                    file.write("Salmon quant: ")
                    print(*salmon_command, sep = " ", file = file)
                subprocess.run(salmon_command) #Run Salmon quant
    
    def salmonPE(work_dir, threads, gtf, salmon_index, reference):
        salmon_output_dir = os.path.join(work_dir,"salmon")
        os.makedirs(salmon_output_dir, exist_ok = True)
        trim_list = glob.glob(os.path.join(work_dir,"trim","*_val_1.fq.gz"))
        salmon_index = settings["salmon_index"][reference] #reload index
        for read1 in trim_list:
            base_read1 = os.path.basename(read1).replace("_val_1.fq.gz", "") + "-quant"
            salmon_folder_test = os.path.join(salmon_output_dir, base_read1)
            if not utils.file_exists(salmon_folder_test):
                print("Mapping sample " + read1.replace("_val_1.fq.gz", ""))
                read2 = read1.replace("_val_1.fq.gz", "_val_2.fq.gz")
                out_file = os.path.basename(read1.replace("_val_1.fq.gz",""))
                salmon_output_file = os.path.join(work_dir,"salmon",out_file)+"-quant"
                salmon_command = [salmon_file,"quant","--index",salmon_index,"-l","A",
                "-g", gtf,"-p",threads,"-1", read1,"-2",read2,"--validateMappings",
                "--gcBias","-o", salmon_output_file]
                with open(os.path.join(work_dir,"commands.log"), "a") as file:
                    file.write("Salmon quant: ")
                    print(*salmon_command, sep = " ", file = file)
                subprocess.run(salmon_command) #Run Salmon quant
            
    
    if utils.getEND(work_dir) == "PE":
        print("Mapping reads with Salmon (paired-end):")
        salmonPE(work_dir, threads, gtf, salmon_index, reference)
    elif utils.getEND(work_dir) == "SE":
        print("Mapping reads with Salmon (single-end):")
        salmonSE(work_dir, threads, gtf, salmon_index, reference)
    
    #merge Salmon quant files
    salmon_list = glob.glob(os.path.join(work_dir,"salmon","*-quant"))
    merge_output = [os.path.basename(i).replace("-quant","") for i in salmon_list]
    out_file = os.path.join(work_dir,"salmon","salmon_merge.txt")
    
    
    if not utils.file_exists(out_file):
        salmon_merge = [salmon_file, "quantmerge", "--quants", "{" +",".join(salmon_list) + "}", "--names", 
                       "{" + ",".join(merge_output) + "}","-o", out_file]
        utils.write2log(work_dir," ".join(salmon_merge),"Salmon merge quant files: ")
        try:    
            subprocess.run(salmon_merge)
        except:
            print("WARNING: Merging of Salmon quant files failed.")

def plotBar(df,y_label,save_file):
    sns.set_style("white")
    sns.set_style("ticks")
    sns.barplot(x=list(df.keys())[0],
                    y=list(df.keys())[1],
                    data=df,
                    color="cornflowerblue",
                    edgecolor="black",
                    linewidth=1)
    plt.ylabel(y_label)
    plt.xticks(rotation = 'vertical')
    plt.xlabel("")
    plt.ylim(0, 100)
    plt.tight_layout()
    sns.despine()
    plt.savefig(save_file)
    plt.close()


def plotMappingRate(work_dir):
    file_list = glob.glob(os.path.join(work_dir, "salmon", "*", "logs", "salmon_quant.log"))
    save_file = os.path.join(work_dir,"salmon","mapping_rates.pdf")
    mapping_rate = []
    samples = []
    df = pd.DataFrame(columns = ["sample","Mapping rate (%)"],
                      index = np.arange(len(file_list)))

    if not utils.file_exists(save_file): 
        for file in file_list:
                sample = os.path.dirname(file)
                sample = sample.replace(os.path.join(work_dir, "salmon"),"")
                sample = sample.replace("/log", "")
                sample = sample.replace("-quants", "")
                sample = sample.replace("/", "")
                samples.append(sample)
                with open(file,"r") as file:
                    for line in file:
                        if "[info] Mapping rate" in line:
                            rate=line.rsplit(" ", 1)[1]
                            rate=rate.replace("%", "")
                            mapping_rate.append(rate)
        
        df["sample"]=samples
        df["Mapping rate (%)"] = mapping_rate
        df["Mapping rate (%)"] = pd.to_numeric(df["Mapping rate (%)"])
        df = df.sort_values(by = ["sample"],
                          ascending = True,
                          inplace = False).reset_index(drop = True)
        
        plotBar(df, "Mapping rate (%)", save_file)


def plotVolcano(work_dir):
    file_list = glob.glob(os.path.join(work_dir,
                                     "DESeq2",
                                     "*",
                                     "DESeq-output.csv"))
    
    for file in file_list:
        base_name=os.path.basename(os.path.dirname(file))
        
        out_dir=os.path.dirname(file)
        out_file=os.path.join(out_dir,base_name+"-volcano.pdf")
        if not utils.file_exists(out_file):
            print("Generating volcano plot for: "+base_name)
            df=pd.read_csv(file)
            df["log.p.value"]=-np.log10(df["padj"])
            df["label"] = ""
            
            #mark genes are upregulated or downregulated for plotting
            conditions=[(df["log2FoldChange"] > 0.5) & (df["log.p.value"] > 3),
                        (df["label"] != "up") & (df["log2FoldChange"] < -0.5) & (df["log.p.value"] > 3),
                        (df["label"] != "up") & (df["label"] != "down")]
            choices=["Upregulated","Downregulated","None"]
            df["label"]=np.select(conditions, 
                                  choices, 
                                  default=np.nan)
            
            #plot
            sns.set_style("white")
            sns.set_style("ticks")
            sns.scatterplot(data=df, 
                            x="log2FoldChange", 
                            y="log.p.value",
                            alpha=0.5,
                            linewidth=0.25,
                            edgecolor="black",
                            hue="label",
                            palette=["red", "blue", "black"],
                            legend= False)
            plt.xlabel("log2(FC)")
            plt.ylabel("-log(adjusted P value)")
            sns.despine()
            plt.savefig(out_file)
            plt.close()


def plotPCA(work_dir,script_dir):
    out_file = os.path.join(work_dir,"salmon","PCAplot.pdf")
    if not utils.file_exists(out_file):
        PCA_command = ["Rscript", 
                       os.path.join(script_dir, 
                                    "R", 
                                    "rna-seq-pcaplot.R"), 
                       work_dir]
        with open(os.path.join(work_dir,"commands.log"), "a") as file:
            file.write("PCA plot (all samples): ")
            print(*PCA_command, 
                  sep = " ", 
                  file = file)
    
    try:
        subprocess.run(PCA_command)
    except:
        print("PCA plot for all samples failed, check log")
        return(None)


def plotAlignmentRatesSTAR(work_dir,slurm):
    pass


def slurmSTAR(work_dir,script_dir,genome,TE=False):
    '''
    Alignment for RNA-Seq with STAR from trimmed paired-end data
    '''
    if TE == False:
        puts(colored.green(f"STAR alignment for RNA-Seq for {genome}"))
    else:
        puts(colored.green(f"STAR alignment for RNA-Seq for {genome} with relaxed multimapping required for TEtranscripts"))
    
    #get trimmed fastq files
    file_list = glob.glob(os.path.join(work_dir,"trim","*_val_1.fq.gz"))
    if len(file_list) == 0:
        print("ERROR: no trimmed paired-end fastq files found")
        return
    
    #CSV files for commands
    if TE == False:
        csv_star = os.path.join(work_dir,"slurm","star.csv")
        csv_sort = os.path.join(work_dir,"slurm","star_sort.csv")
        csv_index = os.path.join(work_dir,"slurm","star_index.csv")
    else:
        csv_star = os.path.join(work_dir,"slurm","te_star.csv")
        csv_sort = os.path.join(work_dir,"slurm","te_star_sort.csv")
        csv_index = os.path.join(work_dir,"slurm","te_star_index.csv")
    
    csv_list = [csv_star,csv_sort,csv_index]
    
    utils.removeFiles(csv_list) #make sure they do not exist already
    
    #load slurm settings    
    with open(os.path.join(script_dir,"yaml","slurm.yaml")) as file:
        slurm_settings = yaml.full_load(file) 
    
    threads = slurm_settings["RNA-Seq"]["STAR"]["cpu"]
    mem = slurm_settings["RNA-Seq"]["STAR"]["mem"]
    time = slurm_settings["RNA-Seq"]["STAR"]["time"]
    account = slurm_settings["groupname"]
    partition = slurm_settings["RNA-Seq"]["partition"]
    
    slurm = {"threads": threads, 
             "mem": mem,
             "time": time,
             "account": account,
             "partition": partition}
    
    #create commands
    for read1 in file_list:
        read2 = read1.replace("_val_1.fq.gz","_val_2.fq.gz")
        sample = os.path.basename(read1).replace("_val_1.fq.gz","")
        
        #each command should have a unique temp dir otherwise parallel alignments cannot be run
        if TE == False:
            temp_dir = os.path.join(work_dir,f"temp_{sample}")
        else:
            temp_dir = os.path.join(work_dir,f"temp_{sample}_te")
        if os.path.isdir(temp_dir) == True:
            shutil.rmtree(temp_dir)
        
        #create output dir
        if TE == False:
            out_dir = os.path.join(work_dir,"bam",genome,sample)
        else:
            out_dir = os.path.join(work_dir,"bam_te",genome,sample)
        os.makedirs(out_dir, exist_ok = True)
    
        #load RNA-Seq settings   
        with open(os.path.join(script_dir,"yaml","rna-seq.yaml")) as file:
            rna_seq_settings = yaml.full_load(file) 
            
        index = rna_seq_settings["STAR_index"][genome]
            
        #create STAR command
        if TE == False:
            prefix = os.path.join(work_dir,"bam",genome,sample,sample)
        else:
            prefix = os.path.join(work_dir,"bam_te",genome,sample,sample)
        
        star = ["STAR", "--runThreadN", threads,"--runMode", "alignReads", "--genomeDir", index,
                "--readFilesIn", read1, read2, "--readFilesCommand", "zcat", "--quantMode",
                "TranscriptomeSAM", "GeneCounts", "--twopassMode", "Basic", "--outSAMunmapped",
                "None", "--outSAMattrRGline","ID:"+sample,"PU:"+sample,"SM:"+sample,"LB:unknown",
                "PL:illumina", "--outSAMtype","BAM", "Unsorted", "--outTmpDir", temp_dir,
                "--outFileNamePrefix", prefix]
        if TE == True: #special settings for multimapping required for TEtranscripts
            extension = ["--outFilterMultimapNmax","100","--winAnchorMultimapNmax","100"]
            star.extend(extension)
        star = " ".join(star)
        utils.appendCSV(csv_star,star)
        
        #create sort command
        unsorted_bam = os.path.join(out_dir,f"{sample}Aligned.out.bam")
        sorted_bam = unsorted_bam.replace("Aligned.out.bam","_sorted.bam")
        sort_bam = ["samtools","sort","-@",threads,"-o",sorted_bam,unsorted_bam]
        sort_bam = " ".join(sort_bam)
        utils.appendCSV(csv_sort,sort_bam)
        
        #create index command
        index_bam = ["samtools","index","-@",threads,sorted_bam]
        index_bam = " ".join(index_bam)
        utils.appendCSV(csv_index,index_bam)
    
    #generate slurm script
    if TE == False:
        slurm_file = os.path.join(work_dir,"slurm","star.sh")
        utils.slurmTemplateScript(work_dir,"star",slurm_file,slurm,None,True,csv_list)
        
        #submit slurm script to HPC
        job_id_star = utils.runSLURM(work_dir, slurm_file, "star")
    else:
        slurm_file = os.path.join(work_dir,"slurm","te_star.sh")
        utils.slurmTemplateScript(work_dir,"te_star",slurm_file,slurm,None,True,csv_list)
        
        #submit slurm script to HPC
        job_id_star = utils.runSLURM(work_dir, slurm_file, "te_star")
    
    return(job_id_star)
    
        

def STAR(work_dir, threads, script_dir, rna_seq_settings, genome, slurm, job_id_trim=None):
    '''
    Alignment for RNA-Seq with STAR from trimmed paired-end data
    '''
    #genome = genome.split("_")[0]
    
    #function for alignment with STAR
    def align(work_dir, index, threads, genome, slurm):
        #get trimmed fastq files
        file_list = glob.glob(os.path.join(work_dir,"trim","*_val_1.fq.gz"))
        
        #create empty command file for SLURM
        if slurm == True:
            os.makedirs(os.path.join(work_dir,"slurm"), exist_ok = True)
            Path(os.path.join(work_dir,"slurm",f"slurm_STAR_{genome}.csv")).touch()

                     
        for read1 in file_list:
            read2 = read1.replace("_val_1.fq.gz","_val_2.fq.gz")
            
            #create sample name
            sample = os.path.basename(read1).replace("_val_1.fq.gz","")
            
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
                threads = str(slurm_settings["RNA-Seq"]["STAR_CPU"])
            
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
                    csv = open(os.path.join(work_dir,"slurm",f"STAR_{genome}.csv"), "a")  
                    csv.write(" ".join(star) +"\n")
                    csv.close()   
                    
        #index all bam files
        if slurm == False:
            utils.indexBam(work_dir, threads, genome)        
        if slurm == True:
            #load slurm settings  
            mem = str(slurm_settings["RNA-Seq"]["STAR_mem"])
            slurm_time = str(slurm_settings["RNA-Seq"]["STAR_time"])
            account = slurm_settings["groupname"]
            partition = slurm_settings["RNA-Seq"]["partition"]
            
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
                print("Submitting SLURM script to cluster")
                job_id_align = subprocess.check_output(f"sbatch {script} | cut -d ' ' -f 4", shell = True)
                print(f"Submitted SLURM script to cluster (job ID {job_id_align})")
            else:
                print("Submitting slurm script to cluster")
                job_id_align = subprocess.check_output(f"sbatch --dependency=afterok:{job_id_trim} {script} | cut -d ' ' -f 4", shell = True)  
                print(f"Submitted SLURM script to cluster (job ID {job_id_align})")
            
            #sort BAM files
            ttseq.bamSortSLURM(work_dir, job_id_align, genome)
              
    
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
    index = rna_seq_settings["STAR_index"][genome]
    puts(colored.green(f"Aligning fastq files to {genome} with STAR"))
    
    if slurm == True:
        if alignPerformedCluster(work_dir, genome) == False:
            align(work_dir, index, threads, genome, slurm)
        else:
            print(f"Skipping STAR alignment, all output BAM files for {genome} already detected")
    else:
        align(work_dir, index, threads, genome, slurm)


def diff_expr(work_dir,gtf,script_dir,species,pvalue,genome, slurm=False):
    '''
    Differential expression analysis using DESeq2
    '''
    puts(colored.green("Differential expression analysis using DESeq2"))

    
    samples_input=os.path.join(work_dir,"samples.csv")
    if not os.path.exists(samples_input):
        sys.exit("ERROR: "+samples_input+" not found")
    
    
    deseq2_command=["Rscript",os.path.join(script_dir, 
                                           "R", 
                                           "rna-seq_deseq2.R"),
                    work_dir,
                    gtf,
                    script_dir,
                    species,
                    str(pvalue),
                    genome]
    deseq2_command = " ".join(deseq2_command)
    
    
    if slurm == False:
        with open(os.path.join(work_dir,"commands.log"), "a") as file:
            file.write("DESeq2: ")
            print(*deseq2_command, 
                  sep = " ", 
                  file = file)
        print("Running differential expression analysis with DESeq2")
        os.makedirs(os.path.join(work_dir, "DESeq2"),
                    exist_ok = True)
        subprocess.run(deseq2_command)
    else:
        #load SLURM settings
        with open(os.path.join(script_dir,"yaml","slurm.yaml")) as file:
            slurm_settings = yaml.full_load(file)        

        threads = slurm_settings["RNA-Seq"]["deseq2"]["cpu"]
        mem = slurm_settings["RNA-Seq"]["deseq2"]["mem"]
        time = slurm_settings["RNA-Seq"]["deseq2"]["time"]
        account = slurm_settings["groupname"]
        partition = slurm_settings["RNA-Seq"]["partition"]

        #generate SLURM script
        print("Generating SLURM script")
        script_ = os.path.join(work_dir,"slurm","slurm_DESeq2.sh")
        script = open(script_, "w")  
        script.write("#!/bin/bash" + "\n")
        script.write("\n")
        script.write(f"#SBATCH -A {account}\n")
        script.write("#SBATCH --mail-type=BEGIN,FAIL,END" + "\n")
        script.write(f"#SBATCH -p {partition}\n")
        script.write(f"#SBATCH -D {work_dir}\n")
        script.write("#SBATCH -o slurm/slurm_DESeq2.log" + "\n")
        script.write(f"#SBATCH -c {threads}\n")
        script.write(f"#SBATCH -t {time}\n")
        script.write(f"#SBATCH --mem={mem}\n")
        script.write("#SBATCH -J DESeq2\n")
        script.write("\n")
        script.write(f"{deseq2_command}\n")
        script.write("\n")
        script.close()
        
        #send job to cluster
        job_id = subprocess.check_output(f"sbatch {script_} | cut -d ' ' -f 4", shell = True)
        job_id = job_id.decode("UTF-8").replace("\n","")
        print(f"Submitted SLURM script to cluster (job ID {job_id})")


def geneSetEnrichment(work_dir,pvalue,gene_sets):
    file_list = glob.glob(os.path.join(work_dir,
                                     "DESeq2",
                                     "*",
                                     "DESeq-output.csv"))
    
    if len(file_list) == 0:
            print("ERROR: no DESeq2 output files found")
            return(None)
    
    
    def doEnrichr(gene_list, gene_sets, out_dir, output_name):
        out_dir = os.path.join(out_dir,"Enrichr",output_name)
        
        enrichr_results = gp.enrichr(gene_list = gene_list,
                                   gene_sets = gene_sets,
                                   outdir = out_dir)
    
   
    def GSEA(df,gene_set, out_dir):#not working properly yet
        rnk = df.dropna()
        rnk = rnk[["SYMBOL","log2FoldChange"]]
        rnk = rnk.drop_duplicates(subset="SYMBOL")
        pre_res = gp.prerank(rnk = rnk,
                            gene_sets = gene_set,
                            processes = 4,
                            permutation_num = 100,
                            format = "pdf",
                            seed = 6,
                            no_plot = True)
        terms=pre_res.res2d.index
        df_out=pre_res.res2d
        df_out.reset_index(inplace = True)
        os.makedirs(os.path.join(out_dir,"GSEA"), exist_ok = True)
        df_out.to_csv(os.path.join(out_dir,"GSEA",'GSEA.csv'), index = False)
              
        for i in range(10): #plot GSEA plots for top 10 terms
            GSEA_dir=os.path.join(out_dir,"GSEA",gene_set) 
            os.makedirs(GSEA_dir, exist_ok=True)
            try:
                gseaplot(rank_metric=pre_res.ranking,
                         term=terms[i],
                         pheno_pos="Upregulated genes",
                         pheno_neg="Downregulated genes",
                         ofname=os.path.join(GSEA_dir,terms[i]+".pdf"),
                         **pre_res.results[terms[i]])
            except FileNotFoundError:
                continue
    
    logPvalue = -math.log10(pvalue)
    
    for file in file_list:
       df = pd.read_csv(file)
       df["log.p.value"] = -np.log10(df["padj"])
       out_dir = os.path.dirname(file)
              
       df_up = df[(df["log2FoldChange"] > 0.5) & (df["log.p.value"] > logPvalue)]
       upregulated_genes = list(df_up["SYMBOL"])
       
       df_down = df[(df["log2FoldChange"] < -0.5) & (df["log.p.value"] > logPvalue)]
       downregulated_genes = list(df_down["SYMBOL"])
       
       input_list = [upregulated_genes, 
                     downregulated_genes]
       output_names = ["upregulated_genes",
                       "downregulated_genes"]
       
       #run Enrichr
       for x,y in zip(input_list,output_names):
           doEnrichr(x,gene_sets,out_dir,y)
           
       #Run GSEA
       for i in gene_sets:
          GSEA(df,i,out_dir)    
    

def retroElementsSLURM(work_dir,script_dir,genome,dependency):
    '''
    Analysis of transposable elements using TEtranscripts
    '''
    puts(colored.green("TE transcript analysis"))
    
    te_dir = os.path.join(work_dir,"TEtranscripts")
    os.makedirs(te_dir,exist_ok=True)
    
    #load RNA-Seq settings   
    with open(os.path.join(script_dir,"yaml","rna-seq.yaml")) as file:
        rna_seq_settings = yaml.full_load(file) 
    
    gtf = rna_seq_settings["gtf"][genome.split("_",2)[0]]
    te_gtf = rna_seq_settings["TEtranscript"][genome.split("_",2)[0]]
    strand = rna_seq_settings["TEtranscript"]["strand"]
    
    #load slurm settings
    with open(os.path.join(script_dir,"yaml","slurm.yaml")) as file:
        slurm_settings = yaml.full_load(file) 
    
    threads = slurm_settings["RNA-Seq"]["TEtranscript"]["cpu"]
    mem = slurm_settings["RNA-Seq"]["TEtranscript"]["mem"]
    time = slurm_settings["RNA-Seq"]["TEtranscript"]["time"]
    account = slurm_settings["groupname"]
    partition = slurm_settings["RNA-Seq"]["partition"]
    
    slurm = {"threads": threads, 
             "mem": mem,
             "time": time,
             "account": account,
             "partition": partition}
    
    #get sample info
    sample_info = pd.read_csv(os.path.join(work_dir,"samples.csv"))
    all_samples = set(utils.getSampleNames(work_dir))
    
    if len(set(sample_info["condition"])) == 1:
        reference_samples = set(sample_info[sample_info["ref"] == "ref"]["genotype"])
    
    #CSV files for commands
    csv_te = os.path.join(work_dir,"slurm","TEtranscripts.csv")
    csv_list = [csv_te]
    
    utils.removeFiles(csv_list) #make sure they do not exist already
    
    #create TEtranscript commands
    for reference in reference_samples:
        #get all test samples for this reference
        test_samples = all_samples - {reference}
        
        #create command for each individual test sample vs reference
        for test_sample in test_samples:
            if len(set(sample_info["condition"])) == 1:
                test_bams = list(sample_info[sample_info["genotype"] == test_sample]["sample"])
            test_bams = [os.path.join(work_dir,"bam_te",genome,x,f"{x}_sorted.bam") for x in test_bams]
            test_bams = " ".join(test_bams)
    
            if len(set(sample_info["condition"])) == 1:        
                ref_bams = list(sample_info[sample_info["genotype"] == reference]["sample"])
            ref_bams = [os.path.join(work_dir,"bam_te",genome,x,f"{x}_sorted.bam") for x in ref_bams]
            ref_bams = " ".join(ref_bams)
            
            project_name = f"{test_sample}_vs_{reference}"
            out_dir_sample = os.path.join(te_dir,project_name)
            os.makedirs(out_dir_sample,exist_ok=True)
            
            command = f"TEtranscripts -t {test_bams} -c {ref_bams} --GTF {gtf} --TE {te_gtf} --sortByPos --project {project_name} --outdir {out_dir_sample} --stranded {strand}"
            utils.appendCSV(csv_te,command)
            
    
    #generate slurm script
    slurm_file = os.path.join(work_dir,"slurm","te.sh")
    utils.slurmTemplateScript(work_dir,"TEtrx",slurm_file,slurm,None,True,csv_list,dependency)
    
    #submit slurm script to HPC
    job_id_te = utils.runSLURM(work_dir, slurm_file, "TEtrx")
    
    

def retroElements(work_dir,script_dir,rna_seq_settings,genome,slurm=False,threads='1'):
    '''
    Analysis of transposable elements using TEtranscripts
    https://hammelllab.labsites.cshl.edu/software/#TEtranscripts

    '''
    puts(colored.green("TE transcript analysis"))
    
    if slurm == False:
        #sort bam files first
        print("Sorting STAR BAM files")
        file_list = glob.glob(os.path.join(work_dir, "bam", genome, "*", "*Aligned.out.bam"))
        sorted_file_list = [x.replace("Aligned.out.bam","Aligned.out.sorted.bam") for x in file_list]
        
        for bam,sorted_bam in zip(file_list,sorted_file_list):
            if not utils.file_exists(sorted_bam):
                print(f"Sorting {os.path.basename(bam)}")
                pysam.sort("-@", threads, "-o", sorted_bam, bam)
        
        #indexing sorted BAM files
        print("Indexing sorted BAM files")
        sorted_index_file_list = [x + ".bai" for x in sorted_file_list]
        for bam,index in zip(sorted_file_list, sorted_index_file_list):
            print(f"Generating index for {os.path.basename(bam)} ")
            if not utils.file_exists(index):
                pysam.index("-@", threads, bam)
    
    #read samples csv file
    samples = pd.read_csv(os.path.join(work_dir,"samples.csv"))
    gtf = rna_seq_settings["gtf"][genome]
    gtf_te = rna_seq_settings["TE-gtf"][genome]
    
    #get number of comparisons
    comparisons = len(samples.columns) - 1 #number of comparisons +1
    
    #create TEtranscript commands for each comparison
    for i in range(1, comparisons):
        #create base df (no comparison)
        df = samples.iloc[:,0:2]
        
        #add experiment column
        warnings.simplefilter('ignore',lineno=616)
        df["exp"] = samples.iloc[:,1+i:2+i]
        
        #remove irrelevant samples (nan)
        df = df.dropna()
        
        #get control samples
        controls = df[df["exp"] == "control"]
        controls = list(controls["sample"])
        
        #get test samples
        test_samples = df[df["exp"] == "test"]
        test_samples = list(test_samples["sample"])

        #create control bam list
        controls = [os.path.join(work_dir,"bam",genome,x,f"{x}Aligned.out.sorted.bam") for x in controls]
        controls = " ".join(controls)
        
        #create test bam list
        test_samples = [os.path.join(work_dir,"bam",genome,x,f"{x}Aligned.out.sorted.bam") for x in test_samples]
        test_samples = " ".join(test_samples)
        
        #generate output directory
        conditions = list(set(df["condition"]))
        dir_name = " ".join(conditions).replace(" ", "_vs_")
        dir_name = os.path.join(work_dir,"TEtranscripts", dir_name)
        os.makedirs(dir_name, exist_ok = True)
        
        #run TEtranscripts command
        command = ["TEtranscripts", "-c", controls, "-t", test_samples, 
                   "--GTF", gtf, "--TE", gtf_te, "--sortByPos", "--project",
                   os.path.basename(dir_name)]
        if not utils.file_exists(os.path.join(dir_name,os.path.basename(dir_name) + "_gene_TE_analysis.txt")):
            utils.write2log(work_dir, " ".join(command))
            if slurm == False:
                os.chdir(dir_name)
                subprocess.call(command)
                os.chdir(work_dir)
            else:
                print(f"Generating SLURM script for {os.path.basename(dir_name)}")
                #csv = open(os.path.join(work_dir,"slurm",f"slurm_TEtranscripts_{genome}.csv"), "a")  
                #csv.write(" ".join(command) +"\n")
                #csv.close()
                
                #loading SLURM settings
                with open(os.path.join(script_dir,"yaml","slurm.yaml")) as file:
                    slurm_settings = yaml.full_load(file)
                threads = str(slurm_settings["RNA-Seq"]["TEtranscript_CPU"])
                mem = slurm_settings["RNA-Seq"]["TEtranscript_mem"]
                time = slurm_settings["RNA-Seq"]["TEtranscript_time"]
                account = slurm_settings["groupname"]
                partition = slurm_settings["partition"]
                
                #generating SLURM script
                base = os.path.basename(dir_name)
                script_ = os.path.join(work_dir,"slurm",f"slurm_TEtranscript_{genome}_{base}.sh")
                script = open(script_, "w")  
                script.write("#!/bin/bash\n")
                script.write("\n")
                script.write(f"#SBATCH -A {account}\n")
                script.write("#SBATCH --mail-type=BEGIN,FAIL,END" + "\n")
                script.write(f"#SBATCH -p {partition}\n")
                script.write(f"#SBATCH -D {dir_name}\n")
                script.write(f"#SBATCH -o slurm/slurm_TEtranscript_{genome}_{base}.log\n")
                script.write(f"#SBATCH -c {threads}\n")
                script.write(f"#SBATCH -t {time}\n")
                script.write(f"#SBATCH --mem={mem}\n")
                script.write(f"#SBATCH -J TEtranscript_{genome}_{base}\n")
                script.write("\n")
                script.write(" ".join(command) +"\n")
                script.close()
                
                #submit job to SLURM
                print("Submitting SLURM script to cluster")
                script_ = os.path.join(work_dir,"slurm",f"slurm_TEtranscript_{genome}_{base}.sh")
                subprocess.call(f"sbatch {script} | cut -d ' ' -f 4", shell = True)
                    

def BigWig(work_dir, threads, genome, rna_seq_settings, slurm=False):
    '''
    Generate BigWig files from BAM files for RNA-Seq using bamCoverage
    '''    
    puts(colored.green("Generating BigWig files from RNA-Seq data using bamCoverage"))
    
    #create BigWig directory
    os.makedirs(os.path.join(work_dir,"bigwig",genome), exist_ok = True)
    
    file_list = glob.glob(os.path.join(work_dir, "bam", genome, "*", "*_sorted.bam"))
    if len(file_list) == 0:
        puts(colored.red("ERROR: no sorted BAM files found"))
        return()

    
    #according to https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html
    effective_genome_sizes = {"hg19":{"50":"2685511504", "75":"2736124973", "100":"2776919808", "150":"2827437033", "200":"2855464000"},
                              "hg38":{"50":"2701495761", "75":"2747877777", "100":"2805636331", "150":"2862010578", "200":"2887553303"}}
    
    #load bamCoverage settings
    read_length = rna_seq_settings["BigWig"]["read_length"]
    genome_size = effective_genome_sizes[genome][read_length]
    normalizeUsing = rna_seq_settings["BigWig"]["normalizeUsing"]
    binSize = rna_seq_settings["BigWig"]["binSize"]
    
    #load scaling factors if they have been generated
    if os.path.exists(os.path.join(work_dir,"scaleFactors.csv")):
        print("Found scaling factors (will be applied when generating BigWig files)")
        df = pd.read_csv(os.path.join(work_dir,"scaleFactors.csv"))
    
    for bam in file_list:
        if "Aligned.out_sorted.bam" in bam:
            extension = "Aligned.out_sorted.bam"
        else:
            extension = "_sorted.bam"
        
        bigwig = os.path.basename(bam).replace(extension, f"_{normalizeUsing}.bigwig")
        bigwig = os.path.join(work_dir,"bigwig", genome, bigwig)
        
        command = ["bamCoverage", "-p", threads, "--binSize", binSize, "--normalizeUsing",
                   normalizeUsing,"--effectiveGenomeSize", genome_size,"-b", bam]
        out_put = ["-o", bigwig]
        if os.path.exists(os.path.join(work_dir,"scaleFactors.csv")):
            sample = os.path.basename(bam).replace(extension,"")
            scale_factor = df[df["sample"] == sample]["scaleFactors"].to_string().split(" ",1)[1].replace(" ","")
            extend_command = ["--scaleFactor", scale_factor]
            command.extend(extend_command)
            bigwig = bigwig.replace(".bigwig","_scaled.bigwig")
            out_put = ["-o", bigwig]
        
        command.extend(out_put)
        if not utils.file_exists(bigwig):
            print(f"Generating {os.path.basename(bigwig)}")
            utils.write2log(work_dir, " ".join(command))
            subprocess.call(command)
    
    if slurm == False:
        pass
    else:
        #load slurm settings 
        with open(os.path.join(os.path.dirname(script_dir),"yaml","slurm.yaml")) as file:
            slurm_settings = yaml.full_load(file)
        threads = str(slurm_settings["RNA-Seq"]["bamCoverage_CPU"])
        mem = str(slurm_settings["RNA-Seq"]["bamCoverage_mem"])
        slurm_time = str(slurm_settings["RNA-Seq"]["bamCoverage_time"])
        account = slurm_settings["groupname"]
        partition = slurm_settings["partition"]
        
        #create CSV file with bamCoverage commands
        os.makedirs(os.path.join(work_dir,"slurm"), exist_ok = True)
        
        csv = os.path.join(work_dir,"slurm",f"slurm_bamCoverage_{genome}_{normalizeUsing}.csv")
        if os.path.exists(csv):
            os.remove(csv)
          
        #create csv file with bamCoverage commands
        csv = open(csv, "a")  
        for bam in file_list:
            bigwig = os.path.basename(bam).replace("Aligned.out.bam", f"_{normalizeUsing}.bigwig")
            if not utils.file_exists(bigwig):
                command = ["bamCoverage", "-p", threads, "--binSize", binSize, "--normalizeUsing",
                           normalizeUsing,"--effectiveGenomeSize", genome_size,"-b", bam]
                out_put = ["-o", bigwig]
                if os.path.exists(os.path.join(work_dir,"scaleFactors.csv")):
                    sample = os.path.basename(bam).replace("Aligned.out.bam","")
                    scale_factor = df[df["sample"] == sample]["scaleFactors"].to_string().split(" ",1)[1].replace(" ","")
                    extend_command = ["--scaleFactor", scale_factor]
                    command.extend(extend_command)
                    bigwig = bigwig.replace(".bigwig","_scaled.bigwig")
                    out_put = ["-o", bigwig]
                
                command.extend(out_put)
                if not utils.file_exists(bigwig): 
                    csv.write(" ".join(command) +"\n")
        csv.close()
        
        #create slurm bash script 
        print("Generating SLURM script for bamCoverage")
        csv = os.path.join(work_dir,"slurm",f"slurm_bamCoverage_{genome}_{normalizeUsing}.csv")
        commands = subprocess.check_output(f"cat {csv} | wc -l", shell = True).decode("utf-8")
        bigwig_dir = os.path.join(work_dir,"bigwig", genome)
        slurm_log = os.path.join(work_dir, "slurm",'slurm_bamCoverage_%a.log')
        script_ = os.path.join(work_dir,"slurm",f"slurm_bamCoverage_{genome}_{normalizeUsing}.sh")
        script = open(script_, "w")  
        script.write("#!/bin/bash" + "\n")
        script.write("\n")
        script.write(f"#SBATCH -A {account}\n")
        script.write("#SBATCH --mail-type=BEGIN,FAIL,END" + "\n")
        script.write(f"#SBATCH -p {partition}\n")
        script.write(f"#SBATCH -D {bigwig_dir}\n")
        script.write(f"#SBATCH -o {slurm_log}\n")
        script.write(f"#SBATCH -c {threads}\n")
        script.write(f"#SBATCH -t {slurm_time}\n")
        script.write(f"#SBATCH --mem={mem}\n")
        script.write("#SBATCH -J bamCoverage\n")
        script.write(f"#SBATCH -a 1-{commands}\n")
        script.write("\n")
        script.write("sed -n ${SLURM_ARRAY_TASK_ID}p " + f"{csv} | bash\n")
        script.close()
        
        #submit job to cluster
        job_id_bigwig = subprocess.check_output(f"sbatch {script_} | cut -d ' ' -f 4", shell = True)
        job_id_bigwig = job_id_bigwig.decode("UTF-8").replace("\n","")
        print(f"Submitted SLURM script to cluster (job ID {job_id_bigwig})")
        

def rsemIndex(work_dir, script_dir, rna_seq_settings, slurm, rsemIndex):
    '''
    Generate index for RSEM STAR alignment
    '''
    genome = rsemIndex[0]
    read_length = rsemIndex[1]
    gtf = rna_seq_settings["gtf"][genome]
    fasta = rna_seq_settings["FASTA"][genome]
    name = genome + "_" + read_length
    index_name = os.path.join(rsemIndex[2], name, name)
    
    puts(colored.green(f"Generating STAR index via RSEM for {genome} at {index_name}"))
    if slurm == True:
        #load slurm settings
        with open(os.path.join(script_dir,"yaml","slurm.yaml")) as file:
            slurm_settings = yaml.full_load(file)
        threads = str(slurm_settings["RNA-Seq"]["rsem_CPU"])
        mem = str(slurm_settings["RNA-Seq"]["rsem_mem"])
        slurm_time = str(slurm_settings["RNA-Seq"]["rsem_time"])
        account = slurm_settings["groupname"]
        partition = slurm_settings["partition"]
        
        command = ["rsem-prepare-reference", "-p", threads, "--star", "--star-sjdboverhang",
                   str(int(read_length)-1), "--gtf", gtf, fasta, index_name]
        command = " ".join(command)
        
       
        #generate slurm script
        script_index = os.path.join(work_dir, "slurm", "rsem_index.sh")
        slurm_log = os.path.join(work_dir, "slurm", 'rsem_index.log')
        script = open(script_index, "w")  
        script.write("#!/bin/bash\n\n")
        
        script.write(f"#SBATCH -A {account}\n")
        script.write("#SBATCH --mail-type=BEGIN,FAIL,END\n")
        script.write(f"#SBATCH -p {partition}\n")
        script.write(f"#SBATCH -o {slurm_log}\n")
        script.write(f"#SBATCH -c {threads}\n")
        script.write(f"#SBATCH -t {slurm_time}\n")
        script.write(f"#SBATCH --mem={mem}\n")
        script.write("#SBATCH -J rsem_index\n")
              
        script.write(f"{command}\n")
        
        script.close()
        
        #submit script to cluster 
        job_id = subprocess.check_output(f"sbatch {script_index} | cut -d ' ' -f 4", shell = True)
        job_id = job_id.decode("UTF-8").replace("\n","")
        print(f"Submitted SLURM script to cluster (job ID {job_id})")
        
        #add index path to rna-seq.yaml
        print("Writing index name to rna-seq.yaml")
        with open(os.path.join(script_dir, "yaml", "rna-seq.yaml")) as f:
            doc = yaml.safe_load(f)
        doc["RSEM_STAR_index"][name] = index_name
        with open(os.path.join(script_dir, "yaml", "rna-seq.yaml"), "w") as f:
            yaml.dump(doc, f)
        print("Done!")
    
    else:
        pass
      
    
def isoformAnalysis(work_dir, script_dir, rna_seq_settings, genome, slurm, isoformAnalysis):
    '''
    Alternative isoform analysis using RSEM/MISO, based on CRICK scripts, or rMATS

    '''
        
    puts(colored.green("Isoform analysis using MISO or rMATS"))
    
    if utils.pairedEnd():
        print("Paired-end data detected")
    else:
        print("Single-end data detected")
    
    #load sample info
    sample_info = pd.read_csv(os.path.join(work_dir,"samples.csv"))
    
    if isoformAnalysis == "miso":
        if slurm == True:
            print("RSEM/MISO selected")
            #load SLURM settings
            with open(os.path.join(script_dir,"yaml","slurm.yaml")) as file:
                slurm_settings = yaml.full_load(file)
            threads = str(slurm_settings["RNA-Seq"]["rsem_CPU"])
            mem = str(slurm_settings["RNA-Seq"]["rsem_mem"])
            time = str(slurm_settings["RNA-Seq"]["rsem_time"])
            account = slurm_settings["groupname"]
            partition = slurm_settings["partition"]
            strand = slurm_settings["RNA-Seq"]["rsem_strand"]
            
            
            slurm = {"threads": threads, 
                     "mem": mem,
                     "time": time,
                     "account": account,
                     "partition": partition
                     }
            
            #load STAR index (IMPORTANT: created via RSEM!)
            star_index = rna_seq_settings["RSEM_STAR_index"][genome]
            
            #create directory for SLURM commands
            os.makedirs(os.path.join(work_dir,"slurm","RSEM"), exist_ok=True)
            
            #create directory for all output
            rsem_dir = os.path.join(work_dir,"rsem", genome)
                    
            #get unique conditions
            conditions = list(set(sample_info["genotype"]))
            
            #remove pre-existing command files
            csv_merge1 = os.path.join(work_dir, "slurm", "RSEM", "merge1.csv")
            if utils.pairedEnd():
                csv_merge2 = os.path.join(work_dir, "slurm", "RSEM", "merge2.csv")
            csv_rsem = os.path.join(work_dir,"slurm", "RSEM", f"RSEM_{genome}.csv")
            csv_move = os.path.join(work_dir,"slurm", "RSEM", "move.csv")
            csv_sort = os.path.join(work_dir,"slurm", "RSEM", "sort.csv")
            csv_index = os.path.join(work_dir,"slurm", "RSEM", "index.csv")
            csv_picard = os.path.join(work_dir,"slurm", "picard.csv")
            csv_miso_compare = os.path.join(work_dir,"slurm", "miso_compare.csv")
            csv_miso = os.path.join(work_dir,"slurm", "miso.csv")
            
            if utils.pairedEnd():
                csv_list = [csv_merge1,csv_merge2,csv_rsem,csv_move,csv_sort,csv_index,
                        csv_picard,csv_miso_compare,csv_miso]
            else:
                csv_list = [csv_merge1,csv_rsem,csv_move,csv_sort,csv_index,
                        csv_picard,csv_miso_compare,csv_miso]
            
            script_rsem = os.path.join(work_dir, "slurm", f"rsem_{genome}.sh")
            #script_miso_compare = os.path.join(work_dir, "slurm", f"miso_compare_{genome}.sh")
            
            utils.removeFiles(csv_list)
            
                               
            #run RSEM/MISO
            for condition in conditions:
                ###create csv files with all commands (for SLURM array script)###
                #merge replicate fastq files by genotype
                               
                if utils.pairedEnd():
                    read1 = glob.glob(os.path.join(work_dir, "trim", f"{condition}*_val_1.fq.gz"))
                    read1.sort()
                    read2 = glob.glob(os.path.join(work_dir, "trim", f"{condition}*_val_2.fq.gz"))
                    read2.sort()
                                
                    
                    read1_merged = os.path.join(work_dir, "trim", f"{condition}_merged_val_1.fq.gz")
                    command = ["cat", " ".join(read1), ">", read1_merged]
                    utils.appendCSV(csv_merge1, command)
                    
                    read2_merged = os.path.join(work_dir, "trim", f"{condition}_merged_val_2.fq.gz")
                    command = ["cat", " ".join(read2), ">", read2_merged]
                    utils.appendCSV(csv_merge2, command)
                else: #data is single-end
                    read1 = glob.glob(os.path.join(work_dir, "trim", f"{condition}*_trimmed.fq.gz"))
                    read1.sort()
                    read1_merged = os.path.join(work_dir, "trim", f"{condition}_merged_trimmed.fq.gz")
                    
                    command = ["cat", " ".join(read1), ">", read1_merged]
                    utils.appendCSV(csv_merge1, command)
                    
                    
                
                #run RSEM
                if utils.pairedEnd():
                    command = ["rsem-calculate-expression", "--paired-end","--star", "-p", threads,
                            "--strandedness", strand, "--star-output-genome-bam",
                            "--estimate-rspd", "--star-gzipped-read-file",
                            "--time", read1_merged, read2_merged, star_index, condition]
                else:
                    command = ["rsem-calculate-expression","--star", "-p", threads,
                                "--strandedness", strand, "--star-output-genome-bam",
                                "--estimate-rspd", "--star-gzipped-read-file",
                                "--time", read1_merged, star_index, condition]
                    
                    
                utils.appendCSV(csv_rsem, command)
                
                #move rsem files to correct directory
                command = ["mv", os.path.join(work_dir, "*.stat",),
                           os.path.join(work_dir, "*.results",),
                           os.path.join(work_dir, "*.bam",),
                           os.path.join(work_dir, "*.time",),
                           os.path.join(work_dir, f"{condition}*.log",),
                           rsem_dir]
                utils.appendCSV(csv_move, command)
                
                #sort BAM file
                rsem_bam = os.path.join(rsem_dir, f"{condition}.STAR.genome.bam")
                sorted_bam = rsem_bam.replace(".STAR.genome.bam","_sorted.bam")
                command = ["samtools", "sort", "--threads", threads, "-o", sorted_bam, rsem_bam]
                utils.appendCSV(csv_sort, command)
                
                #index sorted BAM file
                command = ["samtools", "index", "-@", threads, sorted_bam]
                utils.appendCSV(csv_index, command)
                
                #determine insert length distribution and SD
                picard_dir = os.path.join(work_dir, "picard")
                os.makedirs(picard_dir, exist_ok=True)
                insert_sizes = os.path.join(picard_dir, f"{condition}.insert_sizes.txt")
                insert_sizes_pdf = os.path.join(picard_dir, f"{condition}.insert_sizes.pdf")
                command = ["picard", "CollectInsertSizeMetrics", f"I={sorted_bam}",
                           f"O={insert_sizes}", f"H={insert_sizes_pdf}", "M=0.5",
                           "MAX_RECORDS_IN_RAM=2000000", "VALIDATION_STRINGENCY=LENIENT"]
                
                utils.appendCSV(csv_picard, command)
                           
                #create csv with sample info to run MISO 
                miso_dir = os.path.join(work_dir,"miso", genome, condition)
                gff_index = rna_seq_settings["MISO_index"][genome.split("_")[0]] #(make sure GFF3 file is indexed first)
                command = [sorted_bam, insert_sizes, miso_dir]
                utils.appendCSV(csv_miso, command)
                
            #create SLURM bash script for RSEM/MISO
            print("Generating SLURM script for RSEM and MISO SUMMARY")
            
            commands = subprocess.check_output(f"cat {csv_rsem} | wc -l", shell = True).decode("utf-8")
            os.makedirs(rsem_dir, exist_ok=True)
            slurm_log = os.path.join(work_dir, "slurm", 'alt_spl_%a.log')
            script = open(script_rsem, "w")  
            script.write("#!/bin/bash\n\n")
            
            script.write(f"#SBATCH -A {account}\n")
            script.write("#SBATCH --mail-type=BEGIN,FAIL,END\n")
            script.write(f"#SBATCH -D {rsem_dir}\n")
            script.write(f"#SBATCH -p {partition}\n")
            script.write(f"#SBATCH -o {slurm_log}\n")
            script.write(f"#SBATCH -c {threads}\n")
            script.write(f"#SBATCH -t {time}\n")
            script.write(f"#SBATCH --mem={mem}\n")
            script.write("#SBATCH -J alt_spl\n")
            script.write(f"#SBATCH -a 1-{commands}\n")
            
            script.write("sed -n ${SLURM_ARRAY_TASK_ID}p " + f"{csv_merge1} | bash\n")
            if utils.pairedEnd():
                script.write("sed -n ${SLURM_ARRAY_TASK_ID}p " + f"{csv_merge2} | bash\n")
            script.write("sed -n ${SLURM_ARRAY_TASK_ID}p " + f"{csv_rsem} | bash\n")
            script.write("sed -n ${SLURM_ARRAY_TASK_ID}p " + f"{csv_move} | bash\n")
            script.write("sed -n ${SLURM_ARRAY_TASK_ID}p " + f"{csv_sort} | bash\n")
            script.write("sed -n ${SLURM_ARRAY_TASK_ID}p " + f"{csv_index} | bash\n")
            script.write("sed -n ${SLURM_ARRAY_TASK_ID}p " + f"{csv_picard} | bash\n")
            
            script.write("source ~/.bashrc\n")
            script.write("conda deactivate\n")
            script.write("conda activate miso\n\n")
            
            script.write("SORTED_BAM=$(sed -n ${SLURM_ARRAY_TASK_ID}p " + f"{csv_miso} | awk" + " '{print $1}')\n")
            script.write("INSERT_SIZE_FILE=$(sed -n ${SLURM_ARRAY_TASK_ID}p " + f"{csv_miso} | awk" + " '{print $2}')\n")
            script.write("INSERT_SIZE=$(sed -n 8p $INSERT_SIZE_FILE | " + "awk '{print $1}')\n")
            script.write("INSERT_SIZE=${INSERT_SIZE%.*}\n") #convert to integer
            script.write("SD=$(sed -n 8p $INSERT_SIZE_FILE | awk '{print $7}')\n")
            script.write("SD=${SD%.*}\n\n") #convert to integer
            script.write("MISO_DIR=$(sed -n ${SLURM_ARRAY_TASK_ID}p " + f"{csv_miso} | awk" + " '{print $3}')\n")
            script.write('MISO_SUMMARY_DIR="${MISO_DIR}/summary"\n\n')
            
            script.write(" ".join(["miso", "--run", gff_index, "$SORTED_BAM","--output-dir $MISO_DIR", "-p", threads, 
                       "--paired-end", "$INSERT_SIZE", "$SD" ,"--read-len", genome.split("_")[1], "\n\n"]))
            
            script.write("summarize_miso --summarize-samples $MISO_DIR $MISO_SUMMARY_DIR\n\n")
            
            script.close()
            
            #submit RSEM/MISO job to cluster 
            job_id_miso = utils.runSLURM(work_dir, script_rsem, "rsem/miso")
              
            #run compare_miso (compare to reference sample)
            ref_condition = sample_info[(sample_info["ref"] == "ref" )]
            ref_condition = list(set(ref_condition["genotype"]))[0]
            
            test_conditions = sample_info[(sample_info["ref"] != "ref" )]
            test_conditions = list(set(test_conditions["genotype"]))
                  
            #generate list for commands
            commands = []
            
            for i in test_conditions:
                samples = [ref_condition, i]
                miso_dirs = " ".join([os.path.join(work_dir,"miso", genome, x) for x in samples])
                miso_out_dir = os.path.join(work_dir, "miso", genome, "_vs_".join(samples))
                names = " ".join(samples)
                                
                command = ["compare_miso", "--compare-samples", miso_dirs, miso_out_dir, 
                           "--comparison-labels", names]
                command = " ".join(command)
                commands.append(command)
            
                
            
            #generate slurm script
            slurm_file = os.path.join(work_dir, "slurm", f"miso_compare_{genome}.sh")
            utils.slurmTemplateScript(work_dir,"miso_comp",slurm_file,slurm,command,False,None,job_id_miso,"miso")
            
            #run slurm script
            job_id_miso_compare = utils.runSLURM(work_dir, slurm_file, "miso-compare")
            '''
            #generate slurm script
            print("Generating SLURM script for MISO COMPARE")
            slurm_log = os.path.join(work_dir, "slurm", 'miso_compare_%a.log')
            commands = subprocess.check_output(f"cat {csv_miso_compare} | wc -l", shell = True).decode("utf-8")
            
            script = open(script_miso_compare, "w")  
            script.write("#!/bin/bash\n\n")
            
            script.write(f"#SBATCH -A {account}\n")
            script.write("#SBATCH --mail-type=BEGIN,FAIL,END\n")
            script.write(f"#SBATCH -p {partition}\n")
            script.write(f"#SBATCH -o {slurm_log}\n")
            script.write(f"#SBATCH -c {threads}\n")
            script.write(f"#SBATCH -t {time}\n")
            script.write(f"#SBATCH --mem={mem}\n")
            script.write("#SBATCH -J miso_compare\n")
            if len(commands) > 1:
                script.write(f"#SBATCH -a 1-{commands}\n")
            
            script.write(f"#SBATCH --dependency=afterok:{job_id}\n\n")
                  
            script.write("source ~/.bashrc\n")
            script.write("conda deactivate\n")
            script.write("conda activate miso\n\n")
            
            if len(commands) > 1:
                script.write("sed -n ${SLURM_ARRAY_TASK_ID}p " + f"{csv_miso_compare} | bash\n")
            else:
                script.write("sed -n 1p " + f"{csv_miso_compare} | bash\n")
               
            script.close()
            
            #submit script to cluster 
            job_id_miso = subprocess.check_output(f"sbatch {script_miso_compare} | cut -d ' ' -f 4", shell = True)
            job_id_miso = job_id_miso.decode("UTF-8").replace("\n","")
            print(f"Submitted SLURM script to cluster (job ID {job_id_miso})")
            '''
            
            
    elif isoformAnalysis == "rmats":
        #check if data is single-end
        singleEnd = os.path.exists(os.path.join(work_dir,".single-end"))
        if singleEnd == True:
            end = "single"
        else:
            end = "paired"
        
        print(f"rMATS selected for {end}-end data")
        if slurm == True:
            
            rmats_dir = os.path.join(work_dir, "rmats", genome)
            os.makedirs(rmats_dir, exist_ok=True)
            
            samples = set(sample_info["genotype"])
            conditions = set(sample_info["condition"])
            reference_sample = list(set(sample_info[sample_info["ref"]=="ref"]["genotype"]))[0]
            
            reference_fastq = glob.glob(os.path.join(work_dir, "trim", f"{reference_sample}*.fq.gz"))
            reference_fastq.sort()
            
            #prepare text file for --s1 rmats flag
            if singleEnd == False:
                #paired-end data
                reference_replicates = int(len(reference_fastq) / 2)  
                fastq_text =[]
                j = [-2,0] #indeces for slicing
                i = 1 #counter
                while i < reference_replicates + 1:
                    j = [x + 2 for x in j]
                    n = ":".join(reference_fastq[j[0]:j[1]])
                    fastq_text.append(n)
                    i += 1
                fastq_text = ",".join(fastq_text)
            else:
                #single-end data
                fastq_text = ",".join(reference_fastq)
                
            s1 = os.path.join(rmats_dir, "s1.txt")
            print(fastq_text, file = open(s1, "w"))
            
            #prepare text files for --s2 rmats flag
            test_samples = list(set(sample_info[sample_info["ref"]!="ref"]["genotype"]))
            s2_list = []
            
            for sample in test_samples:
                test_fastq = glob.glob(os.path.join(work_dir, "trim", f"{sample}*.fq.gz"))
                test_fastq.sort()
                if singleEnd == False:
                    test_replicates = int(len(test_fastq) / 2) #paired-end data
                
                    fastq_text =[]
                    j = [-2,0] #indeces for slicing
                    i = 1 #counter
                    while i < test_replicates + 1:
                        j = [x + 2 for x in j]
                        n = ":".join(test_fastq[j[0]:j[1]])
                        fastq_text.append(n)
                        i += 1
                    fastq_text = ",".join(fastq_text)
                else:
                    fastq_text = ",".join(test_fastq)
                
                s2 = os.path.join(rmats_dir, f"s2_{sample}.txt")
                print(fastq_text, file = open(s2, "w"))
                s2_list.append(s2)
                    
            #prepare rmats commands
            star_index = rna_seq_settings["STAR_index"][genome]
            gtf = rna_seq_settings["gtf"][genome.split("_")[0]]
            out_dir = os.path.join(rmats_dir, f"{reference_sample}_vs_{sample}")
            os.makedirs(out_dir, exist_ok=True)
            tmp_dir = os.path.join(work_dir, "temp")
            os.makedirs(tmp_dir, exist_ok=True)
            read_length = genome.split("_")[1]
                        
            with open(os.path.join(script_dir,"yaml","slurm.yaml")) as file:
                slurm_settings = yaml.full_load(file)
            threads = str(slurm_settings["RNA-Seq"]["rMATS"]["CPU"])
            mem = str(slurm_settings["RNA-Seq"]["rMATS"]["mem"])
            slurm_time = str(slurm_settings["RNA-Seq"]["rMATS"]["time"])
            account = slurm_settings["groupname"]
            partition = slurm_settings["partition"]
            strand = slurm_settings["RNA-Seq"]["rMATS"]["strand"]
            
            extension = ["--gtf", gtf, "--bi", star_index, "--od", out_dir,
                         "-t", end, "--readLength", read_length,
                         "--nthread", threads, "--libType", strand,
                         "--tstat", threads, "--tmp", tmp_dir]
            
            csv_rmats = os.path.join(work_dir, "slurm", f"rmats_{genome}.csv")
            os.makedirs(os.path.join(work_dir, "slurm"), exist_ok=True)
            
            try:
                os.remove(csv_rmats)
            except FileNotFoundError:
                pass
            
            for i in s2_list:
                csv_ = open(csv_rmats, "a")
                
                rmats = ["rmats.py", "--s1", s1, "--s2", s2]
                rmats.extend(extension)
                
                csv_.write(" ".join(rmats) + "\n")
                csv_.close()
            
            #generate slurm script
            print("Generating SLURM script for rMATS")
            slurm_log = os.path.join(work_dir, "slurm", 'rmats_%a.log')
            commands = int(subprocess.check_output(f"cat {csv_rmats} | wc -l", shell = True).decode("utf-8"))
            script_rmats = os.path.join(work_dir, "slurm", f"rmats_{genome}.sh")
            
            script = open(script_rmats, "w")  
            script.write("#!/bin/bash\n\n")
            
            script.write(f"#SBATCH -A {account}\n")
            script.write("#SBATCH --mail-type=BEGIN,FAIL,END\n")
            script.write(f"#SBATCH -D {rmats_dir}\n")
            script.write(f"#SBATCH -p {partition}\n")
            script.write(f"#SBATCH -c {threads}\n")
            script.write(f"#SBATCH -t {slurm_time}\n")
            script.write(f"#SBATCH --mem={mem}\n")
            script.write("#SBATCH -J rMATS\n")
            
            if commands > 1:
                script.write(f"#SBATCH -a 1-{commands}\n")
                script.write("#SBATCH -o " + os.path.join(work_dir,'slurm',f'rmats_{genome}_%a.log') + "\n\n")
                
                script.write("source ~/.bashrc\n")
                script.write("conda deactivate\n")
                script.write("conda activate miso\n\n")
                
                script.write("sed -n ${SLURM_ARRAY_TASK_ID}p " + f"{csv_rmats} | bash\n")
            else:
                script.write("#SBATCH -o " + os.path.join(work_dir,'slurm',f'rmats_{genome}.log') + "\n\n")
                
                script.write("source ~/.bashrc\n")
                script.write("conda deactivate\n")
                script.write("conda activate miso\n\n")
                
                script.write("sed -n 1p " + f"{csv_rmats} | bash\n\n")
                  
            script.close()
            
            #submit script to cluster 
            job_id_rmats = subprocess.check_output(f"sbatch {script_rmats} | cut -d ' ' -f 4", 
                                                   shell = True,
                                                   stderr = subprocess.STDOUT)
            job_id_rmats = job_id_rmats.decode("UTF-8").replace("\n","")
            
            try:
                test = int(job_id_rmats)
                print(f"Submitted SLURM script to cluster (job ID {job_id_rmats})")
                
                #log slurm job id
                utils.SLURM_job_id_log(work_dir, "rMATS", job_id_rmats)
            except ValueError:
                print(job_id_rmats)
            
            #submit rMATS plotting script to HPC
            plot_script = os.path.join(script_dir,"R","rna-seq_plot-rmats.R")
            plot = ["Rscript",plot_script,work_dir]
            plot = [" ".join(plot)]
            
            #load slurm settings
            threads = slurm_settings["RNA-Seq"]["rMATS-plot"]["CPU"]
            mem = slurm_settings["RNA-Seq"]["rMATS-plot"]["mem"]
            time = slurm_settings["RNA-Seq"]["rMATS-plot"]["time"]
            account = slurm_settings["groupname"]
            partition = slurm_settings["RNA-Seq"]["partition"]
            
            slurm = {"threads": threads, 
                     "mem": mem,
                     "time": time,
                     "account": account,
                     "partition": partition
                     }
                
            #generate slurm script
            slurm_file = os.path.join(work_dir, "slurm", f"rMATS-plot_{genome}.sh")
            utils.slurmTemplateScript(work_dir,"rMATS-plot",slurm_file,slurm,plot,False,None,job_id_rmats)
            
            #run slurm script
            job_id_rmatsplot = utils.runSLURM(work_dir, slurm_file, "rMATS-plot")
            
            



