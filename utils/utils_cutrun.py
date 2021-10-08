#!/usr/bin/env python3

import os
import subprocess
import sys
import glob


import yaml

script_dir = os.path.abspath(os.path.dirname(__file__))
script_dir = os.path.dirname(script_dir)
sys.path.append(os.path.join(script_dir, "utils"))
import utils_general as utils


def bowtie(work_dir, script_dir, threads, cutrun_settings, genome):
    
                      
    #check for bowtie index
    index = cutrun_settings["bowtie"][genome]
    if index == "":
        print("Bowtie index for " + genome + " not found\nCreating index now")
        fasta = cutrun_settings["fasta"][genome]
        index_dir = os.path.join(script_dir, "index", "bowtie", genome)
        os.makedirs(index_dir, exist_ok = True)
        index = os.path.join(index_dir, genome)
        
        bowtie_index = [os.path.join(bowtie, "bowtie-build"), "--threads", 
                        threads, fasta, index]
        
        #build index
        utils.write2log(work_dir, bowtie_index, "Bowtie index build for " + genome + ": " )
        subprocess.call(bowtie_index)
        
        #write index location to cut-run.yaml
        with open(os.path.join(script_dir, "yaml", "cut-run.yaml")) as f:
            doc = yaml.safe_load(f)
        
        doc["bowtie"][genome] = index
        with open(os.path.join(script_dir,"yaml" ,"cut-run.yaml"), "w") as f:
            yaml.dump(doc,f)
            
        #reload yaml
        with open(os.path.join(script_dir, "yaml", "cut-run.yaml")) as file:
            cutrun_settings = yaml.full_load(file)
        
    #check if data is paired-end
    paired_end = utils.getEND(work_dir)
    
    #get bowtie index
    bowtie_index = cutrun_settings["bowtie"][genome]
    
    #check for bowtie
    bowtie_dir = utils.checkBowtie(script_dir)
    
    #check for samtools
    samtools_bin = utils.checkSamtools(script_dir)
    
    #check for bedtools
    bedtools_bin = utils.checkBedtools(script_dir)
    
    #load other requirements for alignment
    if paired_end == "SE":
        file_list = glob.glob(os.path.join(work_dir, "trim","*_trimmed.fq.gz"))
        
    elif paired_end == "PE":
        print("Aligning reads with Bowtie (paired-end mode):")
        read1_list = glob.glob(os.path.join(work_dir, "trim","*_R1_001_val_1.fq.gz"))
        read2_list = [x.replace("_R1_001_val_1.fq.gz", "_R2_001_val_2.fq.gz") for x in read1_list]
        output_list = [x.replace("trim", "bam").replace("_R1_001_val_1.fq.gz",".bam") for x in read1_list]
        
        out_dir = os.path.join(work_dir, "bam")
        os.makedirs(out_dir, exist_ok = True)
        
        #load blacklist for selected genome
        blacklist = utils.blackList(script_dir, genome)
                      
        for read1, read2, out_file in zip(read1_list, read2_list, output_list):
            
            if not utils.file_exists(out_file):
                print("Aligning " + os.path.basename(read1).replace("_R1_001_val_1.fq.gz", ""))
                samtools = ["|", samtools_bin, "view", "-q", "15", "-F", "260", "-bS", 
                            "-@", threads, "|", bedtools_bin, "intersect", "-v", "-a",
                            "stdin", "-b", blacklist, "-nonamecheck", "|", samtools_bin,
                            "sort", "-@", threads, ">"]
                
                bowtie_command = [os.path.join(bowtie_dir, "bowtie"), "[","-m", "1", "-v", "2", "-S", 
                          "-I", "0", "-X", "2000", "-p", threads, "--sam-nohead", "]","-x", bowtie_index, "{",
                          "-1", read1, "-2", read2, "}"]
                
                bowtie_command.extend(samtools)
                
                bowtie_command.append(out_file)
                print(" ".join(bowtie_command))
                subprocess.call(bowtie_command)
                
               

def bowtie2(work_dir, script_dir, threads, cutrun_settings, genome):
    pass