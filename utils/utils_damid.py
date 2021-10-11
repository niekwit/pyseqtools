#!/usr/bin/env python3

import glob
import os
import subprocess
import sys
from  builtins import any as b_any

import yaml
try:
    import git #module name is GitPython
except ModuleNotFoundError:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "GitPython"])
    
script_dir = os.path.abspath(os.path.dirname(__file__))
script_dir = os.path.dirname(script_dir)
sys.path.append(os.path.join(script_dir, "utils"))
import utils_general as utils


###DamID-Seq ANALYSIS SPECIFIC FUNCTIONS


def damID(script_dir, work_dir, threads, genome, damid_settings):
    #check for damidseq_pipeline in $PATH
    path = os.environ["PATH"].lower().split(":")

    damid = list(filter(lambda x: "damidseq_pipeline" in x, path))
    
    if len(damid) == 0:
    
        #look for damidseq_pipeline elsewhere
        
        damid = subprocess.check_output("find $HOME -type d -wholename *damidseq_pipeline", 
                                        shell = True).decode("utf-8")
        damid = list(damid.split("\n"))
        damid = [i for i in damid if i] #remove empty list element (\n from find command)
        
        if len(damid) == 0:
            print("WARNING: damidseq_pipeline not found\nInstalling now")
            url = "https://github.com/owenjm/damidseq_pipeline.git"
            git.Git(script_dir).clone(url)
        elif len(damid) > 1:
            sys.exit("ERROR: multiple instances found of damidseq_pipeline")
        elif len(damid) == 1:
            damid, = damid
    elif len(damid) == 1:
        damid, = damid
    
    #get all files to be aligned
    file_list = glob.glob(os.path.join(work_dir, "raw-data","*.gz"))
    
    #check whether there is a dam-only control
    extension = utils.get_extension(work_dir)
    dam_control = b_any(["dam." + extension in x for x in file_list])
    
    if not dam_control:
        sys.exit("ERROR: no dam only control defined (dam.fastq.gz) in raw-data directory")
    
    #get bowtie2 dir
    bowtie2_dir = utils.checkBowtie2(script_dir)
    
    #check for bowtie2 index
    index = damid_settings["bowtie2"][genome]
    if index == "":
        print("Bowtie2 index for " + genome + " not found\nCreating index now")
        fasta = damid_settings["fasta"][genome]
        index_dir = os.path.join(script_dir, "index", "bowtie2", genome)
        os.makedirs(index_dir, exist_ok = True)
        index = os.path.join(index_dir, genome)

        bowtie2_index = [os.path.join(bowtie2_dir, "bowtie2-build"), "--threads",
                        threads, fasta, index]

        #build index
        utils.write2log(work_dir, bowtie2_index, "Bowtie2 build index for " + genome + ": " )
        subprocess.call(bowtie2_index)

        #write index location to cut-run.yaml
        with open(os.path.join(script_dir, "yaml", "damid.yaml")) as f:
            doc = yaml.safe_load(f)

        doc["bowtie2"][genome] = index
        with open(os.path.join(script_dir,"yaml" ,"damid.yaml"), "w") as f:
            yaml.dump(doc,f)

        #reload yaml
        with open(os.path.join(script_dir, "yaml", "damid.yaml")) as file:
            damid_settings = yaml.full_load(file)
        index = damid_settings["bowtie2"][genome]
        
    #check for GATC fragment file
    gatc_gff = damid_settings["GATC_gff"][genome]
    
    if gatc_gff == "":
        print("WARNING: GATC fragment file for " + genome + " not found")
        print("Generating GATC fragment file now")
        ## to be finished
    
    #run pipeline
    command = ["perl", os.path.join(damid, "damidseq_pipeline"), "--gatc_frag_file=" + gatc_gff,
               "--bowtie2_genome_dir=" + index]
    
    utils.write2log(work_dir, command, "DamID pipeline run: ")
    os.chdir(os.path.join(work_dir, "raw-data"))
    subprocess.run(command)    
    
    