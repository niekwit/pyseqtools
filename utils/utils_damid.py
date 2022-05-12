#!/usr/bin/env python3

import glob
import os
import subprocess
import sys
from  builtins import any as b_any
import urllib.request
import stat
import pandas as pd
import numpy as np


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
    #check for trimmed fastq files
    files = glob.glob(os.path.join(work_dir, "trim","*.gz"))
    #rename files
    [os.rename(x, x.replace("_trimmed.fq.gz",".fq.gz")) for x in files]
    
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
    file_list = glob.glob(os.path.join(work_dir, "trim","*.gz"))
    
    #check whether there is a dam-only control
    dam_control = b_any(["dam" in x for x in file_list])
    
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

        bowtie2_index = ["python3", os.path.join(bowtie2_dir, "bowtie2-build"), "--threads",
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
        #### to be finished ####
    
    #run pipeline
    command = ["perl", os.path.join(damid, "damidseq_pipeline"), "--gatc_frag_file=" + gatc_gff,
               "--bowtie2_genome_dir=" + index]
    
    utils.write2log(work_dir, " ".join(command), "DamID pipeline run: ")
    os.chdir(os.path.join(work_dir, "trim"))
    subprocess.run(command)    

def bedgraph2BigWig(script_dir, work_dir, damid_settings, genome):
    print("Converting BedGraphs (damid_pipeline output) to BigWig format for plotProfile")
    #check for UCSC bedGraphToBigWig in $PATH
    path = os.environ["PATH"].lower().split(":")

    b2b = list(filter(lambda x: "bedgraphtogigwig" in x.lower(), path))
    
    if len(b2b) == 0:
        #check for UCSC bedGraphToBigWig in $HOME
        b2b = subprocess.check_output("find $HOME -name bedGraphToBigWig", 
                                        shell = True).decode("utf-8")
        
        if b2b == "":
            print("WARNING: UCSC bedgraphToBigWig tool not found\nDownloading now")
            
            #download bedGraphToBigWig
            url = "http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig"
            b2b_dir = os.path.join(script_dir, "UCSC-tools")
            os.makedirs(b2b_dir, exist_ok = True)
            b2b = os.path.join(b2b_dir,url.rsplit("/",1)[1])
            urllib.request.urlretrieve(url, b2b)
            
            #make bedGraphToBigWig executable
            st = os.stat(b2b)
            os.chmod(b2b, st.st_mode | stat.S_IEXEC)
        else:
            b2b = b2b.replace("\n","")
            
    elif len(b2b) > 1:
        print("ERROR: multiple instances of bedgraphToBigWig found")
        print("\n".join(b2b))
        sys.exit()
    elif len(b2b) == 1:
        b2b, = b2b
    
    #get list of damidseq_pipeline output files and generate output files
    bedgraphs = glob.glob(os.path.join(work_dir, "trim", "*.bedgraph"))
    bigwigs = [x.replace(".bedgraph",".bw") for x in bedgraphs]
    bigwigs = [x.replace("trim","bigwig") for x in bigwigs]
    
    #get chrom.sizes file for selected genome
    all_chrom_sizes = glob.glob(os.path.join(script_dir, "chrom.sizes", "*chrom.sizes"))
    for x in all_chrom_sizes:
        if genome in x:
            chrom_sizes = x
    
    #convert bedgraph to bigwig
    os.makedirs(os.path.join(work_dir, "bigwig"), exist_ok = True)
    
    for bigwig, bedgraph in zip(bigwigs, bedgraphs):
        if not utils.file_exists(bigwig):
            b2b_com = [b2b, bedgraph, chrom_sizes, bigwig]
            
            utils.write2log(work_dir, b2b_com, "bedgraphToBigWig: ")
            subprocess.call(b2b_com)

def quantileNormalistion(work_dir, script_dir):
    download_dir = os.path.join(script_dir, "DamID_scripts")
    if not os.path.isdir(download_dir):
        url = "https://github.com/AHBrand-Lab/DamID_scripts.git"
        os.makedirs(download_dir)
        git.Git(script_dir).clone(url)
        
def revLogTrans(work_dir):
    '''
    Creates reverse log transformed bedgraph files for visualisation of tracks
    '''
    file_list = glob.glob(os.path.join(work_dir,"trim","*bedgraph.quant.norm.bedgraph"))
    out_files = [x.replace(".bedgraph.quant.norm.bedgraph",".rev.log.quant.norm.bedgraph") for x in file_list]
    
    for outfile,file in zip(out_files,file_list):
        if not utils.file_exists(outfile):
            df = pd.read_csv(file, skiprows=1, header=None, sep="\t")
            #with np.errstate(divide='ignore'): #ignore divide by zero error 
            df[3] = 2**(df[3])
            #df.fillna(0, inplace=True)
        header = subprocess.check_output(["head", "-1", file])
        header = header.decode("utf-8")
        df.to_csv(outfile, sep="\t", header=False, index=False)
        with open(outfile, 'r') as original: data = original.read()
        with open(outfile, 'w') as modified: modified.write(header + data)