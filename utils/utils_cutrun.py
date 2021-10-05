#!/usr/bin/env python3

import os
import subprocess
from subprocess import CalledProcessError
import sys
import glob
import re
import urllib.request

import yaml

script_dir = os.path.abspath(os.path.dirname(__file__))
script_dir = os.path.dirname(script_dir)
sys.path.append(os.path.join(script_dir, "utils"))
import utils_general as utils


def bowtie(work_dir, script_dir, threads, cutrun_settings, genome):
    #check for bowtie1 (exclude bowtie2) in $PATH
    path = os.environ["PATH"].lower().split(":")
    
    bowtie = list(filter(lambda x: re.search(r"bowtie*[-1][^2]", x), path))
    if len(bowtie) == 1:
        bowtie, = bowtie #unpack list
    elif len(bowtie) > 1:
        print("ERROR: multiple instances of Bowtie found:")
        print(bowtie)
        return
    elif len(bowtie) == 0:
        try:
            bowtie = [line[0:] for line in subprocess.check_output('find $HOME -depth -type d -iname *bowtie* | grep "bowtie*-*1[^2]"', shell = True).splitlines()]
            if len(bowtie) > 1:
                print("ERROR: multiple instances of Bowtie found:")
                bowtie = [i.decode("utf-8") for i in bowtie]
                print(bowtie)
                return
            elif len(bowtie) == 1:
                bowtie = bowtie[0].decode("utf-8")
        except CalledProcessError: #when no instance of bowtie is found
            print("WARING: Bowtie not found\nInstalling Bowtie now")     
            if sys.platform in ["linux", "linux2"]:    
                url = "https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.3.1/bowtie-1.3.1-linux-x86_64.zip/download"
                download_file = os.path.join(script_dir, "bowtie-1.3.1-linux-x86_64.zip")
                bowtie = os.path.join(script_dir, "bowtie-1.3.1-linux-x86_64")
            elif sys.platform == "darwin":
                url = "https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.3.1/bowtie-1.3.1-macos-x86_64.zip/download"
                download_file = os.path.join(script_dir, "bowtie-1.3.1-macos-x86_64.zip")
                bowtie = os.path.join(script_dir, "bowtie-1.3.1-macos-x86_64")
            
            #download bowtie zip file
            urllib.request.urlretrieve(url, download_file)
            
            #unzip bowtie
            unzip = ["unzip", "-qq", download_file, "-d", script_dir]
            subprocess.run(unzip)
           
    #check for bowtie index
    index = cutrun_settings["bowtie"][genome]
    if index == "":
        print("Bowtie index for " + genome + " not found\nCreating index now")
        fasta = cutrun_settings["fasta"][genome]
        index_dir = os.path.join(script_dir, "index", "bowtie", genome)
        os.makedirs(index_dir, exist_ok = True)
        index = os.path.join(index_dir, genome)
        
        bowtie_index = [os.path.join(bowtie, "bowtie-build"), "--threads", threads, fasta, index]
        
        #build index
        utils.write2log()
        subprocess.call(work_dir, bowtie_index, "Bowtie index build for " + genome + ": " )
        
        #write index location to cut-run.yaml
        with open(os.path.join(script_dir, "yaml", "cut-run.yaml")) as f:
            doc = yaml.safe_load(f)
        
        doc["hisat2"][genome] = index_dir
        with open(os.path.join(script_dir,"yaml" ,"cut-run.yaml"), "w") as f:
            yaml.dump(doc,f)
            
        #reload yaml
        with open(os.path.join(script_dir, "yaml", "cut-run.yaml")) as file:
            cutrun_settings = yaml.full_load(file)
        
    def align_se():
        pass
    
    def align_pe():
        pass
            
    #check if data is paired-end
    paired_end = utils.getEND(work_dir)
    
    #load other requirements for alignment
    if paired_end == "SE":
        file_list = glob.glob(os.path.join(work_dir, "trim","*_trimmed.fq.gz"))
    elif paired_end == "PE":
        read1_list = glob.glob(os.path.join(work_dir, "trim","*_R1_001_val_1.fq.gz"))