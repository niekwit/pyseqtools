#!/usr/bin/env python3

import glob
import os
import subprocess
import urllib.request
import sys
from zipfile import ZipFile

script_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(script_dir)
import utils_general as utils


###CHIP-SEQ ANALYSIS SPECIFIC FUNCTIONS


def hisat2(script_dir, work_dir, threads, chip_seq_settings, genome):
    #check for HISAT2
    path = os.environ["PATH"].lower()
    
    if "hisat2" not in path:
        #Check for HISAT2 elsewhere
        hisat2 = [line[0:] for line in subprocess.check_output("find $HOME -name hisat2 ! -path '*/multiqc*'", 
                                                               shell = True).splitlines()]
        
        if hisat2 == []:
            print("WARNING: HISAT2 was not found\nInstalling HISAT2 now")
            url = "https://cloud.biohpc.swmed.edu/index.php/s/hisat2-220-Linux_x86_64/download"
            download_file = os.path.join(script_dir,"hisat2-2.2.0-Linux_x86_64.zip")
            urllib.request.urlretrieve(url,download_file)
            #unzip FastQC file
            with ZipFile(download_file, 'r') as zip_ref:
                zip_ref.extractall(script_dir)
            hisat2 = os.path.join(script_dir, 
                                  "hisat2-2.2.1", 
                                  "hisat2")
            #remove download file
            os.remove(download_file)
    else:
        hisat2 = "hisat2"
    
    #check for HISAT2 index
    index = chip_seq_settings["hisat2"][genome]
    
    if index == "":
        if chip_seq_settings["fasta"][genome] == "":
            try:
                print("WARNING: HISAT2 index was not found")
                print("WARNING: no " + genome + " fasta file was found")
                print("Downloading pre-build index for " + genome)
                
                if genome == "hg19":
                    url = "https://genome-idx.s3.amazonaws.com/hisat/hg19_genome.tar.gz"
                    index_dir = 
                elif genome == "hg38":
                    url = "https://genome-idx.s3.amazonaws.com/hisat/hg38_genome.tar.gz"
                elif genome == "mm9":
                    url = "https://genome-idx.s3.amazonaws.com/hisat/mm10_genome.tar.gz"
                
                download_file = os.path.join(script_dir, os.path.basename(url))
                
                if not utils.file_exists(download_file):
                    urllib.request.urlretrieve(url, 
                                               download_file)
                
                index_dir = os.path.join(script_dir, 
                                         "index", 
                                         "hisat2")
                os.makedirs(index_dir, 
                            exist_ok = True)
                
                #untar index file
                tar_command = "tar -xzf " + download_file + " --directory " + index_dir
                subprocess.run(tar_command,
                               shell = True)
                
                #write index path to yaml
                index_path = glob.glob(index_dir + "*/*/*.ht2")
                index_path = index_path[0].split(".", 1)[0]
                
                with open(os.path.join(script_dir, "yaml", "chip-seq.yaml")) as f:
                doc = yaml.safe_load(f)
                
                doc["hisat2"][genome] = index_dir
                with open(os.path.join(script_dir,"yaml" ,"rna-seq.yaml"), "w") as f:
                    yaml.dump(doc,f)
                
                #remove download file
                os.remove(download_file)
                
            except: #backup method in case index files are offline
                print("WARNING: genome index not available from genome-idx.s3.amazonaws.com")
                print("Downloading and building fasta file needed for building index")
                if genome == "hg19":
                    bash_command = os.path.join(script_dir, 
                                                "bash", 
                                                "build-hg19-fasta.sh")

def bwa():
    pass