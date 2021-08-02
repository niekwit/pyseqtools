#!/usr/bin/env python3

import warnings
import os
import subprocess
import multiprocessing
import yaml
import sys
import csv
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import pandas as pd
import seaborn as sns
from tqdm.auto import tqdm
import gseapy as gp
import hashlib
import pkg_resources


###GENERAL FUNCTIONS


def checkPythonPackages(): #check for required python packages; installs if absent
    required = {"shyaml", "pyyaml", "pandas", "numpy",
                "matplotlib", "seaborn", "multiqc",
                "cutadapt", "tqdm","gseapy",
                "matplotlib-venn"}
    installed = {pkg.key for pkg in pkg_resources.working_set}
    missing = required - installed
    if missing:
        python = sys.executable
        print("Installing missing required Python3 packages")
        try:
            install_command = [python, '-m', 'pip', 'install', *missing]
            subprocess.check_call(install_command, stdout=subprocess.DEVNULL)
        except:
            sys.exit("ERROR: package installation failed, check log")
    else:
        pass

def checkDeps(script_dir):
    
    #Check for FastQC
    fastqc = [line[0:] for line in subprocess.check_output("find $HOME -name run_fastqc.bat", 
                                                           shell = True).splitlines()]
    try:
        fastqc = fastqc[0].decode("utf-8")
        fastqc = os.path.dirname(fastqc)
        #add fastqc path to hidden yaml file
    except IndexError:
        print("ERROR: FastQC was not found\nInstalling FastQC now")
        url = "https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip"
        download_file = os.path.join(script_dir,"fastqc_v0.11.9.zip")
        urllib.request.urlretrieve(url,download_file)
        #unzip FastQC file
        with ZipFile(download_file, 'r') as zip_ref:
            zip_ref.extractall(script_dir)
        #make fastqc file executable by shell
        fastqc_file = os.path.join(fastqc_dir,"fastqc")
        st = os.stat(fastqc_file)
        os.chmod(fastqc_file, st.st_mode | stat.S_IEXEC)
        #remove download file
        os.remove(download_file)



def checkMd5(work_dir):
    md5sum_file = os.path.join(work_dir,"raw-data", "md5sums.csv")     
   
    def md5(file):
        work_dir = os.getcwd()
        file = os.path.join(work_dir, "raw-data", file)
        hash_md5 = hashlib.md5()
        with open(file, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)
        return(hash_md5.hexdigest())

    if not os.path.exists(os.path.join(work_dir, ".md5summscorrect")):
        if os.path.exists(md5sum_file):
            print("Checking MD5 checksums")
            df = pd.read_csv(md5sum_file)
            df["md5sum_new"] = df["file"].apply(md5)
            
            #compare original checksums with calculated ones
            df["md5sumCorrect"] = df["md5sum"] == df["md5sum_new"]
            
            check_list = df[~df["md5sumCorrect"]]
            if len(check_list) > 0:
                print("Calculated MD5 checksums do not match originals:")
                print(check_list["file"])
                sys.exit(1)
            else:
                print("MD5 checksums correct")
                open(".md5summscorrect", 'a').close()

def write2log(work_dir,command,name):
    with open(os.path.join(work_dir,"commands.log"), "a") as file:
        file.write(name)
        print(*command, sep="",file=file)

def set_threads(args):
    max_threads=str(multiprocessing.cpu_count())
    threads=args["threads"]
    if threads == "max":
        threads=max_threads
    threads=str(threads)
    return threads

def rename(work_dir):
    file=open(os.path.join(work_dir,"rename.config"), "r")
    lines=file.readlines()
    count=0
    for line in lines: #removes newline characters
        lines[count]=line.replace("\n","")
        count+=1

    for line in lines:#rename files
        old_name,new_name=line.split(";")
        os.rename(os.path.join(work_dir,
                    "raw-data",
                    old_name),os.path.join(work_dir,
                                "raw-data",
                                new_name))

def get_extension(work_dir):
    file_list=glob.glob(os.path.join(work_dir,"raw-data","*.gz"))
    test_file=file_list[0]
    extension_index=test_file.index(".",0)
    file_extension=test_file[extension_index:]
    return file_extension

def file_exists(file): #check if file exists/is not size zero
    if os.path.exists(file):
            if os.path.getsize(file) > 0:
                print("Skipping "+file+" (already exists/analysed)")
                return(True)
    else:
        return(False)

def fastqc(work_dir,threads,file_extension,exe_dict):
    fastqc_exe=os.path.join(exe_dict["fastqc"],"fastqc")
    if not os.path.isdir(os.path.join(work_dir,"fastqc")) or len(os.listdir(os.path.join(work_dir,"fastqc"))) == 0:
        os.makedirs(os.path.join(work_dir,"fastqc"),exist_ok=True)
        fastqc_command=fastqc_exe+" --threads "+str(threads)+" --quiet -o fastqc/ raw-data/*"+file_extension
        multiqc_command=["multiqc","-o","fastqc/","fastqc/"]
        #log commands
        with open(os.path.join(work_dir,"commands.log"),"w") as file:
            file.write("FastQC: ")
            print(fastqc_command, file=file)
            file.write("MultiQC: ")
            print(*multiqc_command, sep=" ", file=file)
        try:
            print("Running FastQC on raw data")
            subprocess.run(fastqc_command, shell=True)
        except:
            sys.exit("ERROR: FastQC failed, check logs")
        print("Running MultiQC")
        subprocess.run(multiqc_command)
    else:
        print("Skipping FastQC/MultiQC (already performed)")

