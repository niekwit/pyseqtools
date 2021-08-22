#!/usr/bin/env python3

import os
import subprocess
from subprocess import DEVNULL
import multiprocessing
import sys
import glob
import hashlib
import urllib.request
from zipfile import ZipFile
import stat
from  builtins import any as b_any
import tarfile

import yaml
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from tqdm.auto import tqdm
import pysam


###GENERAL FUNCTIONS


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
    max_threads = str(multiprocessing.cpu_count())
    threads = args["threads"]
    if threads == "max":
        threads=max_threads
    threads = str(threads)
    return threads


def rename(work_dir):
    file = open(os.path.join(work_dir,"rename.config"), "r")
    lines = file.readlines()
    count = 0
    for line in lines: #removes newline characters
        lines[count] = line.replace("\n","")
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
    
    if len(file_list) == 0:
        sys.exit("ERROR: no fastq files found")
    
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


def checkSamtools(script_dir):
    #Check for samtools in $PATH
    path = os.environ["PATH"].lower()
    
    if "samtools" not in path:
        #Check for samtools elsewhere
        samtools = [line[0:] for line in subprocess.check_output("find $HOME -name samtools", 
                                                               shell = True).splitlines()]
        
        samtools_file = None
        for i in samtools:
            i = i.decode("utf-8")
            if "bin/samtools" in i:
                samtools_file = i
                return(samtools_file)
            
        if samtools_file == None:
            print("WARNING: samtools was not found\nInstalling samtools now")
            
            url = "https://github.com/samtools/samtools/releases/download/1.13/samtools-1.13.tar.bz2"
            download_file = os.path.join(script_dir,
                                         os.path.basename(url))
            
            urllib.request.urlretrieve(url,
                                       download_file)
            
            #untar download file
            tar = tarfile.open(download_file, "r:bz2")  
            tar.extractall(script_dir)
            tar.close()
            
            os.remove(download_file)
            
            #make samtools
            samtools_dir = os.path.join(script_dir,
                                     "samtools")
            os.makedirs(samtools_dir,
                        exist_ok = True)
            
            os.chdir(os.path.join(script_dir,
                                  os.path.basename(url).replace(".tar.bz2", "")))
            
            subprocess.run(["./configure", "--prefix=" + samtools_dir])
            subprocess.run(["make"])
            subprocess.run(["make", "install"])
            
            samtools_file = os.path.join(script_dir,
                                         "samtools",
                                         "bin",
                                         "samtools")
            
            #remove downloaded files
            samtools_dir = download_file.rsplit(".", 2)[0]
            
            subprocess.run(["rm", "-rf", samtools_dir])
            
            return(samtools_file)

def checkBedtools(script_dir):
    path = os.environ["PATH"].lower()
    
    if "bedtools" not in path:
        #Check for bedtools elsewhere
        bedtools = [line[0:] for line in subprocess.check_output("find $HOME -name bedtools", 
                                                               shell = True).splitlines()]
        
        try:
            bedtools = bedtools[0].decode("utf-8")
            return(bedtools)
        except:
            print("WARNING: Bedtools was not found\nInstalling Bedtools now")
           
            url = "https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary"
            
            os.makedirs(os.path.join(script_dir,
                                     "bedtools-2.30"),
                        exist_ok = True)
            
            download_file = os.path.join(script_dir,
                                         "bedtools-2.30",
                                         os.path.basename(url))
            
            
            urllib.request.urlretrieve(url,
                                       download_file)
            
            bedtools = download_file.rsplit(".", 2)[0]
            os.rename(download_file,
                      bedtools)
            
            #make bedtools executable
            st = os.stat(bedtools)
            os.chmod(bedtools, st.st_mode | stat.S_IEXEC)
            
            return(bedtools)
    

def checkBwa(script_dir):
    pass


def checkPicard(script_dir):
    #check for Java Runtime Environment
    try:
        subprocess.run(["java", "--version",], stdout=DEVNULL)
    except:
        sys.exit("ERROR: No Java Runtime Environment available\nUbuntu installation: https://ubuntu.com/tutorials/install-jre#1-overview")
    
    #check for Picard
    path = os.environ["PATH"].lower()
    
    if "picard" not in path:
        #Check for Picard elsewhere
        picard = [line[0:] for line in subprocess.check_output("find $HOME -name picard.jar ! -path '*/multiqc*'", 
                                                               shell = True).splitlines()]
        try:
            picard = picard[0].decode("utf-8")
            return(picard)
        except:
            print("WARNING: Picard was not found\nInstalling Picard now")
            url = "https://github.com/broadinstitute/picard/releases/download/2.26.0/picard.jar"
            os.makedirs(os.path.join(script_dir,
                                     "picard"),
                        exist_ok = True)
            download_file = os.path.join(script_dir,
                                         "picard",
                                         os.path.basename(url))
            urllib.request.urlretrieve(url, download_file)
            picard = download_file
            return(picard)
            

def deduplicationBam(script_dir, work_dir, threads):
    #check if Picard is installed
    picard = checkPicard(script_dir)
    
    file_list = glob.glob(os.path.join(work_dir,
                                       "bam",
                                       "*-sort-bl.bam"))
    print("Performing deduplication of BAM files")
    for bam in file_list:
        dedup_output = bam.replace("-sort-bl.bam",
                                   "-sort-bl-dedupl.bam")
        if not file_exists(dedup_output):
            
            log_name = os.path.basename(bam).split("-", 1)[0] + "-dedup_metrics.log"
            log_name = os.path.join(work_dir,
                                    "bam",
                                    log_name)
            picard_command = "java -jar " + picard + " MarkDuplicates INPUT=" + bam + " OUTPUT=" + dedup_output + " REMOVE_DUPLICATES=TRUE METRICS_FILE=picard_metrics.log"
            
            write2log(work_dir, picard_command, "Deduplication: ")
                
            subprocess.run(picard_command,
                           shell = True)
        
    #plot number of reads after deduplication
    #get counts from non-deduplicated bam files
    file_list = glob.glob(os.path.join(work_dir,
                                       "bam",
                                       "*-sort-bl.bam"))
    
    column_names = [os.path.basename(bam).replace("-sort-bl.bam","") for bam in file_list]
    df = pd.DataFrame(columns = column_names)
    
    for bam in file_list:
          count = pysam.view("-@", str(threads) ,"-c", "-F" "260", bam)
          column = os.path.basename(bam).replace("-sort-bl.bam","")
          df.loc[1, column] = count
    
    df["condition"] = "pre-deduplication"
    
    #get counts for deduplicated bam files
    file_list = glob.glob(os.path.join(work_dir,
                                       "bam",
                                       "*-sort-bl-dedupl.bam"))
           
    for bam in file_list:
          count = pysam.view("-@", str(threads) ,"-c", "-F" "260", bam)
          column = os.path.basename(bam).replace("-sort-bl-dedupl.bam","")
          df.loc[2, column] = count
    
    df.loc[2, "condition"] = "deduplicated"
    
    #create df for plotting
    df_melt = pd.melt(df, id_vars = ["condition"], 
                      value_vars = column_names)
    df_melt["value"] = pd.to_numeric(df_melt["value"])
    df_melt["value"] = df_melt["value"] / 1000000
    
    #create plot
    save_file = os.path.join(work_dir,
                             "bam",
                             "read-counts-deduplication.pdf")
    
    sns.catplot(x = 'variable', y = 'value', 
               hue = 'condition', 
               data = df_melt, 
               kind = 'bar', 
               legend_out = False,
               edgecolor = "black",)
    plt.ylabel("Uniquely mapped read count (millions)")
    plt.xticks(rotation = 45, ha="right")
    plt.xlabel("")
    plt.ylim((0, df_melt["value"].max() * 1.3))
    plt.legend(title = None,
               frameon = False)
    plt.tight_layout()
    plt.savefig(save_file)
    

def createBigWig(work_dir, threads):
    #creates BigWig files for all existing BAM files
    print("Generating BigWig files")
    os.makedirs(os.path.join(work_dir,
                             "bigwig"),
                exist_ok = True)
    
    file_list = glob.glob(os.path.join(work_dir,
                                       "bam",
                                       "*.bam"))
    
    #index bam files
    #first check if bam files have been indexed
    index_file_list = [bam + ".bai" for bam in file_list]
    
    for bai, bam in zip(index_file_list, file_list):
        if not file_exists(bai):
            pysam.index("-@", str(threads), bam)
    
    
    
    
    #create BigWigs with deeptools
    for bam in file_list:
        if "-sort-bl.bam" in bam:
            bigwig_output = bam.replace("-sort-bl.bam", "-norm.bw")
            bigwig_output = bigwig_output.replace("bam", "bigwig")

            bigwig = "bamCoverage -p " + str(threads) + " --binSize 200 --normalizeUsing RPKM --extendReads 200 --effectiveGenomeSize 2827437033 -b " 
            bigwig = bigwig + bam +" -o " + bigwig_output
        
            subprocess.run(bigwig, 
                           shell = True)
        elif "-sort-bl-dedupl.bam" in bam:
            bigwig_output = bam.replace("-sort-bl-dedupl.bam", "-dedupl-norm.bw")
            bigwig_output = bigwig_output.replace("bam", "bigwig")
            
            bigwig = "bamCoverage -p " + str(threads) + " --binSize 200 --normalizeUsing RPKM --extendReads 200 --effectiveGenomeSize 2827437033 -b " 
            bigwig = bigwig + bam +" -o " + bigwig_output
        
            subprocess.run(bigwig, 
                           shell = True)


def bigwigQC(work_dir, threads):
    
    #generate PCA plot of all BigWig files
    if os.path.exists(os.path.join(work_dir,
                                   "bigwig")):
        file_list = glob.glob(os.path.join(work_dir,
                                           "bigwig",
                                           "*.bw"))
        
        summary_file = os.path.join(work_dir,
                                    "bigwig",
                                    "bigwig-summary.npz")
        summary = "multiBigwigSummary bins --numberOfProcessors " + str(threads)
        summary =  summary + " -b " + " ".join(file_list) + " -o " + summary_file
        
        subprocess.run(summary, 
                           shell = True)
            
def checkFastqc(script_dir):
    
    #Check for FastQC in $PATH
    path = os.environ["PATH"].lower()
    
    if "fastqc" not in path:
        #Check for FastQC elsewhere
        fastqc = [line[0:] for line in subprocess.check_output("find $HOME -name run_fastqc.bat", 
                                                               shell = True).splitlines()]
        try:
            fastqc = fastqc[0].decode("utf-8")
            fastqc_file = os.path.join(os.path.dirname(fastqc), "fastqc")
            return(fastqc_file)
        except:
            print("WARNING: FastQC was not found\nInstalling FastQC now")
            url = "https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip"
            download_file = os.path.join(script_dir,"fastqc_v0.11.9.zip")
            urllib.request.urlretrieve(url,download_file)
            #unzip FastQC file
            with ZipFile(download_file, 'r') as zip_ref:
                zip_ref.extractall(script_dir)
            #make fastqc file executable by shell
            fastqc_file = os.path.join(script_dir,"FastQC","fastqc")
            st = os.stat(fastqc_file)
            os.chmod(fastqc_file, st.st_mode | stat.S_IEXEC)
            #remove download file
            os.remove(download_file)
    else:
        fastqc_file = "fastqc"
        return(fastqc_file)

def fastqc(script_dir, work_dir, threads, file_extension):
    
    #check for FastQC
    fastqc_file = checkFastqc(script_dir)
        
    #run fastqc/multiqc 
    if not os.path.isdir(os.path.join(work_dir,"fastqc")) or len(os.listdir(os.path.join(work_dir,"fastqc"))) == 0:
        os.makedirs(os.path.join(work_dir,"fastqc"),exist_ok=True)
        fastqc_command = fastqc_file + " --threads " + str(threads) + " --quiet -o fastqc/ raw-data/*" + file_extension
        multiqc_command = ["multiqc","-o","fastqc/","fastqc/"]
        #log commands
        with open(os.path.join(work_dir,"commands.log"),"w") as file:
            file.write("FastQC: ")
            print(fastqc_command, file=file)
            file.write("MultiQC: ")
            print(*multiqc_command, sep=" ", file=file)
        try:
            print("Running FastQC on fastq files")
            subprocess.run(fastqc_command, shell=True)
        except:
            sys.exit("ERROR: FastQC failed, check logs")
        print("Running MultiQC")
        subprocess.run(multiqc_command)
    else:
        print("Skipping FastQC/MultiQC (already performed)")

def getEND(work_dir):
    '''
    Determine whether samples are single-end of paired-end
    '''    
    
    file_list = glob.glob(os.path.join(work_dir,
                                       "raw-data",
                                       "*.gz")) 
    
    PE_tag = "R2_001.fastq.gz"
    PE = b_any(PE_tag in x for x in file_list)
    
    ##based on Illumina naming convention: https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/NamingConvention_FASTQ-files-swBS.htm
    
    if PE == True:
        return("PE")
    else:
        return("SE")
    
def trim(script_dir, threads, work_dir):
    
    #check for trim_galore
    path = os.environ["PATH"].lower()
    
    if "trimgalore" not in path:
        #Check for TrimGalore elsewhere
        trimgalore = [line[0:] for line in subprocess.check_output("find $HOME -name trim_galore", 
                                                               shell = True).splitlines()]
        try:
            trimgalore_file = trimgalore[0].decode("utf-8")
        except:
            print("WARNING: TrimGalore was not found\nInstalling TrimGalore now")
            url = "https://github.com/FelixKrueger/TrimGalore/archive/refs/tags/0.6.7.zip"
            download_file = os.path.join(script_dir,"TrimGalore-0.6.7.zip")
            urllib.request.urlretrieve(url,download_file)
            #unzip TrimGalore file
            with ZipFile(download_file, 'r') as zip_ref:
                zip_ref.extractall(script_dir)
            
            #remove download file
            os.remove(download_file)
    else:
        trimgalore_file = "trim_galore"
    
    #tell TrimGalore where FastQC is located
    fastqc_file = checkFastqc(script_dir)
    trimgalore_file = [line[0:] for line in subprocess.check_output("find $HOME -name trim_galore", 
                                                               shell = True).splitlines()]
    trimgalore_file = trimgalore_file[0].decode("utf-8")
    new_line = '"' + "my $path_to_fastqc = " + "q^" + fastqc_file + "^;" + '"'
    awk_command = "awk " + "'NR==456 {$0=" + new_line + "} { print }' " + trimgalore_file + " > " + trimgalore_file + "_temp"

    subprocess.run(awk_command,
                           shell = True)
    #remove original file and rename temp file to original name
    os.remove(trimgalore_file)
    os.rename(trimgalore_file + "_temp", 
              trimgalore_file)
    
    #make TrimGalore file executable by shell
    trimgalore_file = os.path.join(script_dir, 
                                   "TrimGalore-0.6.7", 
                                   "trim_galore")
    st = os.stat(trimgalore_file)
    os.chmod(trimgalore_file, st.st_mode | stat.S_IEXEC)
    
    #cap threads at 4 for trim_galore
    if int(threads) > 4:
        threads = "4"

      
    def trimPE(work_dir, threads):
        print("Trimming paired-end fastq files")
        fastq_list = glob.glob(os.path.join(work_dir,"raw-data","*R1_001.fastq.gz"))
        for read1 in fastq_list:
            out_dir = os.path.dirname(read1)
            out_dir = out_dir.replace("raw-data","trim")
            out_file1 = read1.split(".",1)[0] + "_val_1.fq.gz"
            out_file1 = os.path.basename(out_file1)
            out_file1 = os.path.join(out_dir, out_file1)
            if not file_exists(out_file1):
                read2 = read1.replace("R1","R2")
                trim_galore_command = [trimgalore_file,"-j", threads, "-o", 
                                       "./trim", "--paired", read1, read2]
                #log commands
                with open(work_dir+"/commands.log", "a") as file:
                    file.write("Trim Galore: ")
                    print(*trim_galore_command, sep = " ",file=file)
                subprocess.run(trim_galore_command)
                
    def trimSE(work_dir, threads):
        print("Trimming single-end fastq files")
        fastq_list = glob.glob(os.path.join(work_dir,
                                            "raw-data",
                                            "*.fastq.gz"))
        
        for file in fastq_list:
            out_file = os.path.join(work_dir,
                                    "trim_galore",
                                    file.replace(".fastq.gz", "_trimmed.fq.gz"))
            if not file_exists(out_file):
                trim_galore_command = [trimgalore_file, "-j", threads, 
                                       "-o", "./trim_galore",                                                                                                                                                                                                                                                                                                
                                       file]
                #log command
                with open(os.path.join(work_dir,"commands.log"), "a") as file:
                    file.write("Trim Galore: ")
                    print(*trim_galore_command, sep = " ", file=file)
                try:
                    subprocess.run(trim_galore_command)
                except:
                    sys.exit("ERROR: trimming error. Check commands log.")
    
    #Run appropriate trim function
    if getEND(work_dir) == "PE":
        trimPE(work_dir, threads)
    elif getEND(work_dir) == "SE":
        trimSE(work_dir, threads)
       
        
        