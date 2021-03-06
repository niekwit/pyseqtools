#!/usr/bin/env python3

import os
import subprocess
from subprocess import DEVNULL
from subprocess import CalledProcessError
import multiprocessing
import sys
import glob
import hashlib
import urllib.request
from zipfile import ZipFile
import stat
from builtins import any as b_any
import tarfile
from shutil import copy
import gzip
import shutil
import re
from pathlib import Path



from clint.textui import colored, puts
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
try: 
    import pysam
except ModuleNotFoundError:
    pass
    
import yaml
try:
    import git #module name is GitPython
except ModuleNotFoundError:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "GitPython"])

###GENERAL FUNCTIONS


def check_whitespace(work_dir):
    if " " in work_dir:
        print("ERROR: please remove any whitespace from analysis directory name:")
        print(work_dir)
        sys.exit(1)
    

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

    if not os.path.exists(os.path.join(work_dir, "md5sums_checked.csv")):
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
                df.to_csv("md5sums_checked.csv",index=False)


def logCommandLineArgs(work_dir):
    args = sys.argv
    args = " ".join(args)

    if "-h" not in args:
            if "--help" not in args:
                print("Command line arguments: " + args,
                      file = open(os.path.join(work_dir,"commands.log"), "a"))

def write2log(work_dir,command,name=""):
    with open(os.path.join(work_dir,"commands.log"), "a") as file:
        file.write(name)
        print(*command, sep = "", file = file)


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
    file_extension = ".".join(test_file.rsplit(".",2)[-2:])
    return(file_extension)


def file_exists(file): #check if file exists/is not size zero
    '''Check if file exists or has file size >0.
    '''    

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
        samtools = [line[0:] for line in subprocess.check_output("find $HOME -name samtools", shell = True).splitlines()]

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

    else:
        return("samtools")


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
    else:
        return("bedtools")


def checkBowtie(script_dir):
    #check for bowtie1 (exclude bowtie2) in $PATH
    path = os.environ["PATH"].lower().split(":")

    bowtie = list(filter(lambda x: re.search(r"bowtie*[-1][^2]", x), path))
    if len(bowtie) == 1:
        bowtie, = bowtie #unpack list
        return(bowtie)
    elif len(bowtie) > 1:
        print("ERROR: multiple instances of Bowtie found:")
        print(bowtie)
        sys.exit()
    elif len(bowtie) == 0:
        try:
            bowtie = [line[0:] for line in subprocess.check_output('find $HOME -depth -type d -iname *bowtie* | grep "bowtie*-*1[^2]"', shell = True).splitlines()]
            if len(bowtie) > 1:
                print("ERROR: multiple instances of Bowtie found:")
                bowtie = [i.decode("utf-8") for i in bowtie]
                print(bowtie)
                sys.exit()
            elif len(bowtie) == 1:
                bowtie = bowtie[0].decode("utf-8")
                return(bowtie)
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
            return(bowtie)


def checkBowtie2(script_dir):
    #check for bowtie2 in $PATH
    path = os.environ["PATH"]
    if "bowtie2" in path:
        print("test")
        return(None)
    
    path = path.os.path.lower().split(":")
    bowtie2 = list(filter(lambda x: "bowtie2" in x, path))

    
    if len(bowtie2) > 1:
        print("ERROR: multiple instances of Bowtie found:")
        print(bowtie2)
        sys.exit()
    elif len(bowtie2) == 0:
        try:

            if sys.platform in ["linux", "linux2"]:
                bowtie2 = [line[0:] for line in subprocess.check_output("find $HOME -type d -iname bowtie2*linux* ! -path '*/multiqc*' ! -path '*/pyseqtools/index/*'", shell = True).splitlines()]
            elif sys.platform == "darwin":
                bowtie2 = [line[0:] for line in subprocess.check_output("find $HOME -type d -iname bowtie2*mac* ! -path '*/multiqc*' ! -path '*/pyseqtools/index/*'", shell = True).splitlines()]
            if len(bowtie2) > 1:
                print("ERROR: multiple instances of Bowtie2 found:")
                bowtie2 = [i.decode("utf-8") for i in bowtie2]
                print(bowtie2)
                sys.exit()
            elif len(bowtie2) == 1:
                bowtie2 = bowtie2[0].decode("utf-8")
                return(bowtie2)
        except CalledProcessError: #when no instance of bowtie is found
            print("WARING: Bowtie2 not found\nInstalling Bowtie now")
            if sys.platform in ["linux", "linux2"]:
                url = "https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.4/bowtie2-2.4.4-linux-x86_64.zip/download"
                download_file = os.path.join(script_dir, "bowtie2-2.4.4-linux-x86_64.zip")
                bowtie2 = os.path.join(script_dir, "bowtie2-2.4.4-linux-x86_64")
            elif sys.platform == "darwin":
                url = "https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.4/bowtie2-2.4.4-macos-x86_64.zip/download"
                download_file = os.path.join(script_dir, "bowtie2-2.4.4-macos-x86_64.zip")
                bowtie2 = os.path.join(script_dir, "bowtie2-2.4.4-macos-x86_64")

            #download bowtie zip file
            urllib.request.urlretrieve(url, download_file)

            #unzip bowtie
            unzip = ["unzip", "-qq", download_file, "-d", script_dir]
            subprocess.run(unzip)
            return(bowtie2)

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
        picard = [line[0:] for line in subprocess.check_output("find $HOME -name picard.jar",
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
    else:
        return("picard.jar")


def deduplicationBam(script_dir, work_dir, threads, args):
    #Get Picard location
    picard = checkPicard(script_dir)


    file_list = glob.glob(os.path.join(work_dir,"bam","*-sort-bl.bam"))
    print("Performing deduplication of BAM files")
    for bam in file_list:
        dedup_output = bam.replace("-sort-bl.bam",
                                   "-sort-bl-dedupl.bam")
        if not file_exists(dedup_output):

            log_name = os.path.basename(bam).split("-", 1)[0] + "-dedup_metrics.log"
            log_name = os.path.join(work_dir,
                                    "bam",
                                    log_name)
            picard_command = "java -jar " + picard + " MarkDuplicates INPUT=" + bam + " OUTPUT=" + dedup_output + " REMOVE_DUPLICATES=TRUE METRICS_FILE=" + log_name

            write2log(work_dir, picard_command, "Deduplication: ")

            subprocess.run(picard_command,
                           shell = True)

    #plot number of reads after deduplication
    #get counts from non-deduplicated bam files

    file_list = sorted(glob.glob(os.path.join(work_dir,
                                       "bam",
                                       "*-sort-bl.bam")))

    column_names = [os.path.basename(bam).replace("-sort-bl.bam","") for bam in file_list]
    df = pd.DataFrame(columns = column_names)

    for bam in file_list:
          count = pysam.view("-@", str(threads) ,"-c", "-F" "260", bam)
          column = os.path.basename(bam).replace("-sort-bl.bam","")
          df.loc[1, column] = count

    df["condition"] = "pre-deduplication"

    #get counts for deduplicated bam files
    file_list = sorted(glob.glob(os.path.join(work_dir,
                                       "bam",
                                       "*-sort-bl-dedupl.bam")))

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


def indexBam(work_dir, threads, genome="hg38", slurm=False, script_dir=None):
    '''
    Index BAM files in bam/ directory using samtools

    '''
    
    if slurm == False:
        puts(colored.green("Indexing BAM files"))
        file_list = glob.glob(os.path.join(work_dir,"bam","*.bam"))
        
        #directory structure for TT-Seq experiments is different
        if len(file_list) == 0:
            file_list = glob.glob(os.path.join(work_dir, "bam", genome, "*", "*_sorted.bam"))
            if len(file_list) == 0:
                return("ERROR: no bam files found to be indexed")
    
        #index bam files
        #also check if bam files have been indexed
        index_file_list = [bam + ".bai" for bam in file_list]
    
        for bai, bam in zip(index_file_list, file_list):
            if not file_exists(bai):
                print(os.path.basename(bam))
                pysam.index("-@", str(threads), bam)
    else:
        puts(colored.green("Indexing BAM files (SLURM)"))
        #find all BAM files in working directory
        file_list = list(Path(work_dir).rglob("*.bam"))
        file_list = [x.as_posix() for x in file_list]
        index_list = [x + ".bai" for x in file_list]
        
        #load slurm settings
        with open(os.path.join(script_dir,"yaml","slurm.yaml")) as file:
            slurm_settings = yaml.full_load(file)
        threads = str(slurm_settings["samtools-index_CPU"])
        mem = str(slurm_settings["samtools-index_mem"])
        slurm_time = str(slurm_settings["samtools-index_time"])
        account = slurm_settings["groupname"]
        partition = slurm_settings["partition"]
        
        #set path for csv file that contains index commands and delete pre-existing one
        csv = os.path.join(work_dir,"slurm","slurm_indexBAM.csv")
        os.remove(csv)
        
        #create csv file with all index commands 
        for bam,index in zip(file_list,index_list):
            samtools = ["samtools", "index","-@", threads, bam]
            if not file_exists(index):
                csv = open(csv, "a")  
                csv.write(" ".join(samtools) +"\n")
                csv.close()
        
        #create slurm bash script
        csv = os.path.join(work_dir,"slurm","slurm_indexBAM.csv")
        commands = int(subprocess.check_output(f"cat {csv} | wc -l", shell = True).decode("utf-8"))
        if commands == 0:
            print("All BAM files have already been indexed")
            return()
        script_ = os.path.join(work_dir,"slurm","slurm_indexBAM.sh")
        script = open(script_, "w")  
        script.write("#!/bin/bash" + "\n")
        script.write("\n")
        script.write("#SBATCH -A " + account + "\n")
        script.write("#SBATCH --mail-type=BEGIN,FAIL,END" + "\n")
        script.write("#SBATCH -p " + partition + "\n")
        script.write("#SBATCH -D " + work_dir + "\n")
        script.write("#SBATCH -o slurm/slurm_indexBAM_%a.log" + "\n")
        script.write("#SBATCH -c " + threads + "\n")
        script.write("#SBATCH -t " + slurm_time + "\n")
        script.write("#SBATCH --mem=" + mem + "\n")
        script.write("#SBATCH -J " + "index" + "\n")
        script.write("#SBATCH -a " + "1-" + str(commands) + "\n")
        script.write("\n")
        script.write("sed -n ${SLURM_ARRAY_TASK_ID}p " + csv +" | bash\n")
        script.close()
        
        #run slurm script
        script = os.path.join(work_dir,"slurm","slurm_indexBAM.sh")
        print("Submitting SLURM script to cluster")
        job_id_index = subprocess.check_output(f"sbatch {script} | cut -d ' ' -f 4", shell = True)
        print(f"SLURM job submitted {job_id_index}")

def createBigWig(work_dir, threads):
    #creates BigWig files for all existing BAM files
    print("Generating BigWig files")
    os.makedirs(os.path.join(work_dir,"bigwig"), exist_ok = True)

    file_list = glob.glob(os.path.join(work_dir,"bam","*.bam"))

    
    for bam in file_list:
        bigwig_output = os.path.basename(bam.replace(".bam", "-norm.bw"))
        base_dir = os.path.dirname(os.path.dirname(bam))
        os.makedirs(os.path.join(base_dir,"bigwig"), exist_ok = True)
        bigwig_output = os.path.join(base_dir,"bigwig", bigwig_output)
    
        if not file_exists(bigwig_output):
            bigwig = "bamCoverage -p " + str(threads) + " --binSize 100 --normalizeUsing RPKM --extendReads 200 --effectiveGenomeSize 2827437033 -b "
            bigwig = bigwig + bam +" -o " + bigwig_output
            write2log(work_dir,bigwig, "Create BigWig file: ")
            subprocess.run(bigwig, shell = True)

    '''    
    if b_any("downscaled.bam" in x for x in file_list):
        file_list = glob.glob(os.path.join(work_dir,"bam","*downscaled.bam"))
        for bam in file_list:
            bamCoverage(work_dir, threads, bam)
    elif b_any("dedupl.bam" in x for x in file_list):
        file_list = glob.glob(os.path.join(work_dir,"bam","*dedupl.bam"))
        for bam in file_list:
            bamCoverage(work_dir, threads, bam)
    elif b_any("sort-bl.bam" in x for x in file_list):
        file_list = glob.glob(os.path.join(work_dir,"bam","*sort-bl.bam"))
        for bam in file_list:
            bamCoverage(work_dir, threads, bam)
    '''

def bigwigQC(work_dir, threads):

    #generate PCA plot of all BigWig files
    if os.path.exists(os.path.join(work_dir,
                                   "bigwig")):


        file_list = glob.glob(os.path.join(work_dir,
                                           "bigwig",
                                           "*-dedupl-norm.bw"))

        if len(file_list) == 0:
            file_list = glob.glob(os.path.join(work_dir,
                                           "bigwig",
                                           "*.bw"))

        summary_file = os.path.join(work_dir,
                                    "bigwig",
                                    "bigwig-summary.npz")
        summary = "multiBigwigSummary bins --numberOfProcessors " + str(threads)
        summary =  summary + " -b " + " ".join(file_list) + " -o " + summary_file

        if not file_exists(summary_file):
            subprocess.run(summary, shell = True)

        pca_output = os.path.join(work_dir, "bigwig", "PCA-bigwig.pdf")

        if not file_exists(pca_output):
            pca = "plotPCA -in " + summary_file + "-o " + pca_output + " -T PCA of BigWig files"
            subprocess.run(pca, shell = True)

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
        os.makedirs(os.path.join(work_dir,"fastqc"),exist_ok = True)
        fastqc_command = fastqc_file + " --threads " + str(threads) + " --quiet -o fastqc/ raw-data/*" + file_extension
        multiqc_command = ["multiqc","-o","fastqc/","fastqc/"]
        #log commands
        with open(os.path.join(work_dir,"commands.log"),"a") as file:
            file.write("FastQC: ")
            print(fastqc_command, file=file)
            file.write("MultiQC: ")
            print(*multiqc_command, sep=" ", file=file)
        try:
            print("Running FastQC on fastq files")
            subprocess.run(fastqc_command, shell = True)
        except:
            sys.exit("ERROR: FastQC failed, check logs")
        print("Running MultiQC")
        subprocess.run(multiqc_command)
    else:
        print("Skipping FastQC/MultiQC (already performed)")


def fastqcSLURM(work_dir):
    pass


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
    path = os.environ["PATH"]

    if "TrimGalore" not in path:
        #Check for TrimGalore elsewhere
        trimgalore = [line[0:] for line in subprocess.check_output("find $HOME -wholename *TrimGalore-*/trim_galore",
                                                               shell = True).splitlines()]
        try:
            trimgalore_file = trimgalore[0].decode("utf-8")
        except:
            puts(colored.orange("WARNING: TrimGalore was not found\nInstalling TrimGalore now"))
            url = "https://github.com/FelixKrueger/TrimGalore/archive/refs/tags/0.6.7.zip"
            download_file = os.path.join(script_dir,"TrimGalore-0.6.7.zip")
            urllib.request.urlretrieve(url,download_file)
            #unzip TrimGalore file
            with ZipFile(download_file, 'r') as zip_ref:
                zip_ref.extractall(script_dir)

            #remove download file
            os.remove(download_file)

            #tell TrimGalore where FastQC is located
            fastqc_file = checkFastqc(script_dir)
            trimgalore_file = [line[0:] for line in subprocess.check_output("find $HOME -name trim_galore",
                                                                       shell = True).splitlines()]
            trimgalore_file = trimgalore_file[0].decode("utf-8")
            new_line = '"' + "my $path_to_fastqc = " + "q^" + fastqc_file + "^;" + '"'
            awk_command = "awk " + "'NR==456 {$0=" + new_line + "} { print }' " + trimgalore_file + " > " + trimgalore_file + "_temp"

            subprocess.run(awk_command, shell = True)
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
    else:
        trimgalore_file = "trim_galore"


    #cap threads at 4 for trim_galore
    if int(threads) > 4:
        threads = "4"


    def trimPE(work_dir, threads):
        puts(colored.green("Trimming paired-end fastq files"))
        extension = get_extension(work_dir)
        fastq_list = glob.glob(os.path.join(work_dir,"raw-data","*R1_001." + extension))
        for read1 in fastq_list:
            out_dir = os.path.dirname(read1)
            out_dir = out_dir.replace("raw-data","trim")
            out_file1 = read1.split(".",1)[0] + "_val_1.fq.gz"
            out_file1 = os.path.basename(out_file1)
            out_file1 = os.path.join(out_dir, out_file1)
            if not file_exists(out_file1):
                read2 = read1.replace("R1","R2")
                trim_galore_command = [trimgalore_file,"-j", str(threads), "-o",
                                       os.path.join(work_dir,"trim"), "--paired", read1, read2]
                #log commands
                with open(work_dir+"/commands.log", "a") as file:
                    file.write("Trim Galore: ")
                    print(*trim_galore_command, sep = " ",file=file)
                subprocess.run(trim_galore_command)

    def trimSE(work_dir, threads):
        puts(colored.green("Trimming single-end fastq files"))
        extension = get_extension(work_dir)
        fastq_list = glob.glob(os.path.join(work_dir,"raw-data","*." + extension))

        for file in fastq_list:
            out_file = os.path.join(work_dir,
                                    "trim",
                                    file.replace("." + extension, "_trimmed.fq.gz"))
            out_file = out_file.replace("raw-data", "trim")
            if not file_exists(out_file):
                trim_galore_command = [trimgalore_file, "-j", threads,
                                       "-o", os.path.join(work_dir,"trim"),
                                       file]
                #log command
                with open(os.path.join(work_dir,"commands.log"), "a") as file:
                    file.write("Trim Galore: ")
                    print(*trim_galore_command, sep = " ", file=file)
                try:
                    subprocess.run(trim_galore_command)
                except:
                    puts(colored.red("ERROR: trimming error. Check commands log."))
                    sys.exit()

    #Run appropriate trim function
    if getEND(work_dir) == "PE":
        trimPE(work_dir, threads)
    elif getEND(work_dir) == "SE":
        trimSE(work_dir, threads)


def trimSLURM(script_dir, work_dir):
    """
    Creates SLURM bash script for PE-end quality trimming using Trim_galore

    """
    puts(colored.green("Trimming paired-end fastq files"))
    
    #create output directories
    extension = get_extension(work_dir)
    read1_list = glob.glob(os.path.join(work_dir,"raw-data","*R1_001." + extension))
    os.makedirs(os.path.join(work_dir, "trim"), exist_ok=True)
    os.makedirs(os.path.join(work_dir, "slurm"), exist_ok=True)
    
    #load slurm settings
    with open(os.path.join(script_dir,
                               "yaml",
                               "slurm.yaml")) as file:
            slurm_settings = yaml.full_load(file)
    
    #load slurm parameters
    threads = str(slurm_settings["TT-Seq"]["Trim_galore_CPU"])
    trim_mem = str(slurm_settings["TT-Seq"]["Trim_galore_mem"])
    trim_time = str(slurm_settings["TT-Seq"]["Trim_galore_time"])
    account = slurm_settings["groupname"]
    partition = slurm_settings["TT-Seq"]["partition"]
    #email = slurm_settings["email"]
    
    
    #write trim commands to file for slurm job array
    #trim_galore should be in $PATH
    csv = os.path.join(work_dir,"slurm","slurm_trim.csv")
    if os.path.exists(csv):
        os.remove(csv)
    
    for read1 in read1_list:
        read2 = read1.replace("R1_001." + extension,"R2_001." + extension)
        out_file1 = read1.split(".",1)[0] + "_val_1.fq.gz"
        out_file1 = out_file1.replace("raw-data", "trim")
        
        if not file_exists(out_file1):
            trim_galore = ["trim_galore","-j", str(threads), "-o",
                           os.path.join(work_dir,"trim"), "--paired", read1, read2, "\n"]
            trim_galore = " ".join(trim_galore)
            csv = os.path.join(work_dir,"slurm","slurm_trim.csv")
            csv = open(csv, "a")  
            csv.write(trim_galore)
            csv.close()
            
    #if trimming has already been done return none
    csv = os.path.join(work_dir,"slurm","slurm_trim.csv")
    if not os.path.exists(csv):
        print("Skipping trimming (already performed for all files)")
        return(None)
    
    #generate slurm bash script
    script = os.path.join(work_dir,"slurm","slurm_trim.sh")
    script = open(script, "w")  
    script.write("#!/bin/bash" + "\n")
    script.write("\n")
    script.write("#SBATCH -A " + account + "\n")
    script.write("#SBATCH --mail-type=FAIL" + "\n")
    script.write("#SBATCH --mail-type=END" + "\n")
    script.write("#SBATCH -p " + partition + "\n")
    script.write("#SBATCH -D " + work_dir + "\n")
    script.write("#SBATCH -o slurm/slurm_trim%a.log" + "\n")
    script.write("#SBATCH -c " + threads + "\n")
    script.write("#SBATCH -t " + trim_time + "\n")
    script.write("#SBATCH --mem=" + trim_mem + "\n")
    script.write("#SBATCH -J " + "Trim_galore" + "\n")
    script.write("#SBATCH -a " + "1-" + str(len(read1_list)) + "\n")
    script.write("\n")
    script.write("sed -n ${SLURM_ARRAY_TASK_ID}p slurm/slurm_trim.csv | bash\n")
    script.close()
    
    #run slurm bash script
    script = os.path.join(work_dir,"slurm","slurm_trim.sh")
    #slurm = ["sbatch", script]
    #subprocess.call(slurm)
    
    print("Submitting slurm_trim.sh to cluster")
    job_id = subprocess.check_output(f"sbatch {script} | cut -d ' ' -f 4", shell = True)
    job_id = job_id.decode("utf-8").relpace("\n","")
    print(f"Quality trimming completed (job id {job_id})")
    return(job_id)

def blackList(script_dir, genome):
    with open(os.path.join(script_dir, "yaml", "chip-seq.yaml")) as f:
        doc = yaml.safe_load(f)

    blacklist = doc["blacklist"][genome]

    if blacklist == "":
        print("WARNING: no blacklisted region BED file found")
        print("Downloading blacklist file for " + genome)

        def getBlacklist(script_dir, url, genome):
            os.makedirs(os.path.join(script_dir,
                                     "blacklist",
                                     genome),
                        exist_ok = True)

            download_file = os.path.join(script_dir,
                                         "blacklist",
                                         genome,
                                         os.path.basename(url))

            urllib.request.urlretrieve(url,
                                       download_file)

            #unzip bed file
            with gzip.open(download_file, "rb") as f_in:
                with open(download_file.replace(".gz",""), "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)


            #remove downloaded file
            os.remove(download_file)

            #write black list location to chip-seq.yaml
            blacklist = download_file.replace(".gz", "")

            with open(os.path.join(script_dir, "yaml", "chip-seq.yaml")) as f:
                        doc = yaml.safe_load(f)

            doc["blacklist"][genome] = blacklist

            with open(os.path.join(script_dir,"yaml" ,"chip-seq.yaml"), "w") as f:
                yaml.dump(doc,f)

            return(blacklist)

        if genome == "hg19":
            url = "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDukeMapabilityRegionsExcludable.bed.gz"
            blacklist = getBlacklist(script_dir, url, genome)
            return blacklist
        elif genome == "hg38":
            url = "https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz"
            blacklist = getBlacklist(script_dir, url, genome)
            return blacklist
        elif genome == "mm9":
            url = "http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/mm9-mouse/mm9-blacklist.bed.gz"
            blacklist = getBlacklist(script_dir, url, genome)
            return blacklist
    elif os.path.isfile(blacklist):
        if blacklist.endswith(".bed"):
            if os.path.getsize(blacklist) > 0:
                return blacklist
            else:
                sys.exit("ERROR: blacklist has size zero")
        else:
            sys.exit("ERROR: blacklist has invalid file format (should be .bed)")




def bwa(work_dir, script_dir, args, threads, chip_seq_settings, genome):

    #get $PATH
    path = os.environ["PATH"].lower()

    #check for bwa in $PATH
    if "bwa" not in path:
        #check for bwa in $HOME
        bwa = [line[0:] for line in subprocess.check_output("find $HOME -wholename *bwa/bwa",
                                                            shell = True).splitlines()]
        try:
            bwa = bwa[0].decode("utf-8")
            if len(bwa) > 1:
                print("ERROR: multiple copies of BWA found:")
                sys.exit(bwa)
        except:
            #install bwa from Github
            bwa_dir = os.path.join(script_dir, "bwa")
            url = "https://github.com/lh3/bwa.git"
            git.Git(script_dir).clone(url)
            os.chdir(bwa_dir)
            subprocess.run("make")
            os.chdir(work_dir)
            bwa = os.path.join(script_dir,"bwa","bwa")
    else:
        bwa = "bwa"

    #check for bwa index
    bwa_index = chip_seq_settings["bwa"][genome]
    if bwa_index == "":
        print("WARNING: BWA index not found for " + genome)
        print("Generating BWA index")
        index_dir = os.path.join(script_dir, "index", "bwa", genome)
        os.makedirs(index_dir, exist_ok = True)

        fasta = chip_seq_settings["fasta"][genome]

        copy(fasta, index_dir)
        fasta_index = os.path.join(index_dir,os.path.basename(fasta))
        os.rename(fasta_index,
                  os.path.join(os.path.dirname(fasta_index), genome + ".fasta"))

        #print([bwa, "index", "-p", index_file, fasta])
        subprocess.run([bwa, "index", os.path.join(os.path.dirname(fasta_index), genome + ".fasta")])
        os.remove(os.path.join(os.path.dirname(fasta_index), genome + ".fasta"))

        #add index to yaml
        with open(os.path.join(script_dir, "yaml" ,"chip-seq.yaml")) as f:
            doc = yaml.safe_load(f)
        doc["bwa"][genome] = os.path.join(index_dir, genome + ".fasta")
        with open(os.path.join(script_dir, "yaml","chip-seq.yaml"), "w") as f:
            yaml.dump(doc, f)

        #reload yaml
        with open(os.path.join(script_dir, "yaml" ,"chip-seq.yaml")) as f:
            chip_seq_settings = yaml.safe_load(f)



    def bwaMem(work_dir, threads, chip_seq_settings, genome):
        #check if data is paired-end
        paired_end = getEND(work_dir)

        #load other requirements for alignment
        if paired_end == "SE":
            file_list = glob.glob(os.path.join(work_dir, "trim_galore","*_trimmed.fq.gz"))
        elif paired_end == "PE":
            read1_list = glob.glob(os.path.join(work_dir, "trim_galore","*_R1_001_val_1.fq.gz"))

        index_path = chip_seq_settings["bwa"][genome]
        blacklist = blackList(script_dir, genome)

        common_command = ["2>>", "align.log", "|", "samtools", "view", "-q", "15",
                          "-F", "260", "-bS", "-@", threads, "-", "|", "bedtools",
                          "intersect", "-v", "-a", "stdin", "-b", blacklist,
                          "-nonamecheck", "|", "samtools", "sort", "-@",
                          threads, "-", ">"] #output file will be be added later

        #Run BWA
        os.makedirs(os.path.join(work_dir, "bam"), exist_ok = True)
        if paired_end == "SE":
            print("Aligning fastq files with BWA mem (single-end mode)")
            for file in file_list:
                out_file = file.replace("_trimmed.fq.gz","-sort-bl.bam")
                out_file = out_file.replace("trim_galore","bam")

                bwa_mem = ["bwa", "mem", index_path, "-t", threads] #input file(s) to be specified

                bwa_mem.append(file)
                bwa_mem.extend(common_command)
                bwa_mem.append(out_file)

                if not file_exists(out_file):
                    with open(os.path.join(work_dir, "align.log"), "a") as f:
                        print(os.path.basename(file) + ":", file = f)
                    write2log(work_dir, " ".join(bwa_mem), "BWA mem: ")

                    print("Aligning " + os.path.basename(file))
                    bwa_mem = " ".join(bwa_mem)
                    subprocess.run(bwa_mem, check = True, text = True, shell = True)
        elif paired_end == "PE":
            read2_list = [i.replace("_R1_001_val_1.fq.gz"," _R2_001_val_1.fq.gz") for i in read1_list]
            print("Aligning fastq files with BWA mem (paired-end mode)")
            for read1, read2 in zip(read1_list, read2_list):
                out_file = read1.replace("_R1_001_val_1.fq.gz","*-sort-bl.bam")
                out_file = out_file.replace("trim_galore","bam")
                bwa_mem.append(read1)
                bwa_mem.append(read2)
                bwa_mem.extend(common_command)
                bwa_mem.append(out_file)

                if not file_exists(out_file):
                    write2log(work_dir, " ".join(bwa_mem), "BWA mem: ")
                    print("Aligning " + os.path.basename(read1.replace("_R1_001_val_1.fq.gz","")))
                    subprocess.run(bwa_mem)

    def bwaAln():
        pass

    #run selected BWA aligner
    if args["align"] == "bwa-mem":
        bwaMem(work_dir, threads, chip_seq_settings, genome)
    else:
        pass


def mergeBam(work_dir, threads):
    file = open(os.path.join(work_dir,"merge-bam.csv"), "r")
    lines = file.readlines()
    lines = lines[1:]
    count = 0
    for line in lines: #removes newline characters
        lines[count] = line.replace("\n","")
        count+=1
    
    
    for line in lines:
        out_bam = os.path.join(work_dir, "bam", line.split(",")[0]) + ".bam"
        in_bams = line.split(",")[1].split(";")
        in_bams = [os.path.join(work_dir,"bam",x) for x in in_bams]
        if not file_exists(out_bam):
            try: 
                pysam.merge("-@", str(threads) ,out_bam , in_bams[0], in_bams[1], in_bams[2], in_bams[3], in_bams[4])
            except:
                pysam.merge("-@", str(threads) ,out_bam , in_bams[0], in_bams[1], in_bams[2], in_bams[3])
                
        
def bamToFastq(work_dir):
    pass
        
        
        
        
        