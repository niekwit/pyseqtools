#!/usr/bin/env python3

import os
import subprocess
import sys
import glob
from pathlib import Path

import yaml
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

script_dir = os.path.abspath(os.path.dirname(__file__))
script_dir = os.path.dirname(script_dir)
sys.path.append(os.path.join(script_dir, "utils"))
import utils_general as utils
#plt.style.use(os.path.join(script_dir,"utils", "pyseqtools.mplstyle"))


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

        #reload yaml and index
        with open(os.path.join(script_dir, "yaml", "cut-run.yaml")) as file:
            cutrun_settings = yaml.full_load(file)
        bowtie_index = cutrun_settings["bowtie"][genome]

    #check for bowtie
    bowtie_dir = utils.checkBowtie(script_dir)
    
    with open(os.path.join(script_dir, "yaml", "cut-run.yaml")) as f:
        doc = yaml.safe_load(f)

    bowtie_index = doc["bowtie"][genome]

    #check for samtools
    samtools_bin = utils.checkSamtools(script_dir)

    #check for bedtools
    bedtools_bin = utils.checkBedtools(script_dir)

    #load blacklist for selected genome
    blacklist = utils.blackList(script_dir, genome)

    #check if data is paired-end
    paired_end = utils.getEND(work_dir)

    #paired or single-end alignment:
    samtools = ["|", samtools_bin, "view", "-q", "15", "-F", "260", "-bS",
                "-@", threads, "|", bedtools_bin, "intersect", "-v", "-a",
                "stdin", "-b", blacklist, "-nonamecheck", "|", samtools_bin,
                "sort", "-@", threads, ">"] #common samtools command


    if paired_end == "SE":
        print("Aligning reads with Bowtie (single-end mode):")
        #prepare lists for read1/2 and output files
        file_list = glob.glob(os.path.join(work_dir, "trim","*_trimmed.fq.gz"))
        output_list = [x.replace("trim", "bam").replace("_trimmed.fq.gz","-sort-bl.bam") for x in file_list]

        out_dir = os.path.join(work_dir, "bam")
        os.makedirs(out_dir, exist_ok = True)

        for fq, out_file in zip(file_list,output_list):
            if not utils.file_exists(out_file):
                print("Aligning " + os.path.basename(fq).replace("_trimmed.fq.gz", ""))
                bowtie_command = ["zcat", fq, "|", os.path.join(bowtie_dir, "bowtie"),"-m", "1", "-v", "2", "-S", "-I", "0", "-X", "2000"]

                bowtie_command.extend(samtools)
                bowtie_command.append(out_file)
                print(os.path.basename(fq) + ":", file = open("align.log", "a"))
                utils.write2log(work_dir, bowtie_command, "Bowtie alignment: ")
                subprocess.call(bowtie_command)

    elif paired_end == "PE":
        print("Aligning reads with Bowtie (paired-end mode):")

        #prepare lists for read1/2 and output files
        read1_list = glob.glob(os.path.join(work_dir, "trim","*_R1_001_val_1.fq.gz"))
        read2_list = [x.replace("_R1_001_val_1.fq.gz", "_R2_001_val_2.fq.gz") for x in read1_list]
        output_list = [x.replace("trim", "bam").replace("_R1_001_val_1.fq.gz","-sort-bl.bam") for x in read1_list]

        out_dir = os.path.join(work_dir, "bam")
        os.makedirs(out_dir, exist_ok = True)

        for read1, read2, out_file in zip(read1_list, read2_list, output_list):

            if not utils.file_exists(out_file):
                print("Aligning " + os.path.basename(read1).replace("_R1_001_val_1.fq.gz", ""))
                bowtie_command = [os.path.join(bowtie_dir, "bowtie"),"-x", bowtie_index,
                                  "-1", read1, "-2", read2,"-p", threads,
                                  "-m", "1", "-v", "2", "-S", "-I", "0",
                                  "-X", "2000", "--sam-nohead","2>>",
                                  os.path.join(work_dir, "align.log")]

                bowtie_command.extend(samtools)
                bowtie_command.append(out_file)
                print(os.path.basename(read1) + ":", file = open("align.log", "a"))
                utils.write2log(work_dir, " ".join(bowtie_command), "")
                subprocess.call(bowtie_command)



def bowtie2(work_dir, script_dir, threads, cutrun_settings, genome):
    #check for bowtie2
    bowtie2_dir = utils.checkBowtie2(script_dir)
    if not bowtie2_dir:
        bowtie2_dir = ""

    #check for samtools
    samtools_bin = utils.checkSamtools(script_dir)

    #check for bedtools
    bedtools_bin = utils.checkBedtools(script_dir)

    #check for bowtie2 index
    index = cutrun_settings["bowtie2"][genome]
    if index == "":
        print("Bowtie2 index for " + genome + " not found\nCreating index now")
        fasta = cutrun_settings["fasta"][genome]
        index_dir = os.path.join(script_dir, "index", "bowtie2", genome)
        os.makedirs(index_dir, exist_ok = True)
        index = os.path.join(index_dir, genome)

        bowtie2_index = [os.path.join(bowtie2_dir, "bowtie2-build"), "--threads",
                        threads, fasta, index]

        #build index
        utils.write2log(work_dir, bowtie2_index, "Bowtie2 build index for " + genome + ": " )
        subprocess.call(bowtie2_index)

        #write index location to cut-run.yaml
        with open(os.path.join(script_dir, "yaml", "cut-run.yaml")) as f:
            doc = yaml.safe_load(f)

        doc["bowtie2"][genome] = index
        with open(os.path.join(script_dir,"yaml" ,"cut-run.yaml"), "w") as f:
            yaml.dump(doc,f)

        #reload yaml
        with open(os.path.join(script_dir, "yaml", "cut-run.yaml")) as file:
            cutrun_settings = yaml.full_load(file)
        index = cutrun_settings["bowtie2"][genome]

    #check if data is paired-end
    paired_end = utils.getEND(work_dir)

    if paired_end == "SE":
        print("Aligning reads with Bowtie (single-end mode):")
        file_list = glob.glob(os.path.join(work_dir, "trim","*_trimmed.fq.gz"))
        #to be finished

    elif paired_end == "PE":
        print("Aligning reads with Bowtie2 (paired-end mode):")
        read1_list = glob.glob(os.path.join(work_dir, "trim","*_R1_001_val_1.fq.gz"))
        read2_list = [x.replace("_R1_001_val_1.fq.gz", "_R2_001_val_2.fq.gz") for x in read1_list]
        output_list = [x.replace("trim", "bam").replace("_R1_001_val_1.fq.gz","-sort-bl.bam") for x in read1_list]

        out_dir = os.path.join(work_dir, "bam")
        os.makedirs(out_dir, exist_ok = True)

        #load blacklist for selected genome
        blacklist = utils.blackList(script_dir, genome)

        for read1, read2, out_file in zip(read1_list, read2_list, output_list):

            if not utils.file_exists(out_file):
                print("Aligning " + os.path.basename(read1).replace("_R1_001_val_1.fq.gz", ""))
                bowtie2_command = [os.path.join(bowtie2_dir, "bowtie2"),
                                   "-x", index, "-1", read1, "-2", read2, "-p", threads]

                bowtie2_settings = ['--local', '--very-sensitive-local', '--no-unal',
                                    '--no-mixed', '--no-discordant', '--phred33',
                                    '-I', '10', '-X', '700',
                                    "2>>", os.path.join(work_dir, "align.log")]

                #--dovetail --phred33 #Yuan settings
                #--local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 #Henikoff settings

                samtools = ["|", samtools_bin, "view", "-q", "15", "-F", "260", "-bS",
                            "-@", threads, "-","|", bedtools_bin, "intersect", "-v", "-a",
                            "stdin", "-b", blacklist, "-nonamecheck", "|", samtools_bin,
                            "sort", "-@", threads, ">"]

                bowtie2_command.extend(bowtie2_settings)
                bowtie2_command.extend(samtools)
                bowtie2_command.append(out_file)
                bowtie2_command = " ".join(bowtie2_command)
                print(os.path.basename(read1) + ":", file = open("align.log", "a"))
                utils.write2log(work_dir, bowtie2_command, "Bowtie2 alignment: ")
                subprocess.call(bowtie2_command, shell = True)


def histogramReadLengths(work_dir, threads):
    #check for samtools
    samtools_bin = utils.checkSamtools(script_dir)

    file_list = glob.glob(os.path.join(work_dir, "bam", "*.bam"))

    #df to store frequency of read lengths of each bam file
    df = pd.DataFrame()

    for bam in file_list:
        command = [samtools_bin, "view", "-@", threads, bam, "|", "head", "-n",
                   "1000", "|", "cut", "-f", "10", "|", "perl", "-ne",
                   "'"+r'chomp;print length($_) . "\n"' + "'", "|","sort"]
        command = " ".join(command)

        frequency = subprocess.check_output(command, shell = True)
        frequency = frequency.decode("utf-8")
        frequency = list(frequency.split("\n"))
        del frequency[-1] #removes empty value

        #add frequency of read count to df
        column_name = os.path.basename(bam)
        df[column_name] = frequency

    #create histogram of frequency of read lengths
    #sns.histplot(data = df, x = column_name, bins = 10)
    df_melt = df.melt(var_name = "File", value_name = "Read length")
    sns.displot(data = df_melt,
                x="Read length",
                binwidth = 50,
                col = "File",
                col_wrap = 4)
