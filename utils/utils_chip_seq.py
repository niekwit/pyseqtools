#!/usr/bin/env python3

import glob
import os
import subprocess
import urllib.request
import sys
import gzip
import shutil
from builtins import any as b_any
from itertools import compress

import yaml
import pandas as pd

try:
    import pybedtools
except ModuleNotFoundError:
    pass

try:
    import pysam
except ModuleNotFoundError:
    pass

from clint.textui import colored, puts

script_dir = os.path.abspath(os.path.dirname(__file__))
script_dir = os.path.dirname(script_dir)
sys.path.append(os.path.join(script_dir, "utils"))
import utils_general as utils


###CHIP-SEQ ANALYSIS SPECIFIC FUNCTIONS
 

def hisat2(script_dir, work_dir, threads, chip_seq_settings, genome):
    
    '''
    Alignment of fastq files using HISAT2 for ChIP-Seq
    '''
    ###check for HISAT2###
    path = os.environ["PATH"].lower()
    
    if "hisat2" not in path:
        #Check for HISAT2 elsewhere
        try:
            hisat2 = [line[0:] for line in subprocess.check_output("find $HOME -name hisat2 ! -path '*/multiqc*'",
                                                                   shell = True).splitlines()]
            hisat2 = hisat2[0].decode("utf-8")
        except:
            print("WARNING: HISAT2 was not found\nInstalling HISAT2 now")
            url = "https://cloud.biohpc.swmed.edu/index.php/s/hisat2-220-Linux_x86_64/download"
            download_file = os.path.join(script_dir,"hisat2-2.2.0-Linux_x86_64.zip")
            urllib.request.urlretrieve(url,download_file)
            
            #unzip HISAT2 file
            subprocess.run(["unzip", download_file, "-d", script_dir])
            
           
            
            hisat2 = os.path.join(script_dir, 
                                  "hisat2-2.2.0", 
                                  "hisat2")
            #remove download file
            os.remove(download_file)
            
    else:
        hisat2 = "hisat2"
    
    ###check for HISAT2 index###
    index = chip_seq_settings["hisat2"][genome]
    
    def indexFromFasta(script_dir, hisat2, genome, fasta):
        #build HISAT2 index from genome fasta
        index_location = os.path.join(script_dir,
                                 "index",
                                 "hisat2",
                                 genome,
                                 "index")
        os.makedirs(os.path.join(script_dir,
                                 "index",
                                 "hisat2",
                                 genome), 
                    exist_ok = True)
        
        #get full location of hisat2 to make hisat2-build command
        if hisat2 == "hisat2":
            path = os.environ["PATH"]
            path = path.split(":")
            for i in path:
                if "hisat2" in i:
                    hisat2_dir = i
        
            hisat2_build = os.path.join(hisat2_dir,"hisat2-build")
        else:
            hisat2_build = hisat2 + "-build"
            
        build_command = "python3 " + hisat2_build + " " + fasta + " " + index_location
                    
        utils.write2log(work_dir, build_command, "HISAT2 build index: ")
        subprocess.run(build_command, shell = True)
        
        #add index to chip-seq.yaml
        with open(os.path.join(script_dir, "yaml", "chip-seq.yaml")) as f:
                doc = yaml.safe_load(f)
                
        doc["hisat2"][genome] = index_location
        with open(os.path.join(script_dir,"yaml" ,"chip-seq.yaml"), "w") as f:
            yaml.dump(doc,f)
            
        return(index_location)

    if index == "":
        if chip_seq_settings["fasta"][genome] == "":
            try:
                print("WARNING: HISAT2 index was not found")
                print("WARNING: no " + genome + " fasta file was found")
                print("Downloading pre-build index for " + genome)
                
                if genome == "hg19":
                    url = "https://genome-idx.s3.amazonaws.com/hisat/hg19_genome.tar.gz"
                elif genome == "hg38":
                    url = "https://genome-idx.s3.amazonaws.com/hisat/hg38_genome.tar.gz"
                elif genome == "mm10":
                    url = "https://genome-idx.s3.amazonaws.com/hisat/mm10_genome.tar.gz"
                else:
                    pass
                
                download_file = os.path.join(script_dir, 
                                             os.path.basename(url))
                
                if not utils.file_exists(download_file):
                    urllib.request.urlretrieve(url, 
                                               download_file)
                
                index_dir = os.path.join(script_dir, 
                                         "index", 
                                         "hisat2",
                                         genome)
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
                with open(os.path.join(script_dir,"yaml" ,"chip-seq.yaml"), "w") as f:
                    yaml.dump(doc,f)
                
                #remove download file
                os.remove(download_file)
                
            except: #backup method in case index files are offline
                print("WARNING: genome index not available from genome-idx.s3.amazonaws.com")
                                
                def downloadFasta(script_dir, chromosomes, genome):
                    
                    base_path = "ftp://hgdownload.cse.ucsc.edu/goldenPath"
                    chromosome_dir = "chromosomes"
                    extension = ".fa.gz"
                    
                    #prepare list with urls for each chromosome fasta file
                    download_file_list = []
                    for i in range(1, chromosomes):
                        file = os.path.join(base_path, 
                                            genome,
                                            chromosome_dir,
                                            "chr" + str(i) + extension)
                        download_file_list.append(file)
                    download_file_list.append(os.path.join(base_path, 
                                                  genome,
                                                  chromosome_dir,
                                                  "chr" + "X" + extension))
                    download_file_list.append(os.path.join(base_path, 
                                                  genome,
                                                  chromosome_dir,
                                                  "chr" + "Y" + extension))
                    
                    os.makedirs(os.path.join(script_dir, 
                                             "fasta", 
                                             genome), exist_ok = True)
                    out_put_list = [os.path.join(script_dir, 
                                                 "fasta", 
                                                 genome, 
                                                 os.path.basename(i)) for i in download_file_list]
                    #download fasta files
                    print("Downloading fasta files from UCSC needed for building HISAT2 index")
                  
                    for i,j in zip(download_file_list, out_put_list):
                        urllib.request.urlretrieve(i, j)
                        
                    #concatenate fasta files to build genome fasta
                    print("Building whole genome fasta file")
                    ucsc_fasta = os.path.join(script_dir, "fasta", genome, "ucsc." + genome + ".fasta")
                    zcat_command = "zcat " + " ".join(out_put_list) + " > " + ucsc_fasta
                    utils.write2log(work_dir,zcat_command,"Build fasta: ")
                    subprocess.run(zcat_command, shell = True)
                    
                    #remove downloaded fasta files
                    for i in out_put_list:
                        os.remove(i)
                    
                    #add genome fasta file to chip-seq.yaml
                    with open(os.path.join(script_dir, "yaml", "chip-seq.yaml")) as f:
                        doc = yaml.safe_load(f)
                    
                    doc["fasta"][genome] = ucsc_fasta
                    
                    with open(os.path.join(script_dir,"yaml" ,"chip-seq.yaml"), "w") as f:
                        yaml.dump(doc,f)
                    
                    return(ucsc_fasta)
           
            if genome == "hg19":
                ucsc_fasta = downloadFasta(script_dir, 23, "hg19")
            elif genome == "hg38":
                ucsc_fasta = downloadFasta(script_dir, 23, "hg38")
            elif genome == "mm9":
                ucsc_fasta = downloadFasta(script_dir, 19, "mm9")
            
            #create index
            index_location = indexFromFasta(script_dir, hisat2, genome, ucsc_fasta)
            
            #add index location to chip-seq.yaml
            with open(os.path.join(script_dir, "yaml", "chip-seq.yaml")) as f:
                    doc = yaml.safe_load(f)
                
            doc["hisat2"][genome] = index_location
                
            with open(os.path.join(script_dir,"yaml" ,"chip-seq.yaml"), "w") as f:
                yaml.dump(doc,f)
            
    elif os.path.isfile(chip_seq_settings["fasta"][genome]): 
            fasta = chip_seq_settings["fasta"][genome]
            
            #create index from fasta
            index_location = indexFromFasta(script_dir, hisat2, genome, fasta)
            
            #add index location to chip-seq.yaml
            with open(os.path.join(script_dir, "yaml", "chip-seq.yaml")) as f:
                    doc = yaml.safe_load(f)
                
            doc["hisat2"][genome] = index_location
                
            with open(os.path.join(script_dir,"yaml" ,"chip-seq.yaml"), "w") as f:
                yaml.dump(doc,f)
    else:
        index_location = index
    
    #load blacklist
    blacklist = utils.blackList(script_dir, genome)
    
    
    ###perform alignment with HISAT2###
    os.makedirs(os.path.join(work_dir, "bam"), exist_ok = True)
    
    def alignSE(work_dir, hisat2, blacklist, index_location, threads, genome):
        file_list = glob.glob(os.path.join(work_dir, "trim","*trimmed.fq.gz"))
        
        print("Generating BAM files with HISAT2 (single-end mode)")
        
        for file in file_list:
            
            if genome == "dm3-ic":
                hisat2_output = os.path.basename(file).replace("_trimmed.fq.gz",
                                                               "-dm3ic-sort-bl.bam")
            else:
                hisat2_output = os.path.basename(file).replace("_trimmed.fq.gz",
                                                               "-sort-bl.bam")
                
            hisat2_output = os.path.join(work_dir,"bam",hisat2_output)
            
            if not utils.file_exists(hisat2_output):
                samtools = utils.checkSamtools(script_dir)
                bedtools = utils.checkBedtools(script_dir)
                
                
                align_command = "zcat " + file + " | " + hisat2 + " -p " + str(threads) + " -x " + index_location + " - 2>> align.log | " + samtools + " view -q 15 -F 260 -bS -@ " + str(threads) + " - | " + bedtools + " intersect -v -a 'stdin' -b " + blacklist + " -nonamecheck | " + samtools + " sort -@ " + str(threads) + " - > " + hisat2_output
                
                utils.write2log(work_dir, align_command, "HISAT2: ")
                
                print(os.path.basename(file) + ":", file = open("align.log", "a"))
                print("Aligning " + os.path.basename(file))
                subprocess.run(align_command,
                           shell = True)


    def alignPE(work_dir, hisat2, blacklist, index_location, threads, genome):
        file_list = glob.glob(os.path.join(work_dir, "trim","*1_val_1.fq.gz"))
        
        print("Generating BAM files with HISAT2 (paired-end mode)")
        samtools = utils.checkSamtools(script_dir)
        bedtools = utils.checkBedtools(script_dir)
        
        for read1 in file_list:
            read2 = read1.replace("R1_001_val_1.fq.gz", 
                                  "R2_001_val_2.fq.gz")
            
            
            if genome == "dm3-ic":
                hisat2_output = os.path.basename(read1).replace("_R1_001_val_1.fq.gz",
                                                               "-dm3ic-sort-bl.bam")
            else:
                hisat2_output = os.path.basename(read1).replace("_R1_001_val_1.fq.gz",
                                                           "-sort-bl.bam")
            hisat2_output = os.path.join(work_dir,
                                         "bam",
                                         hisat2_output)
            
            if not utils.file_exists(hisat2_output):
                
                                
                align_command = hisat2 + " -p " + str(threads) + " -x " + index_location
                align_command = align_command + " -1 " + read1 + " -2 " + read2
                align_command = align_command + " 2>> align.log | " + samtools + " view -q 15 -F 260 -bS -@ "
                align_command = align_command + threads + " - | " + bedtools + " intersect -v -a 'stdin' -b "
                align_command = align_command + blacklist + " -nonamecheck | " + samtools + " sort -@ "
                align_command = align_command + threads + " - > " + hisat2_output
                
                utils.write2log(work_dir, align_command, "HISAT2 PE: ")
                
                print(os.path.basename(read1.replace("_R1_001_val_1.fq.gz","")) + ":", file = open("align.log", "a"))
                subprocess.run(align_command, shell = True)
                
           
    if utils.getEND(work_dir) == "PE":
        alignPE(work_dir, hisat2, blacklist, index_location, threads, genome)
    elif utils.getEND(work_dir) == "SE":
        alignSE(work_dir, hisat2, blacklist, index_location, threads, genome)

        #plot alignment rate
        
        
        #plot number of reads before deduplication
        
        
def hisat2SLURM(script_dir, work_dir, threads, chip_seq_settings, genome):
    
    puts(colored.green(f"Aligning quality trimmed files to {genome} using HISAT2"))

    #it is assumed that all data is paired-end, and hisat2/samtools/bedtools are in $PATH
    read1_list = glob.glob(os.path.join(work_dir, "trim","*R1_001_val_1.fq.gz"))
    if len(read1_list) == 0:
        puts(colored.red("ERROR: no fastq files found in trim/"))
    
    #create BAM directory
    os.makedirs(os.path.join(work_dir,"bam", genome), exist_ok = True)
    
    #load HISAT2 index
    index = chip_seq_settings["hisat2"][genome]
    
    #load blacklist
    blacklist = chip_seq_settings["blacklist"][genome]

    #load SLURM settings
    with open(os.path.join(script_dir,"yaml","slurm.yaml")) as file:
        slurm_settings = yaml.full_load(file)        

    threads = slurm_settings["ChIP-Seq"]["hisat2_CPU"]
    mem = slurm_settings["ChIP-Seq"]["hisat2_mem"]
    time = slurm_settings["ChIP-Seq"]["hisat2_time"]
    account = slurm_settings["groupname"]
    partition = slurm_settings["ChIP-Seq"]["partition"]

    for read1 in read1_list:
        read2 = read1.replace("R1_001_val_1.fq.gz","R2_001_val_2.fq.gz")
        out_put_file = os.path.basename(read1.replace("_R1_001_val_1.fq.gz","-sort-bl.bam"))
        out_put_file = os.path.join(work_dir, "bam", genome, out_put_file)
        
        if not utils.file_exists(out_put_file):
            command = ["hisat2","-p", threads, "-x", index, "-1", read1, "-2", read2,
                       "|", "samtools", "view", "-q", "15", "-F", "260", "-bS", "-@",
                       threads, "-", "|", "bedtools", "intersect","-v", "-a", "'stdin'",
                       "-b", blacklist, "-nonamecheck", "|", "samtools", "sort", "-@",
                       threads, "-", ">", out_put_file]
            os.makedirs(os.path.join(work_dir,"slurm"), exist_ok = True)
            csv = os.path.join(work_dir,"slurm",f"slurm_hisat2_{genome}.csv")
            csv = open(csv, "a")  
            csv.write(" ".join(command) +"\n")
            csv.close()  
            
    print(f"Generating slurm_hisat2_{genome}.sh")
    csv = os.path.join(work_dir,"slurm",f"slurm_hisat2_{genome}.csv")
    commands = int(subprocess.check_output(f"cat {csv} | wc -l", shell = True).decode("utf-8"))
    script_ = os.path.join(work_dir,"slurm",f"slurm_hisat2_{genome}.sh")
    script = open(script_, "w")  
    script.write("#!/bin/bash" + "\n")
    script.write("\n")
    script.write(f"#SBATCH -A {account}\n")
    script.write("#SBATCH --mail-type=BEGIN,FAIL,END" + "\n")
    script.write(f"#SBATCH -p {partition}\n")
    script.write(f"#SBATCH -D {work_dir}\n")
    script.write("#SBATCH -o slurm/slurm_hisat2_%a.log" + "\n")
    script.write(f"#SBATCH -c {threads}\n")
    script.write(f"#SBATCH -t {time}\n")
    script.write(f"#SBATCH --mem={mem}\n")
    script.write("#SBATCH -J hisat2\n")
    script.write(f"#SBATCH -a  1-{str(commands)}\n")
    script.write("\n")
    script.write("sed -n ${SLURM_ARRAY_TASK_ID}p " + f"{csv} | bash\n")
    script.close()
    
    print("Submitting SLURM script to cluster")
    job_id_hisat2 = subprocess.check_output(f"sbatch {script_} | cut -d ' ' -f 4", shell = True)
    job_id_hisat2 = job_id_hisat2.decode("UTF-8").replace("\n","")
    print(f"SLURM job submitted successfully (job ID {job_id_hisat2})")       
    

def downsample(script_dir, work_dir, threads, slurm=False):
    '''
    Downsampling of BAM files using samtools
    '''
    #check if bam files have been deduplicated
    file_list = glob.glob(os.path.join(work_dir, "bam", "*.bam"))
    dedup = b_any("dedupl" in x for x in file_list)
    
    if dedup == True:
        file_list = sorted(glob.glob(os.path.join(work_dir,"bam","*-sort-bl-dedupl.bam")))
    
    if slurm == False:
        #get counts from bam files
        df = pd.DataFrame(columns = file_list)
        
        for bam in file_list:
              count = pysam.view("-@", str(threads) ,"-c", "-F" "260", bam)
              df.loc[1, bam] = count.replace("\n","")
        
        #get lowest read count
        df.loc[1] = df.loc[1].astype(int)
        lowest_count = df.min(1).item()
        df.loc[2] = lowest_count
        
        #calculate scaling factor for each bam file
        df.loc[3] = df.loc[2] / df.loc[1]
        df.drop([1,2],0, inplace=True)
        
        #downscale bam files
        for bam,factor in df.iteritems():
            factor = factor.item()
            outfile = bam.replace(".bam","-downscaled.bam")
            if factor == 1.0:
                if not utils.file_exists(outfile):
                    shutil.copyfile(bam, outfile)
            else:
                if not utils.file_exists(outfile):
                    samtools = utils.checkSamtools(script_dir)
                    #seed for random number generator will be 0
                    command = [samtools, "view", "-@", threads,"-bs", str(factor), bam, ">", outfile]
                    utils.write2log(work_dir, " ".join(command), "")
                    subprocess.run(command)  
        
        #export scaling factors
        names = list(df.columns)
        names = [os.path.basename(x) for x in names]
        df.columns = names
        
        df.to_csv(os.path.join(work_dir,"BAM_scaling_factors.txt"), 
                  sep = "\t",
                  index = False)
    else: 
        pass
          
    

    
def dmSpikeIn(work_dir, threads):
    file_list = glob.glob(os.path.join(work_dir,"bam","*-downscaled.bam"))
    dm_file_list = [x.replace("-downscaled.bam","-dm3ic-sort-bl.bam") for x in file_list]
    
    if len(file_list) == 0:
        puts(colored.orange("WARNING: no downscaled bam files found\nDo you want to continue?"))
        
        

def ngsplot(work_dir, genome, feature, window):
    
    file_list = glob.glob(os.path.join(work_dir, "bam", "*.bam"))
    
    os.makedirs(os.path.join(work_dir, "ngsplot"), exist_ok = True)
    
    
    def ngsplotFunction(work_dir, genome, feature, window, extension, bam):
        base_name = os.path.basename(bam).replace(extension, "")
        ngsplot_dir = os.path.join(os.path.dirname(os.path.dirname(bam)), "ngsplot" , base_name)
        os.makedirs(ngsplot_dir, exist_ok = True)
        ngsplot_output = os.path.join(ngsplot_dir, base_name) + "_" + feature
        
        ngsplot = "ngs.plot.r -G " + genome + " -R " + feature + " -C " + bam + \
            " -O " + ngsplot_output + "_"+ feature + " -T " + base_name + " -L " + window
        
        utils.write2log(work_dir, ngsplot, "ngsplot: ")
         
        #subprocess.run(ngsplot,shell = True)
        
    
    for bam in file_list:
        if "-sort-bl-dedupl-downscaled.bam" in bam:
            ngsplotFunction(work_dir, genome, feature, window, "-sort-bl-dedupl-downscaled.bam", bam)
        
            
        elif "-sort-bl-dedupl.bam" in bam:
            ngsplotFunction(work_dir, genome, feature, window, "-sort-bl-dedupl.bam")
        elif "-sort-bl.bam" in bam:
            ngsplotFunction(work_dir, genome, feature, window, "-sort-bl.bam")
     
    
    
    url= "https://drive.google.com/file/d/0B5ldivL0Hd2JN05MOEFuZ0FRQTA/view?usp=sharing&resourcekey=0-Y6Sq22xOYTAb9Yng8ZlmJg"
    
    
def peak(work_dir, threads, genome, chip_seq_settings):
    
    print("Calling peaks with MACS3:")
    
    if "hg" in genome:
        _genome = "hs"
        species = "human"
    elif "mm" in genome:
        _genome = "mm"
        species = "mouse"
    elif "ce" in genome:
        _genome = "ce"
    elif "dm" in genome:
        _genome = "dm"
    
    #load/check MACS3 settings from yaml
    q_value = chip_seq_settings["MACS3"]["q-value"]
    data_format = chip_seq_settings["MACS3"]["format"]
    ip = chip_seq_settings["MACS3"]["IP"]
    broad_cutoff = chip_seq_settings["MACS3"]["broad-cutoff"]
    
    if q_value > 0.05:
        print("WARNING: selected q value higher that default (0.05)")
        q_value = str(q_value)
    else:
        q_value = str(q_value)
        
    format_list = ["ELAND", "BED", "ELANDMULTI", "ELANDEXPORT", "SAM",
                   "BAM", "BOWTIE", "BAMPE", "BEDPE", "AUTO"]    
    if data_format not in format_list:
        print("ERROR: invalid input file format chosen")
        return
    
    if data_format == "AUTO":
        print("WARNING: AUTO format tag selected. Please make sure your data is \nnot BAMPE or BEDPE, as this cannot be detected by MACS3")
    input_format = ["-f", data_format]
    
    peak_settings = ["TF","histone"]
    if ip not in peak_settings:
        print("ERROR: invalid peak calling option chosen \nAvailable options: TF or histone")
        return
    
    if broad_cutoff > 0.1:
        print("WARNING: selected broad cutoff higher that default (0.1)")
    
    if ip == "TF":
        peak_setting = "--call-summits"
    elif ip == "histone":
        peak_setting = ["--broad", "--broad-cutoff", str(broad_cutoff)]
        
    #check if there are any deduplicated bam files
    bam_list = sorted(glob.glob(os.path.join(work_dir, "bam", "*.bam")))
    dedup = b_any("dedupl" in x for x in bam_list)
    if dedup == True:
        bam_list = glob.glob(os.path.join(work_dir, "bam", "*sort-bl-dedupl*.bam"))
        
    #import sample info
    sample_df = pd.read_csv(os.path.join(work_dir, "peak_samples.csv"))
    conditions = set(sample_df["condition"])
    
    #run MACS3 for each sample and input combination
    for i in list(conditions):
        df = sample_df[sample_df["condition"] == i]
        replicates = max(set(df["Replicate"]))
        for j in range(1,replicates + 1):
            df_rep = df[df["Replicate"] == j]
            
            df_input = df_rep[df_rep["input"] == "yes"]
            input_sample = df_input["sample"].values[0]
            df_sample = df_rep[df_rep["input"] == "no"]
            chip_sample = df_sample["sample"].values[0]
            
            input_bam = os.path.join(work_dir,"bam",input_sample)
            sample_bam = os.path.join(work_dir,"bam",chip_sample)
                        
            #run MACS3
            out_dir = os.path.join(work_dir, "peaks", chip_sample)
            os.makedirs(out_dir, exist_ok = True)
            out_file = os.path.join(work_dir, "peaks", chip_sample, chip_sample + "_peaks.xls")
            #out_file = os.path.join(work_dir, "peaks", chip_sample, os.path.basename(chip_sample) + "_peaks.xls")
            
            if not utils.file_exists(out_file):
                macs3 = ["macs3", "callpeak", "-t", sample_bam, "-c", input_bam,
                         "-n", chip_sample, "--outdir", out_dir, "-q", str(q_value),
                         "-g", _genome]
                macs3.extend(input_format )
                macs3.extend(peak_setting) # add peak settings
                utils.write2log(work_dir, " ".join(macs3), "MACS3: ")
                subprocess.call(macs3)
                
                #convert MACS3 output to BED format for bedtools 
                bed_file = os.path.join(out_dir, chip_sample + ".bed")
                
                if not utils.file_exists(bed_file):
                    df_bed = pd.read_csv(os.path.join(out_dir, chip_sample + "_peaks.broadPeak"),
                                         sep = "\t",
                                         header = None)
                    df_bed = df_bed.drop(range(4,9), axis = 1)
                    #new_order = [0,1,2]
                    #df_bed = df_bed[df_bed.columns[new_order]]
                    df_bed.to_csv(bed_file, 
                                     sep = "\t", 
                                     index = False, 
                                     header = False)
                    
    #group replicates together in list and get bed files
    peak_dirs = sorted(glob.glob(os.path.join(work_dir,"peaks","*")))
    df_repl = pd.DataFrame(peak_dirs, columns = ["Dir"])
    df_repl["base"] = df_repl["Dir"].str.rsplit("_", n=1, expand=True)[0]
    
    unique_repl = set(df_repl["base"])
    replicates = []
  
    for i in unique_repl:
        j = glob.glob(os.path.join(work_dir,i + "*"))
        replicates.append(j)
        
    bed_list = []
    for i in replicates:
        bed = []
        repl_range = range(0, len(i))
        
        for j in repl_range:
            _bed = glob.glob(os.path.join(i[j], "*.bed"))
        
            bed.append(_bed)
        bed_list.append(bed)

    #put intersecting peaks from each replicate in new BED file
    if len(unique_repl) > 1:
        for bed in bed_list:
            a_bed, = bed[0]
            b_bed, = bed[1]

            os.makedirs(os.path.dirname(a_bed).rsplit("_", 1)[0], exist_ok = True)
            output_bed = os.path.join(os.path.dirname(a_bed).rsplit("_", 1)[0],
                                      os.path.basename(a_bed).rsplit("_", 1)[0]) + ".bed"

            if not utils.file_exists(output_bed):
                a_bed = pybedtools.BedTool(a_bed)
                b_bed = pybedtools.BedTool(b_bed)

                intersect_bed = a_bed.intersect(b_bed)
                out = intersect_bed.saveas(output_bed)

    #check for HOMER
    print("Annotating peaks with HOMER:")
    path = os.environ["PATH"].lower()

    if "homer" not in path:
        # Check for HOMER elsewhere
        try:
            homer = [line[0:] for line in subprocess.check_output("find $HOME -iname homer -type d ! -path '*/multiqc*'",
                                                                   shell = True).splitlines()]
            homer = homer[0].decode("utf-8")
            configure_file = os.path.join(homer, "configureHomer.pl")
            homer_dir = homer
            print("Please add the following line to your ~/.bashrc file and run pyseqtools again:")
            print("export PATH=" + os.path.join(script_dir,"homer","bin") + ":$PATH")
            return
        except:
            print("WARNING: HOMER was not found\nInstalling HOMER now")
            url = "http://homer.ucsd.edu/homer/configureHomer.pl"
            os.makedirs(os.path.join(script_dir, "homer"), exist_ok = True)
            configure_file = os.path.join(script_dir, "homer", "configureHomer.pl")
            urllib.request.urlretrieve(url, configure_file)
            subprocess.run(["perl", configure_file, "-install"])

            homer = os.path.join(script_dir, "homer")
            configure_file = os.path.join(homer, "configureHomer.pl")
            homer_dir = homer
            print("Please add the following line to your ~/.bashrc file and run pyseqtools again:")
            print("export PATH=" + os.path.join(script_dir,"homer","bin") + ":$PATH")
            return

    else:
        homer = ""
        path = os.environ["PATH"].split(":")
        homer_dir = list(compress(path, ["homer" in x.lower() for x in path]))
        
        if len(homer_dir) > 1:
            print("ERROR: multiple instances of HOMER found:" )
            print(homer_dir)
            return
        elif len(homer_dir) == 1:
            homer_dir, = homer_dir
            configure_file = os.path.join(os.path.dirname(homer_dir),"configureHomer.pl")

    #install HOMER package for chosen genome if unavailable
    homer_genomes = glob.glob(os.path.join(os.path.dirname(homer_dir), "data", "genomes", "*"))
    if not b_any(genome in x for x in homer_genomes):
        print("WARNING: HOMER genome files not found\nInstalling now for " +genome)
        subprocess.call(["perl", configure_file, "-install", genome])

    #annotate peaks with HOMER
    bed_list = glob.glob(os.path.join(work_dir, "peaks", "*", "*.bed"))
    #bed_list = glob.glob(os.path.join(work_dir, "peak", "*.bed"))
    
    for bed in bed_list:
        out_file = bed.replace(".bed", "_annotated-peaks.txt")
        
        if not utils.file_exists(out_file):
            #further annotate BED file (required by HOMER)
            df_bed = pd.read_csv(bed, sep = "\t", header = None)
            column_number = len(df_bed.columns)
            if column_number == 4:
                
                df_bed["empty"] = ""
                df_bed["strand"] = "+"
                df_bed.to_csv(bed, sep = "\t", index = False, header = False)
            elif column_number > 6:
                print("ERROR: too many columns found in .bed file\n" + bed)
                return
            
            #run HOMER
            command = ["perl", os.path.join(homer_dir, "annotatePeaks.pl"),
                       bed, genome, ">", out_file]
            command = " ".join(command)
            utils.write2log(work_dir, command, "Peak annotation (HOMER): " )
            print("Annotating peaks for " + os.path.basename(bed))
            subprocess.run(command, shell = True, check = True)
    
    #merge peak annotation with MACS3 output
    macs3_list = glob.glob(os.path.join(work_dir, "peak", "*_peaks.xls"))
    homer_list = [x.replace("_peaks.xls", "_annotated-peaks.txt") for x in macs3_list]
    
    for macs3, homer in zip(macs3_list, homer_list):
        out_file = homer.replace("_annotated-peaks.txt","_annotated_macs3_peaks.csv")
        if not utils.file_exists(out_file):
            df_macs3 = pd.read_csv(macs3, skiprows=(29), sep = "\t")
            df_homer = pd.read_csv(homer, sep = "\t")
            df_homer = df_homer.rename(columns = {df_homer.columns[0]: "name"})
            df_merge = pd.merge(df_macs3, df_homer[["name","Gene Name", "Gene Alias", 
                                                    "Gene Description", "Gene Type"]], 
                                on = "name", 
                                how = "left" )
            df_merge.to_csv(out_file, index = False)
   

def overlappingPeaks(work_dir):
    #load settings
    df_settings = pd.read_csv(os.path.join(work_dir, "peak_settings.txt"))
    
    for index, row in df_settings.iterrows():
        control_bed = row["control"]
        sample_bed = row["sample"]
        output_bed = os.path.join(work_dir,"peak", "overlap-" + control_bed.replace(".bed","") + "_vs_" +sample_bed)
        
        if not utils.file_exists(output_bed):
            a_bed = pybedtools.BedTool(os.path.join(work_dir, "peak", control_bed))
            b_bed = pybedtools.BedTool(os.path.join(work_dir, "peak", sample_bed))

            intersect_bed = a_bed.intersect(b_bed)
            out = intersect_bed.saveas(output_bed)
    
def bam_bwQC(work_dir, threads):
    #import sample info
    sample_df = pd.read_csv(os.path.join(work_dir, "samples.csv"))
    conditions = set(sample_df["condition"])
        
    #get bam and bigwig files
    bam_list = sorted(glob.glob(os.path.join(work_dir, "bam", "*.bam")))
    bw_list = sorted(glob.glob(os.path.join(work_dir, "bigwig", "*.bw")))
    
    if len(bw_list) == 0:
        bam_list = " ".join(bam_list)
        os.makedirs(os.path.join(work_dir, "chip-qc"), exist_ok= True)
        
        #create multiBamSummary:
        sum_file = os.path.join(work_dir, "chip-qc", "multibamsummary.npz")
        
        if not utils.file_exists(sum_file):
            
            bam_sum = "multiBamSummary bins --numberOfProcessors " + threads + " -b " + \
                bam_list + " -o " + sum_file
                
            utils.write2log(work_dir, bam_sum, "multiBamSummary: ")
             
            subprocess.run(bam_sum, shell = True)
        
        #create PCA plot
        pca_file = os.path.join(work_dir, "chip-qc", "PCA_bam.pdf")
        
        if not utils.file_exists(pca_file):
            
            pca = "plotPCA -in " + sum_file + " -o " + pca_file + " -T 'Principle component analysis'"
            
            utils.write2log(work_dir, pca, "PCA BAM files: ")
            
            subprocess.run(pca, shell = True)
        
        #create Peasrson correlation plot
        pearson_file = os.path.join(work_dir, "chip-qc", "scatterplot_PearsonCorr_bam-Scores.pdf")
        
        if not utils.file_exists(pearson_file):
            pearson = "plotCorrelation -in " + sum_file + " --corMethod pearson --skipZeros --plotTitle 'Pearson Correlation of Average Scores Per Read' --whatToPlot scatterplot -o " \
                + pearson_file + " --outFileCorMatrix PearsonCorr_bam-Scores.tab"
            utils.write2log(work_dir, bam_sum, "Pearson correlation BAM files: ")
            subprocess.run(pearson, shell = True)
            
        #create Spearman correlation heatmap
        spearman_file = os.path.join(work_dir, "chip-qc", "heatmap_SpearmanCorr_readCounts_bam.pdf")
        if not utils.file_exists(pearson_file):
            spearman = "plotCorrelation -in " + sum_file + " --corMethod spearman --skipZeros --plotTitle 'Spearman Correlation of Read Counts' --whatToPlot scatterplot --colorMap viridis --whatToPlot heatmap -o " + spearman_file \
             + " --outFileCorMatrix" + os.path.join(work_dir, "chip-qc","SpearmanCorr_readCounts_bam.tab")   
             
            utils.write2log(work_dir, spearman, "Spearman correlation BAM files: ")
            subprocess.run(spearman, shell = True)
    
    else:
        #create multiBigwigSummary:
        
        os.makedirs(os.path.join(work_dir, "chip-qc"), exist_ok= True)
        
        dedup = b_any("dedupl" in x for x in bw_list)
    
        if dedup == True:
            bw_list = sorted(glob.glob(os.path.join(work_dir, "bigwig", "*dedupl-norm.bw")))
            labels = [i.replace("-dedupl-norm.bw","") for i in bw_list]
        else:
            labels = [i.replace("-norm.bw","") for i in bw_list]
        
        markers = len(conditions) * "o s "
        all_colours = ["forestgreen","navy","red","orange","salmon","yellow",
                       "brown","black","cyan","orchid","chocolate","palegreen",
                       "silver","gray","teal","beige","skyblue","gold",]
        _colours = all_colours[0:len(conditions)]
        colours = []
        for i in _colours:
            colours.append(i)
            colours.append(i)
        
        colours = " ".join(colours)
        
        labels = [os.path.basename(i) for i in labels]
        labels = " ".join(labels)
                
        bw_list = " ".join(bw_list)
        sum_file = os.path.join(work_dir, "chip-qc", "multibigwigsummary.npz")

        if not utils.file_exists(sum_file):
            print("Generating multiBigWigSummary file")
            bw_sum = "multiBigwigSummary bins --numberOfProcessors " + threads + " -b " + \
                bw_list + " -o " + sum_file + " --labels " + labels
                
            utils.write2log(work_dir, bw_sum, "multiBigwigSummary: ")
            
            subprocess.run(bw_sum, shell = True)   
        
        #create PCA plot
        pca_file = os.path.join(work_dir, "chip-qc", "PCA_bigwig.pdf")
        
        if not utils.file_exists(pca_file):
            print("Generating PCA plot")
            
            pca = "plotPCA -in " + sum_file + " -o " + pca_file + " -T 'Principle component analysis'" \
                + " --markers " + markers + " --colors " + colours
            
            utils.write2log(work_dir, pca, "PCA BigWig files: ")
            
            subprocess.run(pca, shell = True)
    
        #create Peasrson correlation plot
        pearson_file = os.path.join(work_dir, "chip-qc", "scatterplot_PearsonCorr_bigwig-Scores.pdf")
        
        if not utils.file_exists(pearson_file):
            print("Generating Pearson correlation plot")
            pearson = "plotCorrelation -in " + sum_file + " --corMethod pearson --skipZeros --plotTitle 'Pearson Correlation of Average Scores Per Read' --whatToPlot scatterplot -o " \
                + pearson_file + " --outFileCorMatrix " + os.path.join(work_dir, "chip-qc","PearsonCorr_bigwig-Scores.tab")
            utils.write2log(work_dir, pearson, "Pearson correlation BigWig files: ")
            subprocess.run(pearson, shell = True)
            
        #create Spearman correlation heatmap
        spearman_file = os.path.join(work_dir, "chip-qc", "heatmap_SpearmanCorr_readCounts_bigwig.pdf")
        if not utils.file_exists(spearman_file):
            print("Generating Spearman correlation plot")
            spearman = "plotCorrelation -in " + sum_file + " --corMethod spearman --skipZeros --plotTitle 'Spearman Correlation of Read Counts' --whatToPlot heatmap --colorMap viridis -o " + spearman_file \
             + " --outFileCorMatrix " + os.path.join(work_dir,"chip-qc","SpearmanCorr_readCounts_bigwig.tab")
             
            utils.write2log(work_dir, spearman, "Spearman correlation BigWig files: ")
            subprocess.run(spearman, shell = True)
            
    #create bam Fingerprint
    #check if deduplicated BAM files are present
    dedup = b_any("dedupl" in x for x in bam_list) 
    
    if dedup == True:
        bam_list = sorted(glob.glob(os.path.join(work_dir, "bam", "*dedupl.bam")))
    
    out_file = os.path.join(work_dir,"chip-qc", "bam_fingerprints.pdf")
    
    if not utils.file_exists(out_file):
        print("Generating BAM fingerprint plot")
        labels = [os.path.basename(i.replace("-sort-bl-dedupl.bam","")) for i in bam_list]
        labels = " ".join(labels)
        bam_list = " ".join(bam_list)
        
        command = "plotFingerprint -b " + bam_list + " --labels " + labels \
            + " --skipZeros -p " + threads + " --plotFile " + out_file
        subprocess.run(command, shell = True)

        
def plotProfile(work_dir, chip_seq_settings, genome, threads):
    #check for gtf file
    gtf = chip_seq_settings["gtf"][genome]

    #download gtf file if not available
    if gtf == "":
        if genome == "hg19":
            url = "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz"
        elif genome == "hg38":
            url = "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz"
        elif genome == "mm10":
            url = "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M27/gencode.vM27.annotation.gtf.gz"
        elif genome == "mm9":
            url= "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M1/gencode.vM1.annotation.gtf.gz"
    
        print("WARNING: gtf file not found for " + genome + ", downloading now")    
        os.makedirs(os.path.join(script_dir, "gtf", genome), exist_ok = True)
        download_file = os.path.join(script_dir, "gtf", genome, os.path.basename(url))
    
        urllib.request.urlretrieve(url, download_file)
        
        with gzip.open(download_file, "rb") as f_in:
                    with open(download_file.replace(".gz",""), "wb") as f_out:
                        shutil.copyfileobj(f_in, f_out)
                        
        os.remove(download_file)
        
        #write gtf file location to chip-seq.yaml
        gtf = download_file.replace(".gz", "")
        
        with open(os.path.join(script_dir, "yaml", "chip-seq.yaml")) as f:
                    doc = yaml.safe_load(f)
                
        doc["gtf"][genome] = gtf
                
        with open(os.path.join(script_dir,"yaml" ,"chip-seq.yaml"), "w") as f:
            yaml.dump(doc,f)
        
    
    #check if deduplicated BigWig files are present
    file_list = glob.glob(os.path.join(work_dir, "bigwig", "*bw"))
    dedup = b_any("dedupl" in x for x in file_list)
    
    if dedup == True:
        file_list = glob.glob(os.path.join(work_dir, "bigwig", "*dedupl-norm.bw"))
        samples_label = [i.replace("-dedupl-norm.bw", "") for i in file_list]
        samples_label = [os.path.basename(i) for i in samples_label]
        samples_label = " ".join(samples_label)
        file_list = " ".join(file_list)
        
    else:
        file_list = " ".join(file_list)
        samples_label = [i.replace("-norm.bw", "") for i in file_list]
        samples_label = [os.path.basename(i) for i in samples_label]
        samples_label = " ".join(samples_label)
        if len(file_list) == 0:
            print("ERROR: no BigWig files found for metagene plot creation")
            return None
    
    #generate compute matrix needed for metagene plot generation
    print("Generating compute matrix for metagene plots")
    os.makedirs(os.path.join(work_dir, "metagene_plots"), exist_ok = True)
    matrix = os.path.join(work_dir, "metagene_plots", "compute_matrix.mar.gz")
    
    if not utils.file_exists(matrix):
        cm = "computeMatrix scale-regions -S " + file_list + " -R " + gtf +\
            " --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 -p " + \
                str(threads) + " -o " + matrix + " --samplesLabel " + samples_label
        utils.write2log(work_dir, cm, "Compute matrix generation: ")
        subprocess.run(cm, shell = True)

    #generate metagene plots
    print("Generating metagene plots")
    metagene_file = os.path.join(work_dir, "metagene_plots","metagene_plot_separate.pdf")
    
    if not utils.file_exists(metagene_file):
        if os.path.exists(matrix):
            pp = "plotProfile -m " + matrix + " -out " + metagene_file 
            utils.write2log(work_dir, pp, "Generate metagene plot (overview): ")
            subprocess.run(pp, shell = True)
    
    metagene_file = os.path.join(work_dir, "metagene_plots","metagene_plot_overlay.pdf")
    
    if not utils.file_exists(metagene_file):
        if os.path.exists(matrix):
            pp = "plotProfile -m " + matrix + " -out " + metagene_file + \
                " --perGroup --plotType=se --legendLocation upper-left" + \
                    " --samplesLabel " + samples_label
            utils.write2log(work_dir, pp, "Generate metagene plot (overlay): ")
            subprocess.run(pp, shell = True)

def bamCompare(work_dir, threads):
    file = open(os.path.join(work_dir,"config_bamcompare.txt"), "r")
    lines = file.readlines()
    count = 0
    for line in lines: #removes newline characters
        lines[count] = line.replace("\n","")
        count+=1
    
    os.makedirs(os.path.join(work_dir, "bigwig"), exist_ok = True)
    
    for line in lines:
        chip_bam = line.split(":")[0]
        input_bam = line.split(":")[1].split(",")[0]
        output_bw = os.path.join(work_dir,"bigwig",line.split(":")[1].split(",")[1] + "-log2fc.bw")
    
        command = ["bamCompare", "-p", threads, "-b1", chip_bam, "-b2", input_bam, "-o", output_bw]
        command = " ".join(command)
        
        if not utils.file_exists(output_bw):
            utils.write2log(work_dir, command, "" )
            subprocess.run(command, shell = True, check = True)
            
