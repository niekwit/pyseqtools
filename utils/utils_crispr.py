#!/usr/bin/env python3

import warnings
import os
import shutil
import subprocess
from subprocess import CalledProcessError
import urllib.request
import yaml
import sys
import glob

import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import pandas as pd
import seaborn as sns
from tqdm.auto import tqdm
import gseapy as gp

script_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.join(script_dir,"utils"))
import utils_general as utils


###CRISPR SCREEN ANALYSIS SPECIFIC FUNCTIONS

   
def csv2fasta(csv,script_dir):
    df_CSV = pd.read_csv(csv)
    line_number_fasta = len(df_CSV) * 2

    df = pd.DataFrame(columns = ["column"],
                    index = np.arange(line_number_fasta))

    #create fasta df
    df["column"] = df_CSV.stack().reset_index(drop = True)
    df.iloc[0::2, :] = ">"+df.iloc[0::2, :]

    library_name = os.path.basename(csv)
    library_name = library_name.replace(".csv","")
    fasta_base = library_name + ".fasta"
    fasta_file = os.path.join(script_dir,"index",library_name,fasta_base)
    os.makedirs(os.path.join(script_dir,"index",library_name),
                exist_ok = True)
    df.to_csv(fasta_file,index = False, 
              header = False)

    #add new CRISPR library and fasta file location to library.yaml
    yaml_list = ["clip_seq","fasta","index_path","read_mod","sg_length","species"]
    with open(os.path.join(script_dir, "yaml" ,"crispr-library.yaml")) as f:
        doc = yaml.safe_load(f)
        doc[library_name] = {}
        for i in yaml_list:
            doc[library_name][i] = ""
        doc[library_name]["fasta"] = fasta_file
        with open(os.path.join(script_dir, "yaml" , "crispr-library.yaml"), "w") as f:
            yaml.dump(doc,f)

    #exit message
    sys.exit("Fasta file created and added to library.yaml\nPlease provide more CRISPR library information in this file before first run.")


def check_index(library, crispr_library, script_dir, work_dir):
    
    #Check for Bowtie2 in $PATH
    path = os.environ["PATH"].lower()
    
    if "bowtie2" not in path:
        #Check for Bowtie2
        bowtie2 = [line[0:] for line in subprocess.check_output("find $HOME -name bowtie2-build", 
                                                               shell = True).splitlines()]
        try:
            bowtie2 = bowtie2[0].decode("utf-8")
            bowtie2_dir = os.path.dirname(bowtie2)
        except IndexError:
            print("WARNING: Bowtie2 was not found\nInstalling Bowtie2 now")
            if sys.platform in ["linux", "linux2"]:
                url="https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.3/bowtie2-2.4.3-linux-x86_64.zip/download"
                download_file = os.path.join(script_dir,
                                             "bowtie2-2.4.3-linux-x86_64.zip")
            elif sys.platform == "darwin":
                url = "https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.3/bowtie2-2.4.3-macos-x86_64.zip/download"
                download_file = os.path.join(script_dir, 
                                             "bowtie2-2.4.3-macos-x86_64.zip")
            print("Installing Bowtie2")
            urllib.request.urlretrieve(url, download_file)
            #unzip Bowtie2 file
            unzip_command = "unzip -qq " + download_file+" -d " + script_dir
            subprocess.run(unzip_command,
                           shell = True)
            #remove downloaded file
            os.remove(download_file)
    else:
        bowtie2_dir = ""

                         
    try:
        index_path = library[crispr_library]["index_path"]
        fasta = library[crispr_library]["fasta"]
        print(crispr_library + " library selected")

        if index_path in ["", None]:
            print("No index file found for " + crispr_library)
            if fasta == "":
                sys.exit("ERROR: No fasta file found for " + crispr_library)
            else:
                index_base=os.path.join(script_dir,
                                "index",
                                "bowtie2",
                                crispr_library,
                                crispr_library + "-index")
                index_dir=os.path.join(script_dir, 
                                       "index", 
                                       "bowtie2",
                                       crispr_library)
                os.makedirs(index_dir,
                            exist_ok = True)
                bowtie2_build_command = os.path.join(bowtie2_dir, 
                                                     "bowtie2-build") + " " + fasta + " " + index_base
                utils.write2log(work_dir, 
                          bowtie2_build_command, 
                          "Bowtie2-build: ")
                print("Building Bowtie2 index for " + crispr_library + " library")
                try:
                    subprocess.run(bowtie2_build_command, shell=True) #build index
                    #Write bowtie2 index file location to library.yaml
                    with open(os.path.join(script_dir, "yaml" ,"crispr-library.yaml")) as f:
                        doc = yaml.safe_load(f)
                        doc[crispr_library]["index_path"] = index_base
                        with open(os.path.join(script_dir, "yaml","crispr-library.yaml"), "w") as f:
                            yaml.dump(doc, f)
                except:
                    sys.exit("ERROR: bowtie2-build failed, check logs")
    except KeyError:
        sys.exit("ERROR: CRISPR library not specified in command line")


def guide_names(library, crispr_library):
    try:
        fasta = library[crispr_library]["fasta"]
    except KeyError:
        sys.exit("ERROR: CRISPR library not specified in command line")
    output_name = fasta.replace(".fasta",
                                "-guide_names.csv")

    if not os.path.exists(output_name):
        library = pd.read_csv(fasta, names = ['guide'])
        #creates new dataframe with only guide names:
        library = library[library['guide'].str.contains('>')]
        #removes '<' character from each row:
        library = library['guide'].str.strip('>')
        #saves guide names to a .csv file:
        library.to_csv(output_name, 
                       index=False, 
                       header=False)


def count(library, 
          crispr_library, 
          mismatch, 
          threads, 
          script_dir, 
          work_dir):
    
    #reload library yaml
    with open(os.path.join(script_dir, "yaml" ,"crispr-library.yaml")) as file:
        library = yaml.full_load(file)

    os.makedirs(os.path.join(work_dir,"count"),exist_ok=True)
    try:
        read_mod = library[crispr_library]["read_mod"]
        sg_length = library[crispr_library]["sg_length"]
        sg_length = str(sg_length)
        index_path = library[crispr_library]["index_path"]
        clip_seq = library[crispr_library]["clip_seq"]
        mismatch = str(mismatch)
    except KeyError:
        sys.exit("ERROR: CRISPR library not specified in command line")

    file_extension = utils.get_extension(work_dir)

    print("Aligning reads to reference (mismatches allowed: "+mismatch+")")
    
    '''
    #bowtie2 and bash commands (common to both trim and clip)
    bowtie2 = [line[0:] for line in subprocess.check_output("find $HOME -name bowtie2-build", 
                                                               shell = True).splitlines()]

    bowtie2 = bowtie2[0].decode("utf-8")
    bowtie2_dir = os.path.dirname(bowtie2)
    
    #bowtie2 = [os.path.join(bowtie2_dir, "bowtie2"), "--no-hd", "-p", threads,
    #           "-t", "-N", mismatch, "-x", index_path, "-", "2>>", "crispr.log",
    #           "|"]   
    '''
    
    #trim
    if read_mod == "trim":
        file_list = [x for x in glob.glob(os.path.join(work_dir,"raw-data","*"+file_extension)) if not "_trim" in x]
        for file in tqdm(file_list, position = 0, leave = True):
            base_file = os.path.basename(file)
            if ".fq.gz" in file:
                trim_file = file.replace(".fq.gz","_trim.fq.gz")
            elif ".fastq.gz" in file:
                trim_file = file.replace(".fastq.gz","_trim.fastq.gz")
            
            if not utils.file_exists(trim_file):
                tqdm.write("Trimming " + base_file)
                print(base_file + ":", file = open("crispr.log", "a"))
                
                cutadapt = ["cutadapt", "-j", threads, "--quiet","--quality-base", "33",
                            "-l", sg_length, "-o", trim_file, file, "2>>", "crispr.log"]
                #cutadapt = ["cutadapt", "-j", threads,"--quality-base", "33",
                 #           "-l", sg_length, "-o", "-", file, "2>>", "crispr.log", "|"]
                
                utils.write2log(work_dir," ".join(cutadapt),"Trimming: ")
                
                try:
                    subprocess.call(" ".join(cutadapt), shell = True)
                except:
                    sys.exit("ERROR: trimming of fastq files failed, check log")
    '''
    elif read_mod == "clip":
        file_list = glob.glob(os.path.join(work_dir,"raw-data","*"+file_extension))
        for file in tqdm(file_list, position=0, leave=True):
            base_file = os.path.basename(file)
            out_file = os.path.join(work_dir,"count",base_file.replace(file_extension,".guidecounts.txt"))

            if not utils.file_exists(out_file):
                print("Aligning "+base_file)
                print(base_file+":", file=open("crispr.log", "a"))
                cutadapt = "cutadapt -j "+threads+" --quality-base 33 -a "+clip_seq+" -o - "+file+" 2>> crispr.log | "
                cutadapt = str(cutadapt)
                bowtie2 = str(bowtie2)
                count_command = cutadapt+bowtie2+bash+out_file
                utils.write2log(work_dir,count_command,"Count: ")
                try:
                    subprocess.run(count_command, 
                                   shell = True)
                except:
                    sys.exit("ERROR: read count failed, check logs")
    '''
                    
    #count reads from trimmed fastq files
    file_list = glob.glob(os.path.join(work_dir, "raw-data", "*_trim*"))
    
    for file in tqdm(file_list, position = 0, leave = True):
        base_file = os.path.basename(file).replace("_trim","")
        out_file = os.path.join(work_dir,"count",base_file.replace(file_extension,"guidecounts.txt"))
        if not utils.file_exists(out_file):
            tqdm.write("Aligning " + base_file)
            #print("Aligning " + base_file)
            print(base_file+":", file=open("crispr.log", "a"))
            hisat2 = ["zcat", file,"|", "hisat2", 
                       "--no-hd", "-p", threads, "-t", "-N", mismatch, "-x", index_path, 
                       "-", "2>>", "crispr.log", "|"]   
            bash = ["sed", "'/XS:/d'", "|", "cut", "-f3","|", "sort", "|", "uniq", 
                    "-c", "|", "sed", '"s/^ *//"', "|", "sed", "'1d'", ">" ]
            
            hisat2.extend(bash)
            hisat2.append(out_file)
            utils.write2log(work_dir," ".join(hisat2),"Count sgRNAs: ")
            
            #cutadapt.extend(bowtie2)
            #cutadapt.extend(bash)
            #cutadapt.append(out_file)
            #utils.write2log(work_dir," ".join(cutadapt),"Count: ")
                        
            try:
                subprocess.call(" ".join(hisat2), shell = True)
            except:
                sys.exit("ERROR: read counting failed, check logs")


def plot(df,y_label,save_file):
    sns.set_style("white")
    sns.set_style("ticks")
    sns.barplot(x=list(df.keys())[0],
                    y=list(df.keys())[1],
                    data=df,
                    color="royalblue",
                    edgecolor="black",
                    linewidth=1)
    plt.ylabel(y_label)
    plt.xticks(rotation = 'vertical')
    plt.xlabel("")
    plt.tight_layout()
    sns.despine()
    plt.savefig(save_file)
    plt.close()


def plot_alignment_rate(work_dir):
    plot_file=os.path.join(work_dir,"count","alignment-rate.pdf")
    if not utils.file_exists(plot_file):
        open(os.path.join(work_dir,"files.txt"),"w").writelines([ line for line in open(os.path.join(work_dir,"crispr.log")) if ".gz:" in line ])
        open(os.path.join(work_dir,"alignment-rate.txt"),"w").writelines([ line for line in open(os.path.join(work_dir,"crispr.log")) if "overall alignment rate" in line ])

        line_number=len(open(os.path.join(work_dir,"files.txt")).readlines())
        df=pd.DataFrame(columns=["file","alignment_rate"],index=np.arange(line_number))

        counter=0
        for line in open(os.path.join(work_dir,"files.txt")):
            line=line.replace(":","")
            line=line.replace("\n","")
            df.iloc[counter,0]=line
            counter+=1

        counter=0
        for line in open(os.path.join(work_dir,"alignment-rate.txt")):
            line=line.replace("% overall alignment rate","")
            line=line.replace("\n","")
            df.iloc[counter,1]=line
            counter+=1

        df["alignment_rate"]=pd.to_numeric(df["alignment_rate"])

        os.remove(os.path.join(work_dir,"files.txt"))
        os.remove(os.path.join(work_dir,"alignment-rate.txt"))

        #plot alignment rate
        plot(df,"Alignment rate (%)",plot_file)


def plot_coverage(work_dir,library,crispr_library): #plots coverage per sample after alignment
    plot_file = os.path.join(work_dir,"count","coverage.pdf")
    if not utils.file_exists(plot_file):
        #get number of sgRNAs in CRISPR library
        fasta = library[crispr_library]["fasta"]
        fasta = pd.read_table(fasta, header = None)
        lib_size = len(fasta) / 2

        #extract number of single mapped aligned reads from counts-aggregated.tsv
        df = pd.read_table(os.path.join(work_dir,"count","counts-aggregated.tsv")
                           , sep = "\t")
        
        column_names = list(df.columns)
        del column_names[0:2] #delete sgRNA and gene columns
        
        counts = {}
        
        for i in column_names:
            count_sum = []
            count_sum.append(df[i].sum())
            counts[i] = count_sum
        
        #convert counts to sequence coverage
        df = pd.DataFrame(counts)
        df = df / lib_size
        
        #order columns alphabetically
        df = df.reindex(sorted(df.columns), axis=1)
        
        #transpose data frame
        df = df.transpose()
        
        df["samples"] = df.index
        
        names = ["coverage", "samples"]
        df.columns = names
        df = df[["samples","coverage"]]
        
        #plot coverage per sample
        plot(df,"Fold sequence coverage per sample",plot_file)


def normalise(work_dir):
    df = pd.read_table(os.path.join(work_dir,"count","counts-aggregated.tsv"))
    column_range = range(2,len(df.columns))
    for i in column_range:
        column_sum = df.iloc[:,i].sum()
        df.iloc[:,i] = df.iloc[:,i] / column_sum * 1E8
        df.iloc[:,i] = df.iloc[:,i].astype(int)
    df.to_csv(os.path.join(work_dir,"count","counts-aggregated-normalised.csv"),
              index = False,
              header = True)


def join_counts(work_dir,library,crispr_library):
    file_list = glob.glob(os.path.join(work_dir,"count","*guidecounts.txt"))
    
    #dictionary to store all guide counts
    counts = {}
    
    for file in file_list:
       #sed_command = sed + file 
       #subprocess.run(sed_command, shell = True)
       #add counts to counts dict
       df = pd.read_csv(file,
                        sep=" ",
                        header = None)
       key = os.path.basename(file).replace(".guidecounts.txt", "")
       df.columns = [key, "sgRNA"]
       df[key].astype(int)
       counts.update({key:df})
       
    #prepare data frame for left join
    fasta = library[crispr_library]["fasta"]

    df = pd.read_csv(fasta, header = None)
    df.columns = ["sgRNA"]
    
    df = df[df["sgRNA"].str.contains(">")]
    df = df.reset_index(drop = True)
    df["sgRNA"] = df["sgRNA"].str.replace(">", "")
    
    df["gene"] = df["sgRNA"].str.split(pat = "_",
                                       n = 1,
                                       expand = True)[0]
    
    #perform left join
    for key, value in counts.items():
        df = pd.merge(df, value, on='sgRNA', how='left')
        
    df["sgRNA"] = df["sgRNA"].str.split(pat = "_",
                                   n = 1,
                                   expand = True)[1]
    
    #replace nan with zero
    df = df.fillna(0)
    df = df.sort_values(by = ["sgRNA"])
    
    #convert floats back to int after pandas merge (bug in pandas)
    index_range = range(2, len(df.columns))
    index_list = []
    for i in index_range:
        index_list.append(i)
        
    df.iloc[:, index_list] = df.iloc[:, index_list].astype(int)
    
    #save data frame to file
    df.to_csv(os.path.join(work_dir,"count",'counts-aggregated.tsv'), 
              sep = '\t',
              index = False)
    

def mageck(work_dir,script_dir,cnv,fdr):
    
    #Check for MAGeCK
    try:
        subprocess.check_output("which mageck", shell = True)
    except CalledProcessError:
        print("WARNING: MAGeCK was not found\nInstalling MAGeCK now")
        
        url = "https://sourceforge.net/projects/mageck/files/0.5/mageck-0.5.9.4.tar.gz/download"
        download_file = os.path.join(script_dir,"mageck-0.5.9.4.tar.gz")
        urllib.request.urlretrieve(url, download_file)
        #unpack MAGeCK file
        tar_command = "tar -xzf " + download_file
        subprocess.run(tar_command, shell = True)
        os.chdir(os.path.join(script_dir, "mageck-0.5.9.4"))
        mageck_install_command = "python3 " + os.path.join(script_dir ,
                                                           "mageck-0.5.9.4", 
                                                           "setup.py" + " install --user")
        subprocess.run(mageck_install_command, shell = True)
                
        #remove download file and folder
        os.remove(download_file)
        shutil.rmtree(os.path.join(script_dir, "mageck-0.5.9.4"))
    
        
    if fdr > 0.25:
        print("WARNING: MAGeCK FDR cut off set higher than default 0.25")

    #determine number of samples in count table
    header = subprocess.check_output(["head", "-1", os.path.join(work_dir,
                                                                 "count",
                                                                 "counts-aggregated.tsv")])
    header = header.decode("utf-8")
    sample_count = header.count("\t") - 1
    if sample_count == 2:
        if "pre" and "post" in header:
            print("Skipping MAGeCK (only CRISPR library samples present)")
            return(None)

    #check for stats.config
    stats_config = os.path.join(work_dir,"stats.config")
    if not os.path.exists(stats_config):
        print("ERROR: stats.config not found (MAGeCK comparisons)")
        return(None)

    #create MAGeCK dir
    os.makedirs(os.path.join(work_dir,"mageck"),
                exist_ok = True)

    #load MAGeCK comparisons and run MAGeCK
    try:
        df = pd.read_csv(os.path.join(work_dir,"stats.config"), 
                         sep = ";")
    except pd.errors.ParserError:
        print("ERROR: multiple semi-colons per line detected in stats.config")
        return(None)
    sample_number = len(df)
    sample_range = range(sample_number)

    def cnv_com(script_dir,cnv,ccle_ref): #generate MAGeCK command for CNV correction
        #check if specified cell line is in CCLE data list
        cell_line_list=subprocess.check_output(["head",
                                                "-1",os.path.join(script_dir,
                                                "CCLE","CCLE_copynumber_byGene_2013-12-03.txt")])
        cell_line=cnv
        cnv_command=" --cnv-norm "+ccle_ref+" --cell-line "+cell_line
        return(cnv_command,cell_line_list)

    def CCLE_cell_line_exists(script_dir,cnv,ccle_ref):
        cell_line_list=cnv_com(script_dir,cnv,ccle_ref)[1]
        cell_line_list=cell_line_list.decode("utf-8")
        cell_line=cnv
        if not cell_line in cell_line_list:
            print("ERROR: specified cell line not found in CCLE reference file")
            print("Skipping CNV correction for MAGeCK")
            return(False)
        else:
            return(True)

    for i in sample_range:
        test_sample=df.loc[i]["t"]
        control_sample=df.loc[i]["c"]
        mageck_output=test_sample+"_vs_"+control_sample

        if not utils.file_exists(os.path.join(work_dir,"mageck",mageck_output)):
            os.makedirs(os.path.join(work_dir,"mageck",mageck_output),exist_ok=True)
            prefix=os.path.join(work_dir,"mageck",mageck_output,mageck_output)
            input=os.path.join(work_dir,"count","counts-aggregated.tsv")
            log=" 2>> "+os.path.join(work_dir,"crispr.log")
            mageck_command="mageck test -k "+input+" -t "+test_sample+" -c "+control_sample+" -n "+prefix+log
            utils.write2log(work_dir,mageck_command,"MAGeCK: ")
            try:
                print("Running MAGeCK without CNV correction: " + mageck_output.replace("_vs_", " vs "))
                subprocess.run(mageck_command, shell=True)
            except:
                print("MAGeCK failed. Check log")
                return(None)

        #check if CNV correction is requested and perform checks
        if cnv != None:
            ccle_ref=os.path.join(script_dir,"CCLE","CCLE_copynumber_byGene_2013-12-03.txt")
            if not os.path.exists(ccle_ref):
                print("WARNING: no CCLE copy number file found")
                print("Downloading CCLE copy number file from https://data.broadinstitute.org")
                url=" https://data.broadinstitute.org/ccle_legacy_data/dna_copy_number/CCLE_copynumber_byGene_2013-12-03.txt"
                download_command="wget --directory-prefix="+os.path.join(script_dir,"CCLE")+url
                utils.write2log(work_dir,download_command,"Download CCLE file: ")
                try:
                    subprocess.run(download_command, shell=True)
                except:
                    sys.exit("ERROR: download failed, check log and url")

                if not CCLE_cell_line_exists(script_dir,cnv,ccle_ref):
                    return
                else:
                    cnv_command=cnv_com(script_dir,cnv,ccle_ref)[0]
            else:
                if not CCLE_cell_line_exists(script_dir,cnv,ccle_ref):
                    return
                else:
                    cnv_command=cnv_com(script_dir,cnv,ccle_ref)[0]

        if cnv != None:
            cnv_dir=os.path.join(work_dir,"mageck-cnv",mageck_output)
            if not utils.file_exists(cnv_dir):
                os.makedirs(cnv_dir, exist_ok=True)
                prefix=os.path.join(work_dir,"mageck-cnv",mageck_output,mageck_output)
                input=os.path.join(work_dir,"count","counts-aggregated.tsv")
                log=" 2>> "+os.path.join(work_dir,"crispr.log")
                mageck_command="mageck test -k "+input+" -t "+test_sample+" -c "+control_sample+" -n "+prefix+log
                mageck_command=mageck_command+cnv_command
                utils.write2log(work_dir,mageck_command,"MAGeCK: ")
                subprocess.run(mageck_command, shell=True)

    #plot MAGeCK hits
    file_list=glob.glob(os.path.join(work_dir,"mageck*","*","*gene_summary.txt"))

    for file in file_list:
        save_path=os.path.dirname(file)
        out_put_file=os.path.join(save_path,"logFC-plot.pdf")

        if not utils.file_exists(out_put_file):
            plot_command="Rscript "+os.path.join(script_dir,"R","crispr-plot-hits.R ")+ \
                        work_dir+" "+file+" mageck "+save_path+" "+mageck_output+ \
                        " "+script_dir+" "+str(fdr)
            utils.write2log(work_dir,plot_command,"Plot hits MAGeCK: ")
            try:
                subprocess.run(plot_command, shell=True)
            except:
                sys.exit("ERROR: plotting hits failed, check log")


def convert4bagel(work_dir,library,crispr_library): #convert MAGeCK formatted count table to BAGEL2 format
    count_table_bagel2=os.path.join(work_dir,"bagel",'counts-aggregated-bagel2.tsv')

    if not utils.file_exists(count_table_bagel2):
        #obtain sequences of each guide
        try:
            fasta=library[crispr_library]["fasta"]
        except KeyError:
            sys.exit("ERROR: CRISPR library not specified in command line")
        df_fasta=pd.read_csv(fasta, header=None)

        #place sgRNA name and sequence is separate columns
        df_name=df_fasta[df_fasta[0].str.contains(">")]
        names=df_name.squeeze()#convert to series
        names=names.reset_index(drop=True)#reset index
        names=names.str.replace(">","")#remove >
        names.name="sgRNA"
        df_seq=df_fasta[~df_fasta[0].str.contains(">")]
        seq=df_seq.squeeze()#convert to series
        seq=seq.reset_index(drop=True)#reset index
        seq.name="SEQUENCE"
        df_join=pd.concat([names,seq],axis=1)#create df with names and sequences

        #create gene column
        df_join["gene"]=df_join["sgRNA"].str.split(pat="_").str[0]
        #df_join.rename(columns = {"sgRNA":"gene"}, inplace = True)

        #open MAGeCK formatted count table
        count_file=os.path.join(work_dir,"count","counts-aggregated.tsv")
        df_master=pd.read_csv(count_file, sep="\t")

        #format sgRNA column df_join to same format as df_master
        df_join["sgRNA"]=df_join["sgRNA"].str.split(pat="_",n=1).str[1]

        #sort data frames
        df_join=df_join.sort_values(by="sgRNA")
        df_master=df_master.sort_values(by="sgRNA")

        #merge data frames and format data frame for BAGEL2
        df_join=df_join.drop(["gene"], axis=1)#remove gene column
        df_merge=pd.merge(df_master,df_join, on="sgRNA", how="left")
        df_merge=df_merge.drop(["sgRNA"], axis=1)#remove gene column
        df_merge.rename(columns={"gene":"GENE"}, inplace=True)
        cols=list(df_merge)
        cols.insert(0,cols.pop(cols.index("SEQUENCE")))
        df_merge=df_merge.loc[:,cols]

        #save df to file
        os.makedirs(os.path.join(work_dir,"bagel"),exist_ok=True)
        df_merge.to_csv(count_table_bagel2, sep='\t',index=False)


def remove_duplicates(work_dir):
    df = pd.read_table(os.path.join(work_dir,"bagel","counts-aggregated-bagel2.tsv"))
    df.drop_duplicates(subset = "SEQUENCE", 
                       keep = "first", 
                       inplace = True)
    df.to_csv(os.path.join(work_dir, "bagel","counts-aggregated-bagel2.tsv"),
              index = False,
              header = True,
              sep = "\t")


def bagel2(work_dir, script_dir, fdr, crispr_settings, crispr_library):
    
    #Check for BAGEL2:
    bagel2 = [line[0:] for line in subprocess.check_output("find $HOME -name BAGEL.py", 
                                                           shell = True).splitlines()]
    try: 
        bagel2_exe = bagel2[0].decode("utf-8")
        bagel2_dir = os.path.dirname(bagel2_exe)
        
    except IndexError:
        print("WARNING: BAGEL2 was not found\nInstalling BAGEL2 now")
        clone_command = "git " + "clone --quiet " + "https://github.com/hart-lab/bagel.git " + os.path.join(script_dir, "bagel2")
        subprocess.run(clone_command, 
                       shell = True)
        #removes example directory (contains old BAGEL2 script)
        shutil.rmtree(os.path.join(script_dir, 
                                   "bagel2", 
                                   "pipeline-script-example"))
        
        bagel2_dir = os.path.join(script_dir, "bagel2")
        bagel2_exe = os.path.join(script_dir, "bagel2", "BAGEL.py")
    
   
        
    #get sample names from BAGEL2 count table
    header = subprocess.check_output(["head", "-1",os.path.join(work_dir,
                                                "bagel",
                                                "counts-aggregated-bagel2.tsv")])
    header = header.decode("utf-8")
    header = header.replace("\n","")
    header = list(header.split("\t"))# convert string into list

    #create dictionary that holds column name (key) and column index
    column_dict = {key: i for i, key in enumerate(header)}
    column_dict = {key: column_dict[key] - 1 for key in column_dict} #first sample column should have value 1

    count_table = os.path.join(work_dir,
                               "bagel",
                               'counts-aggregated-bagel2.tsv')

    #reference genes files for Bayes Factor calculation
    species = crispr_settings[crispr_library]["species"]
    
    if species.lower() == "human": 
        essential_genes = os.path.join(bagel2_dir,"CEGv2.txt")
        nonessential_genes = os.path.join(bagel2_dir,"NEGv1.txt")
    elif species.lower() == "mouse":
        essential_genes = os.path.join(bagel2_dir,"CEG_mouse.txt")
        nonessential_genes = os.path.join(bagel2_dir,"NEG_mouse.txt")
        
    #load stats.config for sample comparisons
    df = pd.read_csv(os.path.join(work_dir,"stats.config"),
                   sep=";")

    #remove rows with multiple samples (for MAGeCK) per condition (temp fix)
    df = df[~df["t"].str.contains(",")]
    df = df[~df["c"].str.contains(",")]
    df = df.reset_index()

    sample_number = len(df)
    sample_range = range(sample_number)

    #run BAGEL2 for each comparison in stats.config
    for i in sample_range:
        
        test_sample = df.loc[i]["t"]
        control_sample = df.loc[i]["c"]
        
        #generate fold change table
        if "," not in test_sample and "," not in control_sample: #check for replicate comparisons
            bagel2_output = test_sample+"_vs_"+control_sample
            bagel2_output_base = os.path.join(work_dir,"bagel",bagel2_output,bagel2_output)
            control_sample_column = column_dict[control_sample]
    
            if not utils.file_exists(os.path.join(work_dir,
                                                  "bagel",
                                                  bagel2_output)):
                os.makedirs(os.path.join(work_dir,
                                         "bagel",
                                         bagel2_output),
                            exist_ok = True)
    
            print("Generatig fold change table for " + bagel2_output)
            fc_file = os.path.join(bagel2_output_base + ".foldchange")
            if not utils.file_exists(fc_file):
                bagel2fc_command = "python3 " + bagel2_exe +" fc" +" -i " + \
                                    count_table + " -o "+bagel2_output_base + \
                                    " -c " + str(control_sample_column)
                utils.write2log(work_dir,bagel2fc_command,"BAGEL2 fc: ")
                try:
                    subprocess.run(bagel2fc_command, shell = True)
                except:
                    sys.exit("ERROR: generation of BAGEL2 fc file failed, check log")
        
            #calculate Bayes factors
            print("Calculating Bayes Factors for "+bagel2_output)
            bf_file = os.path.join(bagel2_output_base+".bf")
            if not utils.file_exists(bf_file):
                #get sample names from BAGEL2 foldchange table
                header2 = subprocess.check_output(["head", "-1",os.path.join(bagel2_output_base+".foldchange")])
                header2 = header2.decode("utf-8")
                header2 = header2.replace("\n","")
                header2 = list(header2.split("\t"))# convert string into list
    
                #create dictionary that holds column name (key) and column index
                column_dict2 = {key: i for i, key in enumerate(header2)}
                column_dict2 = {key: column_dict2[key] - 1 for key in column_dict2} #first sample column should have value 1
                test_sample_column2 = column_dict2[test_sample]
    
                bagel2bf_command = "python3 " + bagel2_exe+" bf" + " -i " + fc_file + \
                                " -o " + bf_file + " -e " + essential_genes + " -n " + \
                                nonessential_genes + " -c " + str(test_sample_column2)
                utils.write2log(work_dir, bagel2bf_command, "BAGEL2 bf: ")
                try:
                    subprocess.run(bagel2bf_command, 
                                   shell = True)
                except:
                    sys.exit("ERROR: Calculation of Bayes Factors failed, check log")
            
            #calculate precision-recall
            print("Calculating precision-recall for " + bagel2_output)
            pr_file = os.path.join(bagel2_output_base + ".pr")
            if not utils.file_exists(pr_file):
                bagel2pr_command = "python3 " + bagel2_exe + " pr" + " -i " + bf_file + \
                                " -o " + pr_file + " -e " + essential_genes +" -n " + \
                                nonessential_genes
                utils.write2log(work_dir, bagel2pr_command, "BAGEL2 pr: ")
                try:
                    subprocess.run(bagel2pr_command, 
                                   shell = True)
                except:
                    sys.exit("ERROR: Calculation of precision-recall failed, check log")
        
            #plot results                    
            print("Plotting BAGEL2 results for " + bagel2_output)
            plot_file =  os.path.join(work_dir,"bagel", bagel2_output, "PR-" + bagel2_output + ".pdf")
            if not utils.file_exists(plot_file):
                plot_script = os.path.join(script_dir,
                                           "R",
                                           "crispr-plot-hits.R")
                plot_command = "Rscript "+plot_script + " " + work_dir + " " + pr_file + \
                            " bagel2 " + os.path.join(work_dir, "bagel", bagel2_output)+ \
                            " " + bagel2_output
                utils.write2log(work_dir, plot_command, "BAGEL2 plot: ")
                try:
                    subprocess.run(plot_command, shell=True)
                except:
                    sys.exit("ERROR: Calculation of precision-recall failed, check log")
                
        else: # perform BAGEL2 with replicates
            test_samples = test_sample.split(",")
            control_samples = control_sample.split(",")
            
            bagel2_output = test_sample+"_vs_"+control_sample
            bagel2_output_base = os.path.join(work_dir,"bagel",bagel2_output,bagel2_output)
            
            control_sample_columns = []
            for i in control_samples:
                control_sample_column = column_dict[i]
                control_sample_columns.append(str(control_sample_column))
            
            control_sample_columns = ",".join(control_sample_columns)
                
            if not utils.file_exists(os.path.join(work_dir,
                                                  "bagel",
                                                  bagel2_output)):
                os.makedirs(os.path.join(work_dir,
                                         "bagel",
                                         bagel2_output),
                            exist_ok = True)
                
            print("Generatig fold change table for " + bagel2_output)
            fc_file = os.path.join(bagel2_output_base + ".foldchange")
            if not utils.file_exists(fc_file):
                bagel2fc_command = "python3 " + bagel2_exe +" fc" +" -i " + \
                                    count_table + " -o "+bagel2_output_base + \
                                    " -c " + str(control_sample_columns)
                utils.write2log(work_dir,bagel2fc_command,"BAGEL2 fc: ")
                try:
                    subprocess.run(bagel2fc_command, shell = True)
                except:
                    sys.exit("ERROR: generation of BAGEL2 fc file failed, check log")
            
            #calculate Bayes factors
            print("Calculating Bayes Factors for "+bagel2_output)
            bf_file = os.path.join(bagel2_output_base+".bf")
            if not utils.file_exists(bf_file):
                #get sample names from BAGEL2 foldchange table
                header2 = subprocess.check_output(["head", "-1",os.path.join(bagel2_output_base+".foldchange")])
                header2 = header2.decode("utf-8")
                header2 = header2.replace("\n","")
                header2 = list(header2.split("\t"))# convert string into list
    
                #create dictionary that holds column name (key) and column index
                column_dict2 = {key: i for i, key in enumerate(header2)}
                column_dict2 = {key: column_dict2[key] - 1 for key in column_dict2} #first sample column should have value 1
                
                #get column numbers of test samples
                test_sample_columns = []
                for i in test_samples:
                    test_sample_column = column_dict[i]
                    test_sample_columns.append(str(test_sample_column))
                
                test_sample_columns = ",".join(test_sample_columns)
                
                #perform Bayes factor calculation
                bagel2bf_command = "python3 " + bagel2_exe+" bf" + " -i " + fc_file + \
                                " -o " + bf_file + " -e " + essential_genes + " -n " + \
                                nonessential_genes + " -c " + test_sample_columns
                utils.write2log(work_dir, bagel2bf_command, "BAGEL2 bf: ")
                try:
                    subprocess.run(bagel2bf_command, 
                                   shell = True)
                except:
                    sys.exit("ERROR: Calculation of Bayes Factors failed, check log")
            
            #calculate precision-recall
            print("Calculating precision-recall for " + bagel2_output)
            pr_file = os.path.join(bagel2_output_base + ".pr")
            if not utils.file_exists(pr_file):
                bagel2pr_command = "python3 " + bagel2_exe + " pr" + " -i " + bf_file + \
                                " -o " + pr_file + " -e " + essential_genes +" -n " + \
                                nonessential_genes
                utils.write2log(work_dir, bagel2pr_command, "BAGEL2 pr: ")
                try:
                    subprocess.run(bagel2pr_command, 
                                   shell = True)
                except:
                    sys.exit("ERROR: Calculation of precision-recall failed, check log")
            
            #plot results                    
            print("Plotting BAGEL2 results for " + bagel2_output)
            plot_file =  os.path.join(work_dir,"bagel", bagel2_output, "PR-" + bagel2_output + ".pdf")
            if not utils.file_exists(plot_file):
                plot_script = os.path.join(script_dir,
                                           "R",
                                           "crispr-plot-hits.R")
                plot_command = "Rscript "+plot_script + " " + work_dir + " " + pr_file + \
                            " bagel2 " + os.path.join(work_dir, "bagel", bagel2_output)+ \
                            " " + bagel2_output
                utils.write2log(work_dir, plot_command, "BAGEL2 plot: ")
                try:
                    subprocess.run(plot_command, shell=True)
                except:
                    sys.exit("ERROR: Calculation of precision-recall failed, check log")

    def histogramBF(df,out_file):
        sns.set_style("white")
        sns.set_style("ticks")
        sns.despine()
        sns.histplot(data=df,
                    x="BF",
                    bins=50)
        plt.xlabel('Bayes Factor')
        plt.ylabel('Number of Genes')
        plt.savefig(out_file)
        plt.clf()

    file_list=glob.glob(os.path.join(work_dir, "bagel*", "*", "*.pr"))

    for file in file_list:
        #plot histogram BF
        out_file = os.path.join(os.path.dirname(file), "histogram-BF.pdf")
        if not utils.file_exists(out_file):
            df = pd.read_table(file)
            histogramBF(df, out_file)


def lib_analysis(work_dir, crispr_settings, crispr_library, script_dir):
    #determine whether count file contains library samples pre and post
    header = subprocess.check_output(["head", "-1",
                                    os.path.join(work_dir,
                                                "count",
                                                "counts-aggregated.tsv")])
    if "pre" and "post" in str(header):
        #check if analysis has been performed already
        out_file = os.path.join(work_dir,
                                "library-analysis",
                                "normalised-guides-frequency.pdf")
        if not os.path.exists(out_file):
            print("Analysing CRISPR library quality")
        else:
            print("CRISPR library quality already analysed")
            return(None)

        os.makedirs(os.path.join(work_dir,"library-analysis"),exist_ok=True)
        warnings.filterwarnings("ignore")

        fasta = crispr_settings[crispr_library]["fasta"]
        
        subprocess.call(["Rscript", 
                   os.path.join(script_dir, "R", "crispr-library_analysis.R"), 
                   work_dir,
                   fasta])
        
        '''
        
        
        
        
        
        
        

        df = pd.read_csv(os.path.join(work_dir,
                                        "count",
                                        'counts-aggregated.tsv'),sep='\t')

        #determines total read count per column
        pre_lib_sum = df['pre'].sum()
        pre_lib_sum = int(pre_lib_sum)
        post_lib_sum = df['post'].sum()
        post_lib_sum = int(post_lib_sum)

        #normalises guide counts to total guide count
        #X: pre-amplification library, Y: post-amplification library
        datax = df['pre']
        datax = datax.sort_values(ascending=True)
        X = datax.to_numpy()
        X = X / pre_lib_sum

        datay = df['post']
        datay = datay.sort_values(ascending=True)
        Y = datay.to_numpy()
        Y = Y / pre_lib_sum

        index_len = len(df.index)

        #plots data:
        ax = sns.lineplot(x=range(index_len),
                        y=X,
                        color='navy',
                        label='Pre-amplification library')
        ax = sns.lineplot(x=range(index_len),
                        y=Y,
                        color='green',
                        label='Post-amplification library')
        ax.set_yscale('log')
        ax.legend(loc='lower right')
        ax.set(ylabel='Normalised sgRNA count', xlabel='sgRNA')
        plt.savefig(os.path.join(work_dir,"library-analysis","normalised-guides-frequency.pdf"))
        plt.close()

        #
        datax2 = df['pre']
        X2 = datax2.to_numpy()
        X2 = X2 / pre_lib_sum

        datay2 = df['post']
        Y2 = datay2.to_numpy()
        Y2 = Y2 / pre_lib_sum

        data2 = X2 / Y2
        data2 = np.sort(data2)

        ax = sns.lineplot(x=range(index_len),y=data2, color='navy')
        ax.set(ylabel='Normalised \n pre-amplification/post-amplification', xlabel='sgRNA')
        ax.set_yscale('log')
        plt.tight_layout()
        plt.savefig(os.path.join(work_dir,
                    "library-analysis",
                    "normalised-pre-amplification-post-amplification.pdf"))
        plt.close()

        #Calculates Gini index of data sets
        
        ##Code taken and adapted from https://zhiyzuo.github.io/Plot-Lorenz/

        #function to calculate Gini index
        def gini(arr):
            ## first sort
            sorted_arr = arr.copy()
            sorted_arr.sort()
            n = arr.size
            coef_ = 2. / n
            const_ = (n + 1.) / n
            weighted_sum = sum([(i+1)*yi for i, yi in enumerate(sorted_arr)])
            return coef_*weighted_sum/(sorted_arr.sum()) - const_

        pre_gini_index = gini(X)
        pre_gini_index = round(pre_gini_index, 3)
        post_gini_index = gini(Y)
        post_gini_index = round(post_gini_index, 3)

        X_lorenz = X.cumsum() / X.sum()
        X_lorenz = np.insert(X_lorenz, 0, 0)
        X_lorenz[0], X_lorenz[-1]

        Y_lorenz = Y.cumsum() / Y.sum()
        Y_lorenz = np.insert(Y_lorenz, 0, 0)
        Y_lorenz[0], Y_lorenz[-1]

        #plots Lorenz curve
        fig, ax = plt.subplots(figsize=[6,6])
        ax.plot(np.arange(X_lorenz.size)/(X_lorenz.size-1),
                X_lorenz, label='Library pre-amplification',
                color='green')
        ax.plot(np.arange(Y_lorenz.size)/(Y_lorenz.size-1),
                Y_lorenz,
                label='Library post-amplification',
                color='red')
        ax.plot([0,1], [0,1], color='k', label='Ideal library')#line plot of equality
        ax.set(ylabel='Cumulative fraction of reads represented',
               xlabel='sgRNAs ranked by abundance')
        plt.text(0.075, 0.9, 'pre-amplification Gini index = '+str(pre_gini_index))
        plt.text(0.075, 0.85, 'post-amplification Gini index = '+str(post_gini_index))
        ax.legend(loc='lower right')
        plt.tight_layout()
        plt.savefig(os.path.join(work_dir,"library-analysis","lorenz-curve.pdf"))
        plt.close()
        '''
    else:
        return None


def gcBias(work_dir,library,crispr_library):
    #determine whether count file contains library samples pre and post
    header=subprocess.check_output(["head", "-1",
                                    os.path.join(work_dir,
                                                "count",
                                                "counts-aggregated.tsv")])
    if "pre" and "post" in str(header):

        out_file=os.path.join(work_dir,
                                "library-analysis",
                                "gc-bias.pdf")

        if not utils.file_exists(out_file):
            fasta=library[crispr_library]["fasta"]

            #get sgRNA counts
            df=pd.read_csv(os.path.join(work_dir,
                                    "count",
                                    "counts-aggregated.tsv"),sep='\t')

            df=df[["sgRNA","gene","pre","post"]]

            #get sgRNA sequences from fasta file
            df_fasta=pd.read_csv(fasta, header=None)
            df_fasta=df_fasta.rename(columns={0:"sgRNA"})
            df_seq=df_fasta[df_fasta["sgRNA"].str.contains(">")].reset_index(drop=True)
            df_seq["sequence"]=df_fasta[~df_fasta["sgRNA"].str.contains(">")].reset_index(drop=True)

            #get sgRNA sequence in sgRNA count df
            df["pre"]=df["pre"].astype(str)
            df["post"]=df["post"].astype(str)
            df_seq["sequence"]=df_seq["sequence"].astype(str)
            df_seq["sgRNA"]=df_seq["sgRNA"].astype(str)
            df_seq["sgRNA"]=df_seq["sgRNA"].str.replace(">","")
            df_seq["sgRNA"]=df_seq["sgRNA"].str.split("_",
                                                        n=1,
                                                        expand=True)[1]
            df = pd.merge(df, df_seq, on='sgRNA', how='left')
            #df = df.drop(columns="sequence_x")
            #df=df.rename(columns={"sequence_y":"sequence"})

            df["pre"]=df["pre"].astype(int)
            df["post"]=df["post"].astype(int)

            #calculate %GC for each sgRNA
            def calculateGC(x):
                total_N=len(x)
                total_G=x.count("G")
                total_C=x.count("C")
                GC_content=(total_G + total_C) / total_N *100
                return(GC_content)

            df["%GC"]=df["sequence"].apply(calculateGC)

            #select sgRNAs with 10% highest and 10% lowest abundances
            index_range=list(range(int(len(df)*0.1)))
            df_pre_low=df.sort_values(by=["pre"],
                                    ascending=True,
                                    inplace=False).reset_index(drop=True)
            df_pre_low=df_pre_low.iloc[index_range]

            df_post_low=df.sort_values(by=["post"],
                                    ascending=True,
                                    inplace=False).reset_index(drop=True)
            df_post_low=df_post_low.iloc[index_range]

            df_bottom=df_pre_low["%GC"]
            df_bottom=df_bottom.to_frame()
            df_bottom=df_bottom.rename(columns={df_bottom.columns[0]:"pre"})
            df_bottom["post"]=df_post_low["%GC"]


            df_pre_high=df.sort_values(by=["pre"],
                                    ascending=False,
                                    inplace=False).reset_index(drop=True)
            df_pre_high=df_pre_high.iloc[index_range]

            df_post_high=df.sort_values(by=["post"],
                                    ascending=False,
                                    inplace=False).reset_index(drop=True)
            df_post_high=df_post_high.iloc[index_range]

            df_top=df_pre_high["%GC"]
            df_top=df_top.to_frame()
            df_top=df_top.rename(columns={df_top.columns[0]:"pre"})
            df_top["post"]=df_post_high["%GC"]
            #####
            #df_all=df[["%GC","%GC"]]
            #df_all.columns.values[0]="pre"
            #df_all.columns.values[1]="post"

            def meltDf(df):
                df["id_var"]=range(len(df))
                df_melt=pd.melt(df,id_vars=["id_var"],value_vars=["pre","post"])
                df_melt=df_melt.rename(columns={"value":"%GC"})
                df_melt=df_melt.rename(columns={"variable":"sample"})
                return(df_melt)

            df_melt_top=meltDf(df_top)
            df_melt_top["group"]="top_10pc"
            df_melt_bottom=meltDf(df_bottom)
            df_melt_bottom["group"]="bottom_10pc"
            #df_melt_all=meltDf(df_all)
            #df_melt_all["group"]="all"

            df_melt=df_melt_bottom.append(df_melt_top,ignore_index=True)

            #plot data
            sns.set_style("white")
            sns.set_style("ticks")
            sns.despine()

            g=sns.displot(df_melt,x="%GC",
                            hue="sample",
                            fill=True,
                            col="group",
                            multiple="stack",
                            stat="density",
                            bins=15)

            for ax in g.axes.flat:
                ax.axvline(x=df["%GC"].median(),color='r',ls='--')

            plt.savefig(out_file)
            plt.close()


def essentialGenes(work_dir, analysis, essential_genes, fdr):
    '''
    Determines overlap with core essential genes from MAGeCK/BAGEL2gene lists,
    plots a Venn diagram and returns a list of genes that contains overlapping
    and non-overlapping genes.
    '''
    essential_genes = {line.strip() for line in open(essential_genes)}


    def plotVenn(test_genes_only,essential_genes_only,overlapping_genes,title,out_file):
        if len(overlapping_genes) > 0:
            venn=venn2(subsets=(len(test_genes_only),
                                        len(essential_genes_only),
                                        len(overlapping_genes)),
                               set_labels=("Dropouts","Core essential genes"),
                               set_colors=("darkorange","dodgerblue"),
                               alpha=0.7)
            venn.get_patch_by_id('11').set_edgecolor('black')
            venn.get_patch_by_id("11").set_linestyle("--")
            venn.get_patch_by_id("11").set_linewidth(1.25)
            venn.get_patch_by_id('10').set_edgecolor('black')
            venn.get_patch_by_id("10").set_linewidth(1)
            venn.get_patch_by_id('01').set_edgecolor('black')
            venn.get_patch_by_id("01").set_linewidth(1)
            plt.title(title)
            plt.savefig(out_file)
            plt.close()
        else:
            print("No overlapping genes found for "+title)
            return(None)


    if analysis == "mageck":
        print("Running essential gene analysis")
        #get list og MAGeCK gene summary file_exists
        file_list=glob.glob(os.path.join(work_dir,
                                "mageck*",
                                "*",
                                "*gene_summary.txt"))
        if len(file_list) == 0:
            print("ERROR: no MAGeCK output files found")
            return(None)

        for file in file_list:
            out_file=os.path.join(os.path.dirname(file),"essential_genes_venn.pdf")
            if not utils.file_exists(out_file):
                df=pd.read_csv(file,sep="\t")
                df_fdr=df.loc[df["neg|fdr"] < fdr ]
                if len(df_fdr) != 0:
                    test_genes=set(df_fdr["id"])

                    overlapping_genes=test_genes & essential_genes
                    essential_genes_only=essential_genes - test_genes
                    test_genes_only=test_genes - essential_genes

                    title=os.path.basename(file)
                    title=title.replace(".gene_summary.txt","")

                    plotVenn(test_genes_only,essential_genes_only,overlapping_genes,title,out_file)
    elif analysis == "bagel2":
        file_list=glob.glob(os.path.join(work_dir,
                                "bagel*",
                                "*",
                                "*.pr"))

        if len(file_list) == 0:
            print("ERROR: no BAGEL2 output files found")
            return(None)

        for file in file_list:
            #plot venn
            out_file=os.path.join(os.path.dirname(file),"essential_genes_venn.pdf")
            if not utils.file_exists(out_file):
                prefix=os.path.basename(os.path.normpath(file))
                prefix=prefix.replace(".pr","")
                df=pd.read_table(file)
                df_bf=df[(df["BF"] > 0)]
                test_genes=set(df_bf["Gene"])

                overlapping_genes=test_genes & essential_genes
                essential_genes_only=essential_genes - test_genes
                test_genes_only=test_genes - essential_genes

                title=os.path.basename(file)
                title=title.replace(".pr","")

                plotVenn(test_genes_only,essential_genes_only,overlapping_genes,title,out_file)


def goPython(work_dir,fdr,library,crispr_library,analysis,gene_sets):
    species=library[crispr_library]["species"]
    print(gene_sets)
    #check if chosen gene sets are valid
    default=["GO_Molecular_Function_2021",
            "GO_Cellular_Component_2021",
            "GO_Biological_Process_2021"]

    if not gene_sets == default:
        gene_sets=list(gene_sets.split(","))
        all_enrichr_gene_sets=gp.get_library_name()
        for i in gene_sets:
            if i not in all_enrichr_gene_sets:
                print("ERROR: invalid enrichR gene set chosen")
                print("List of valid gene sets: "+all_enrichr_gene_sets)
                return(None)

    if analysis == "mageck":
        #get list og MAGeCK gene summary file_exists
        file_list = glob.glob(os.path.join(work_dir,
                                "mageck*",
                                "*",
                                "*gene_summary.txt"))
        if len(file_list) == 0:
            print("ERROR: no MAGeCK output files found")
            return(None)

        input_list=[["neg|fdr","depletion"],["pos|fdr","enrichment"]]
        def enrichrMAGeCK(file_list,fdr,i):
            for file in file_list:
                prefix = os.path.basename(os.path.normpath(file))
                prefix = prefix.replace(".gene_summary.txt","")
                print("Performing gene set enrichment analysis with enrichR for "+prefix+": "+i[1])
                for set in gene_sets:
                    save_path = os.path.dirname(file)
                    df = pd.read_table(file)
                    df_fdr = df[(df[i[0]] < fdr)]
                    geneList = df_fdr["id"].to_list()
                    enrichr_results = gp.enrichr(gene_list=geneList,
                        gene_sets = set,
                        organism = species,
                        no_plot = False,
                        outdir = os.path.join(save_path,"enrichR",i[1]))

        for i in input_list:
            enrichrMAGeCK(file_list,fdr,i)

    elif analysis == "bagel2":
        file_list=glob.glob(os.path.join(work_dir,
                                "bagel*",
                                "*",
                                "*.pr"))

        if len(file_list) == 0:
            print("ERROR: no BAGEL2 output files found")
            return(None)

        for file in file_list:
            prefix=os.path.basename(os.path.normpath(file))
            prefix=prefix.replace(".pr","")
            print("Performing gene set enrichment analysis with enrichR for: "+prefix+"(depletion)")

            for set in gene_sets:
                save_path=os.path.dirname(file)
                df=pd.read_table(file)
                df_bf=df[(df["BF"] > 0)]
                geneList=df_bf["Gene"].to_list()
                enrichr_results=gp.enrichr(gene_list=geneList,
                    gene_sets=set,
                    organism=species,
                    no_plot=False,
                    outdir=os.path.join(save_path,"enrichR",i[1]))

