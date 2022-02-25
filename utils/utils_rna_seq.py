#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import glob
import os
import subprocess
import yaml
import sys
import pkg_resources
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import gseapy as gp
from gseapy.plot import gseaplot
import math
import urllib.request

script_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.join(script_dir, "utils"))
import utils_general as utils

def install_packages(): #check for required python packages; installs if absent
    required = {"pyyaml, cutadapt, multiqc"}
    installed = {pkg.key for pkg in pkg_resources.working_set}
    missing = required - installed
    if missing:
        python = sys.executable
        print("Installing missing required Python3 packages")
        subprocess.check_call([python, '-m', 'pip3', 'install', 
                               *missing], stdout = subprocess.DEVNULL)


def salmon(salmon_index, threads, work_dir, gtf, fasta, script_dir, settings, reference):
    
    #check for salmon
    path = os.environ["PATH"].lower()
    
    if "salmon" not in path:
        #Check for Salmon elsewhere
        salmon = [line[0:] for line in subprocess.check_output("find $HOME -name salmon ! -path '*/multiqc*'", 
                                                               shell = True).splitlines()]
        
        salmon_file = None
        for i in salmon:
            i = i.decode("utf-8")
            if "bin/salmon" in i:
                salmon_file = i
            
        if salmon_file == None:
            print("WARNING: Salmon was not found\nInstalling Salmon now")
            url = "https://github.com/COMBINE-lab/salmon/releases/download/v1.5.2/salmon-1.5.2_linux_x86_64.tar.gz"
            download_file = os.path.join(script_dir,"salmon-1.5.2_linux_x86_64.tar.gz")
            urllib.request.urlretrieve(url, download_file)
            #untar Salmon download file
            tar_command = "tar -xzf " + download_file + " --directory " + script_dir
            subprocess.run(tar_command,
                           shell = True)
            
            #remove download file
            os.remove(download_file)
            
            salmon_file = os.path.join(script_dir, 
                                       "salmon-1.5.2_linux_x86_64", 
                                       "bin", 
                                       "salmon")
    else:
        salmon_file = "salmon"
    
    
    if salmon_index == "": #Salmon index not found, make it
        print("No Salmon index found: generating Salmon index")
        if os.path.isfile(fasta):
            index_dir = os.path.join(script_dir,"index", "salmon", reference)
            os.makedirs(index_dir, 
                        exist_ok = True)
            salmon_index_command = [salmon_file, "index", "-t", fasta, "-i", index_dir, "--gencode"]
            #log commands
            with open(os.path.join(work_dir, "commands.log"), "a") as file:
                file.write("Salmon index: ")
                print(*salmon_index_command, 
                      sep = " ", 
                      file = file)

            subprocess.run(salmon_index_command) #run Salmon index

            #Write salmon index file location to rna-seq.yaml
            with open(os.path.join(script_dir, "yaml", "rna-seq.yaml")) as f:
                doc = yaml.safe_load(f)
            doc["salmon_index"][reference] = index_dir
            with open(os.path.join(script_dir,"yaml" ,"rna-seq.yaml"), "w") as f:
                yaml.dump(doc,f)
        else:
            print("ERROR: no FASTA file specified in rna-seq.yaml")
            sys.exit()

    
    
    
    def salmonSE(work_dir, threads, gtf, salmon_index, reference):
        salmon_output_dir = os.path.join(work_dir,"salmon")
        os.makedirs(salmon_output_dir, 
                exist_ok = True)
        trim_list = glob.glob(os.path.join(work_dir,"trim/*_trimmed.fq.gz"))
        for read1 in trim_list:
            base_read1 = os.path.basename(read1).replace("_trimmed.fq.gz", "") + "-quant"
            salmon_folder_test = os.path.join(salmon_output_dir, base_read1)
            if not utils.file_exists(salmon_folder_test):
                print("Mapping sample " + read1.replace("_trimmed.fq.gz", ""))
                
                out_file = os.path.basename(read1.replace("_trimmed.fq.gz",""))
                salmon_output_file = os.path.join(work_dir, "salmon", out_file)+"-quant"
                salmon_index = settings["salmon_index"][reference] #reload index
                salmon_command = [salmon_file, "quant", "--index", salmon_index, "-l", "A",
                "-g", gtf, "-p", threads, "-r", read1, "--validateMappings",
                "--gcBias", "-o", salmon_output_file]
                with open(os.path.join(work_dir,"commands.log"), "a") as file:
                    file.write("Salmon quant: ")
                    print(*salmon_command, sep = " ", file = file)
                subprocess.run(salmon_command) #Run Salmon quant
    
    def salmonPE(work_dir, threads, gtf, salmon_index, reference):
        salmon_output_dir = os.path.join(work_dir,"salmon")
        os.makedirs(salmon_output_dir, 
                exist_ok = True)
        trim_list = glob.glob(os.path.join(work_dir,"trim/*_R1_001_val_1.fq.gz"))
        for read1 in trim_list:
            base_read1 = os.path.basename(read1).replace("_R1_001_val_1.fq.gz", "") + "-quant"
            salmon_folder_test = os.path.join(salmon_output_dir, base_read1)
            if not utils.file_exists(salmon_folder_test):
                print("Mapping sample " + read1.replace("_R1_001_val_1.fq.gz", ""))
                read2 = read1.replace("R1_001_val_1.fq.gz", "R2_001_val_2.fq.gz")
                out_file = os.path.basename(read1.replace("_R1_001_val_1.fq.gz",""))
                salmon_output_file = os.path.join(work_dir,"salmon",out_file)+"-quant"
                salmon_index = settings["salmon_index"][reference] #reload index
                salmon_command = [salmon_file,"quant","--index",salmon_index,"-l","A",
                "-g", gtf,"-p",threads,"-1", read1,"-2",read2,"--validateMappings",
                "--gcBias","-o", salmon_output_file]
                with open(os.path.join(work_dir,"commands.log"), "a") as file:
                    file.write("Salmon quant: ")
                    print(*salmon_command, sep = " ", file = file)
                subprocess.run(salmon_command) #Run Salmon quant
            
    
    if utils.getEND(work_dir) == "PE":
        print("Mapping reads with Salmon (paired-end):")
        salmonPE(work_dir, threads, gtf, salmon_index, reference)
    elif utils.getEND(work_dir) == "SE":
        print("Mapping reads with Salmon (single-end):")
        salmonSE(work_dir, threads, gtf, salmon_index, reference)
  
            
def plotBar(df,y_label,save_file):
    sns.set_style("white")
    sns.set_style("ticks")
    sns.barplot(x=list(df.keys())[0],
                    y=list(df.keys())[1],
                    data=df,
                    color="cornflowerblue",
                    edgecolor="black",
                    linewidth=1)
    plt.ylabel(y_label)
    plt.xticks(rotation = 'vertical')
    plt.xlabel("")
    plt.ylim(0, 100)
    plt.tight_layout()
    sns.despine()
    plt.savefig(save_file)
    plt.close()


def plotMappingRate(work_dir):
    file_list = glob.glob(os.path.join(work_dir, "salmon", "*", "logs", "salmon_quant.log"))
    save_file = os.path.join(work_dir,"salmon","mapping_rates.pdf")
    mapping_rate = []
    samples = []
    df = pd.DataFrame(columns = ["sample","Mapping rate (%)"],
                      index = np.arange(len(file_list)))

    if not utils.file_exists(save_file): 
        for file in file_list:
                sample = os.path.dirname(file)
                sample = sample.replace(os.path.join(work_dir, "salmon"),"")
                sample = sample.replace("/log", "")
                sample = sample.replace("-quants", "")
                sample = sample.replace("/", "")
                samples.append(sample)
                with open(file,"r") as file:
                    for line in file:
                        if "[info] Mapping rate" in line:
                            rate=line.rsplit(" ", 1)[1]
                            rate=rate.replace("%", "")
                            mapping_rate.append(rate)
        
        df["sample"]=samples
        df["Mapping rate (%)"] = mapping_rate
        df["Mapping rate (%)"] = pd.to_numeric(df["Mapping rate (%)"])
        df = df.sort_values(by = ["sample"],
                          ascending = True,
                          inplace = False).reset_index(drop = True)
        
        plotBar(df, "Mapping rate (%)", save_file)


def plotVolcano(work_dir):
    file_list = glob.glob(os.path.join(work_dir,
                                     "DESeq2",
                                     "*",
                                     "DESeq-output.csv"))
    
    for file in file_list:
        base_name=os.path.basename(os.path.dirname(file))
        
        out_dir=os.path.dirname(file)
        out_file=os.path.join(out_dir,base_name+"-volcano.pdf")
        if not utils.file_exists(out_file):
            print("Generating volcano plot for: "+base_name)
            df=pd.read_csv(file)
            df["log.p.value"]=-np.log(df["padj"])
            df["label"] = ""
            
            #mark genes are upregulated or downregulated for plotting
            conditions=[(df["log2FoldChange"] > 0.5) & (df["log.p.value"] > 3),
                        (df["label"] != "up") & (df["log2FoldChange"] < -0.5) & (df["log.p.value"] > 3),
                        (df["label"] != "up") & (df["label"] != "down")]
            choices=["Upregulated","Downregulated","None"]
            df["label"]=np.select(conditions, 
                                  choices, 
                                  default=np.nan)
            
            #plot
            sns.set_style("white")
            sns.set_style("ticks")
            sns.scatterplot(data=df, 
                            x="log2FoldChange", 
                            y="log.p.value",
                            alpha=0.5,
                            linewidth=0.25,
                            edgecolor="black",
                            hue="label",
                            palette=["red", "blue", "black"],
                            legend= False)
            plt.xlabel("log2(FC)")
            plt.ylabel("-log(adjusted P value)")
            sns.despine()
            plt.savefig(out_file)
            plt.close()


def plotPCA(work_dir,script_dir):
    out_file = os.path.join(work_dir,"salmon","PCAplot.pdf")
    if not utils.file_exists(out_file):
        PCA_command = ["Rscript", 
                       os.path.join(script_dir, 
                                    "R", 
                                    "rna-seq-pcaplot.R"), 
                       work_dir]
        with open(os.path.join(work_dir,"commands.log"), "a") as file:
            file.write("PCA plot (all samples): ")
            print(*PCA_command, 
                  sep = " ", 
                  file = file)
    
    try:
        subprocess.run(PCA_command)
    except:
        print("PCA plot for all samples failed, check log")
        return(None)
  
    
def hisat2():
    print("Mapping reads with HISAT2 (UNDER CONSTRUCTION)")
    sys.exit()


def diff_expr(work_dir,gtf,script_dir,species,pvalue):
    samples_input=os.path.join(work_dir,"samples.csv")
    if not os.path.exists(samples_input):
        sys.exit("ERROR: "+samples_input+" not found")
    
    
    deseq2_command=["Rscript",os.path.join(script_dir, 
                                           "R", 
                                           "rna-seq_deseq2.R"),
                    work_dir,
                    gtf,
                    script_dir,
                    species,
                    str(pvalue)]
    
    with open(os.path.join(work_dir,"commands.log"), "a") as file:
        file.write("DESeq2: ")
        print(*deseq2_command, 
              sep = " ", 
              file = file)
    print("Running differential expression analysis with DESeq2")
    os.makedirs(os.path.join(work_dir, "DESeq2"),
                exist_ok = True)
    subprocess.run(deseq2_command)


def geneSetEnrichment(work_dir,pvalue,gene_sets):
    file_list = glob.glob(os.path.join(work_dir,
                                     "DESeq2",
                                     "*",
                                     "DESeq-output.csv"))
    
    if len(file_list) == 0:
            print("ERROR: no DESeq2 output files found")
            return(None)
    
    
    def doEnrichr(gene_list, gene_sets, out_dir, output_name):
        out_dir = os.path.join(out_dir,"Enrichr",output_name)
        
        enrichr_results = gp.enrichr(gene_list = gene_list,
                                   gene_sets = gene_sets,
                                   outdir = out_dir)
    
   
    def GSEA(df,gene_set, out_dir):#not working properly yet
        rnk = df.dropna()
        rnk = rnk[["SYMBOL","log2FoldChange"]]
        rnk = rnk.drop_duplicates(subset="SYMBOL")
        pre_res = gp.prerank(rnk = rnk,
                            gene_sets = gene_set,
                            processes = 4,
                            permutation_num = 100,
                            format = "pdf",
                            seed = 6,
                            no_plot = True)
        terms=pre_res.res2d.index
        df_out=pre_res.res2d
        df_out.reset_index(inplace = True)
        os.makedirs(os.path.join(out_dir,"GSEA"), 
                    exist_ok = True)
        df_out.to_csv(os.path.join(out_dir,"GSEA",'GSEA.csv'), 
                      index = False)
        
        
        
        for i in range(10): #plot GSEA plots for top 10 terms
            GSEA_dir=os.path.join(out_dir,"GSEA",gene_set) 
            os.makedirs(GSEA_dir,
                        exist_ok=True)
            try:
                gseaplot(rank_metric=pre_res.ranking,
                         term=terms[i],
                         pheno_pos="Upregulated genes",
                         pheno_neg="Downregulated genes",
                         ofname=os.path.join(GSEA_dir,terms[i]+".pdf"),
                         **pre_res.results[terms[i]])
            except FileNotFoundError:
                continue
    
    logPvalue = -math.log10(pvalue)
    
    for file in file_list:
       df = pd.read_csv(file)
       df["log.p.value"] = -np.log10(df["padj"])
       out_dir = os.path.dirname(file)
              
       df_up = df[(df["log2FoldChange"] > 0.5) & (df["log.p.value"] > logPvalue)]
       upregulated_genes = list(df_up["SYMBOL"])
       
       df_down = df[(df["log2FoldChange"] < -0.5) & (df["log.p.value"] > logPvalue)]
       downregulated_genes = list(df_down["SYMBOL"])
       
       input_list = [upregulated_genes, 
                     downregulated_genes]
       output_names = ["upregulated_genes",
                       "downregulated_genes"]
       
       #run Enrichr
       for x,y in zip(input_list,output_names):
           doEnrichr(x,gene_sets,out_dir,y)
           
       #Run GSEA
       for i in gene_sets:
          GSEA(df,i,out_dir)    
    
