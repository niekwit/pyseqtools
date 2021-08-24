#!/usr/bin/env python3

import glob
import os
import subprocess
import urllib.request
import sys
import gzip
import shutil

from tqdm.auto import tqdm
import yaml

script_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.join(script_dir, "utils"))
import utils_general as utils


###CHIP-SEQ ANALYSIS SPECIFIC FUNCTIONS
 

def hisat2(script_dir, work_dir, threads, chip_seq_settings, genome):
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
                                 "HISAT2",
                                 genome,
                                 "index")
        os.makedirs(os.path.join(script_dir,
                                 "index",
                                 "HISAT2",
                                 genome), 
                    exist_ok = True)
        build_command = "python3 " + hisat2 + "-build " + ucsc_fasta + " " + index_location
        subprocess.run(build_command,
                           shell = True)
        
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
                elif genome == "mm9":
                    url = "https://genome-idx.s3.amazonaws.com/hisat/mm10_genome.tar.gz"
                
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
                                             genome),
                                exist_ok = True)
                    out_put_list = [os.path.join(script_dir, 
                                                 "fasta", 
                                                 genome, 
                                                 os.path.basename(i)) for i in download_file_list]
                    #download fasta files
                    print("Downloading fasta files from UCSC needed for building HISAT2 index")
                  
                    for i,j in tqdm(zip(download_file_list, out_put_list),position = 0, leave = True):
                        urllib.request.urlretrieve(i, j)
                        
                    #concatenate fasta files to build genome fasta
                    print("Building whole genome fasta file")
                    ucsc_fasta = os.path.join(script_dir, "fasta", genome, "ucsc." + genome + ".fasta")
                    zcat_command = "zcat " + " ".join(out_put_list) + " > " + ucsc_fasta
                    subprocess.run(zcat_command,
                                   shell = True)
                    
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
    
    ###look for blacklist###
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
        elif genome == "hg38":
            url = "https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz"
            blacklist = getBlacklist(script_dir, url, genome)
        elif genome == "mm9":
            url = "http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/mm9-mouse/mm9-blacklist.bed.gz"
            blacklist = getBlacklist(script_dir, url, genome)
    elif os.path.isfile(blacklist):
        if blacklist.endswith(".bed"):
            if os.path.getsize(blacklist) > 0:
                pass
    
    ###perform alignment with HISAT2###
    os.makedirs(os.path.join(work_dir, "bam"),
                    exist_ok = True)
    
    def alignSE(work_dir, hisat2, blacklist, index_location, threads):
        file_list = glob.glob(os.path.join(work_dir, "trim_galore","*trimmed.fq.gz"))
        
        print("Generating BAM files with HISAT2 (single-end mode)")
        
        for file in file_list:
            hisat2_output = os.path.basename(file).replace("_trimmed.fq.gz",
                                                           "-sort-bl.bam")
            hisat2_output = os.path.join(work_dir,
                                         "bam",
                                         hisat2_output)
            
            if not utils.file_exists(hisat2_output):
                samtools = utils.checkSamtools(script_dir)
                bedtools = utils.checkBedtools(script_dir)
                
                
                align_command = "zcat " + file + " | " + hisat2 + " -p " + str(threads) + " -x " + index_location + " - 2>> align.log | " + samtools + " view -q 15 -F 260 -bS -@ " + str(threads) + " - | " + bedtools + " intersect -v -a 'stdin' -b " + blacklist + " -nonamecheck | " + samtools + " sort -@ " + str(threads) + " - > " + hisat2_output
                
                
                utils.write2log(work_dir, align_command, "HISAT2: ")
                
                subprocess.run(align_command,
                           shell = True)


    def alignPE(work_dir, hisat2, blacklist, index_location, threads):
        file_list = glob.glob(os.path.join(work_dir, "trim_galore","*1_val_1.fq.gz"))
        
        print("Generating BAM files with HISAT2 (paired-end mode)")
        
        for read1 in file_list:
            read2 = read1.replace("R1_001_val_1.fq.gz", 
                                  "R2_001_val_2.fq.gz")
            
            
            hisat2_output = os.path.basename(read1).replace("_R1_001_val_1.fq.gz",
                                                           "-sort-bl.bam")
            hisat2_output = os.path.join(work_dir,
                                         "bam",
                                         hisat2_output)
            
            if not utils.file_exists(hisat2_output):
                samtools = utils.checkSamtools(script_dir)
                bedtools = utils.checkBedtools(script_dir)
                                
                align_command = hisat2 + " -p " + str(threads) + " -x " + index_location
                align_command = align_command + " -1 " + read1 + " -2 " + read2
                align_command = align_command + " 2>> align.log | " + samtools + " view -q 15 -F 260 -bS -@"
                align_command = align_command + threads + " - | " + bedtools + " intersect -v -a 'stdin' -b "
                align_command = align_command + blacklist + " -nonamecheck | " + samtools + " sort -@ "
                align_command = align_command + threads + " - > " + hisat2_output
                
                utils.write2log(work_dir, align_command, "HISAT2 PE: ")
                
                subprocess.run(align_command,
                           shell = True)
                
           
    if utils.getEND(work_dir) == "PE":
        alignPE(work_dir, hisat2, blacklist, index_location, threads)
    elif utils.getEND(work_dir) == "SE":
        alignSE(work_dir, hisat2, blacklist, index_location, threads)
    
    #plot alignment rate
    
    
    #plot number of reads before deduplication
    
    
def bwa():
    pass


def ngsplot(work_dir, genome, feature, window):
    
    file_list = glob.glob(os.path.join(work_dir, "bam", "*.bam"))
    
    os.makedirs(os.path.join(work_dir, "ngsplot"),
                exist_ok = True)
    
    
    def ngsplotFunction(work_dir, genome, feature, window, extension):
        base_name = os.path.basename(bam).replace(extension, "")
        ngsplot_dir = os.path.join(os.path.dirname(bam).replace("bam",""), "ngsplot" , base_name)
        os.makedirs(ngsplot_dir,
            exist_ok = True)
        ngsplot_output = os.path.join(ngsplot_dir, base_name) + "_" + feature
        
        ngsplot = "ngs.plot.r -G " + genome + " -R " + feature + " -C " + bam + " -O " + ngsplot_output + "_"+ feature + " -T " + base_name + " -L " + window

        subprocess.run(ngsplot,
                       shell = True)
        
    
    for bam in file_list:
        if "-sort-bl-dedupl.bam" in bam:
            ngsplotFunction(work_dir, genome, feature, window, "-sort-bl-dedupl.bam")
            
        elif "-sort-bl.bam" in bam:
            ngsplotFunction(work_dir, genome, feature, window, "-sort-bl.bam")
     
    
    
    url= "https://drive.google.com/file/d/0B5ldivL0Hd2JN05MOEFuZ0FRQTA/view?usp=sharing&resourcekey=0-Y6Sq22xOYTAb9Yng8ZlmJg"
    
    
    
    
    
    