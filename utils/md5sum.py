#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Checks md5sums of fastqc files
"""
import sys
import os
import hashlib
import glob
import pandas as pd

args = sys.argv

#get parsed arguments
work_dir = args[0]

#get md5sum files
md5sum_files = glob.glob(os.path.join(work_dir,"raw-data","*.md5sums.txt"))

#df to store md5sums
df = pd.DataFrame(columns=["fastq","md5sum","md5sum_calculated","correct"],
                  index = range(0,len(md5sum_files)*2))

#function to cacluate md5sum
def md5(work_dir,file):
    file = os.path.join(work_dir, "raw-data", file)
    hash_md5 = hashlib.md5()
    with open(file, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return(hash_md5.hexdigest())

#for index,row in df.iterrows():
for i in range(0,len(md5sum_files)):
    #read file
    file = open(os.path.join(work_dir,"raw-data",md5sum_files[i]), "r")
    lines = file.readlines()
        
    for l in range(0,len(lines)):
        #add all fastq files and their md5sums to df
        md5sum = lines[l].split("  ")[0]
        df.at[i*2+l,"md5sum"] = md5sum
        
        fastq = lines[l].split("  ")[1].replace("\n","")
        df.at[i*2+l,"fastq"] = fastq

        #calculate md5sums of actual fastq files and add to df
        md5sum_calc = md5(work_dir,fastq)
        df.at[i*2+l,"md5sum_calculated"] = md5sum_calc
        
#compare original checksums with calculated ones
df["correct"] = df["md5sum"] == df["md5sum_calculated"]

#save df to csv
df.to_csv(os.path.join(work_dir,"raw-data","md5sums_checked.csv"),index=False)      
     
#check for different check sums and save to separate csv
df_fail = df[df["correct"] == False]

if len(df_fail.axes[0]) > 0:
    df_fail.to_csv(os.path.join(work_dir,"raw-data","md5sums_failed.csv"),index=False)
    
