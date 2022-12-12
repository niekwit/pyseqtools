#!/usr/bin/env python3

"""
Counts reads from pre and post deduplication bam files with pysamtools and saves to csv
"""
import sys
import glob
import os
import pandas as pd

import pysam

args = sys.argv

#get parsed arguments
work_dir = args[1]
script_dir = args[2]
genome = args[3]
threads = args[4]

#get bam files
non_dedup_bam = glob.glob(os.path.join(work_dir,"bam",genome,"*-sort-bl.bam"))

#create df to store read counts
sample_names = [os.path.basename(bam).replace("-sort-bl.bam","") for bam in non_dedup_bam]
df = pd.DataFrame(columns = ["sample","pre-dedup_count","post-dedup_count"])
df["sample"] = sample_names

#count reads
for bam in non_dedup_bam:
    dedup_bam = bam.replace("-sort-bl.bam","-sort-bl-dedupl.bam") 
    sample = os.path.basename(bam).replace("-sort-bl.bam","")
    
    #count reads
    count_non_dedup = pysam.view("-@", str(threads) ,"-c", "-F" "260", bam)
    count_non_dedup = int(count_non_dedup.replace("\n",""))
    
    count_dedup = pysam.view("-@", str(threads) ,"-c", "-F" "260", dedup_bam)
    count_dedup = int(count_dedup.replace("\n",""))
    
    #add counts to df
    index = df[df["sample"] == sample].index[0]
    df.loc[index,"pre-dedup_count"] = count_non_dedup
    df.loc[index,"post-dedup_count"] = count_dedup
    
#save df to csv for plotting with R
csv = os.path.join(work_dir, "bam", "read-counts.csv")
df.to_csv(index = False)




