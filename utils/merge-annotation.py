#!/usr/bin/env python3

"""
Merges MACS3 xls file with HOMER peak annotations
"""

import pandas as pd
import sys
import subprocess

#get parsed arguments
args = sys.argv

macs3 = args[1]
homer = args[2]

#remove all comments and empty rows from peak data
macs3_clean = macs3.replace(".xls","_cleaned.xls")
bash = f"grep -v '#' {macs3} | sed '/^$/d' > {macs3_clean}"
subprocess.run(bash,shell=True)

#load data
df_macs3 = pd.read_csv(macs3_clean,sep = "\t")

df_homer = pd.read_csv(homer, sep = "\t")
df_homer = df_homer.rename(columns = {df_homer.columns[0]: "name"})
df_homer = df_homer.drop(["Chr","Start","End","Strand","Peak Score","Focus Ratio/Region Size"],axis=1)

#merge HOMER peak annotation with MACS3 output
out_file = homer.replace("_annotated-peaks.txt","_annotated_macs3_peaks.csv")
df_merge = pd.merge(df_macs3, df_homer, on = "name", how = "left" )

#save to csv file
df_merge.to_csv(out_file, index = False)




