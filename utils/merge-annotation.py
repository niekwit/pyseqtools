#!/usr/bin/env python3

"""
Merges MACS3 xls file with HOMER peak annotations
"""

import pandas as pd
import sys

#get parsed arguments
args = sys.argv

macs3 = args[1]
homer = args[2]

#load data
df_macs3 = pd.read_csv(macs3, skiprows=29,sep = "\t")

df_homer = pd.read_csv(homer, sep = "\t")
df_homer = df_homer.rename(columns = {df_homer.columns[0]: "name"})
df_homer = df_homer.drop(["Chr","Start","End","Strand","Peak Score","Focus Ratio/Region Size"],axis=1)

#merge HOMER peak annotation with MACS3 output
out_file = homer.replace("_annotated-peaks.txt","_annotated_macs3_peaks.csv")
df_merge = pd.merge(df_macs3, df_homer, on = "name", how = "left" )

#save to csv file
df_merge.to_csv(out_file, index = False)




