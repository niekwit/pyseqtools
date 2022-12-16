#!/usr/bin/env python3

"""
Created on Fri Dec 16 13:10:59 2022

@author: niek
"""

import pandas as pd
import sys

#get parsed arguments
args = sys.argv

macs3 = args[1]
homer = args[2]

#merge HOMER peak annotation with MACS3 output
out_file = homer.replace("_annotated-peaks.txt","_annotated_macs3_peaks.csv")
    
df_macs3 = pd.read_csv(macs3, skiprows=(29), sep = "\t")
df_homer = pd.read_csv(homer, sep = "\t")
df_homer = df_homer.rename(columns = {df_homer.columns[0]: "name"})
df_merge = pd.merge(df_macs3, df_homer[["name","Gene Name", "Gene Alias", 
                                        "Gene Description", "Gene Type"]], 
                    on = "name", 
                    how = "left" )
df_merge.to_csv(out_file, index = False)




